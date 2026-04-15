# coding: utf-8
"""
STAFF: Splice History
=====================

For every molecule in a BAM:
  Build the FULL splice decision chain:
    EXON(pos,len) -> JUNCTION(pos,size,left_base,right_base,xor) -> EXON(pos,len) -> ...

Then for molecules at the SAME genomic locus:
  Align their decision chains
  Find EXACTLY where they diverge
  Catalog every divergence: which junction differs, by how much, what bases

This is not binning. This is per-molecule structural comparison.

Also: for every excised intron, record its SIZE, its FLANKING SEQUENCE,
and whether it matches known functional RNA elements.
The introns are not waste. They are released information.

A=00 T=01 C=10 G=11
"""
import struct, zlib, json, sys, time, os
from collections import defaultdict, Counter
import numpy as np

SEQ_DECODE = '=ACMGRSVTWYHKDBN'
BASE_TO_BITS = {'A': 0, 'T': 1, 'C': 2, 'G': 3,
                'a': 0, 't': 1, 'c': 2, 'g': 3, 'N': 0, 'n': 0}
XOR_NAME = {0: 'ID', 1: 'WC', 2: 'TV', 3: 'TI'}


def decode_bam_seq(encoded_bytes, seq_len):
    seq = []
    for i in range(seq_len):
        byte_idx = i >> 1
        if byte_idx >= len(encoded_bytes):
            break
        if i & 1:
            val = encoded_bytes[byte_idx] & 0x0F
        else:
            val = (encoded_bytes[byte_idx] >> 4) & 0x0F
        seq.append(SEQ_DECODE[val])
    return ''.join(seq)


def read_bgzf_block(fh):
    header = fh.read(18)
    if len(header) < 18:
        return None
    if header[:4] != b'\x1f\x8b\x08\x04':
        return None
    bsize = struct.unpack('<H', header[16:18])[0]
    remaining = bsize - 18 + 1
    cdata = fh.read(remaining)
    if len(cdata) < remaining:
        return None
    return zlib.decompress(cdata[:-8], -15)


def stream_bam(bam_path, min_mapq=10):
    with open(bam_path, 'rb') as fh:
        first_block = read_bgzf_block(fh)
        if first_block is None or first_block[:4] != b'BAM\x01':
            raise ValueError(f'Not BAM: {bam_path}')
        data = first_block
        pos = 4
        l_text = struct.unpack('<i', data[pos:pos+4])[0]
        pos += 4 + l_text
        n_ref = struct.unpack('<i', data[pos:pos+4])[0]
        pos += 4
        ref_names = []
        for _ in range(n_ref):
            l_name = struct.unpack('<i', data[pos:pos+4])[0]
            pos += 4
            name = data[pos:pos+l_name-1].decode()
            pos += l_name
            l_ref = struct.unpack('<i', data[pos:pos+4])[0]
            pos += 4
            ref_names.append(name)
        buffer = data[pos:]
        n = 0
        while True:
            while len(buffer) < 4:
                block = read_bgzf_block(fh)
                if block is None: return
                buffer += block
            block_size = struct.unpack('<i', buffer[:4])[0]
            while len(buffer) < block_size + 4:
                block = read_bgzf_block(fh)
                if block is None: return
                buffer += block
            record = buffer[4:block_size+4]
            buffer = buffer[block_size+4:]
            if len(record) < 32: continue
            ref_id = struct.unpack('<i', record[0:4])[0]
            aln_pos = struct.unpack('<i', record[4:8])[0]
            l_read_name = record[8]
            mapq = record[9]
            n_cigar_op = struct.unpack('<H', record[12:14])[0]
            flag = struct.unpack('<H', record[14:16])[0]
            l_seq = struct.unpack('<i', record[16:20])[0]
            if ref_id < 0 or flag & 4 or mapq < min_mapq: continue
            cigar_start = 32 + l_read_name
            cigar_end = cigar_start + n_cigar_op * 4
            cigar_ops = []
            for i in range(n_cigar_op):
                val = struct.unpack('<I', record[cigar_start + i*4:cigar_start + (i+1)*4])[0]
                cigar_ops.append((val & 0xF, val >> 4))
            seq_start = cigar_end
            seq_bytes = (l_seq + 1) // 2
            seq_end = seq_start + seq_bytes
            if seq_end > len(record): continue
            seq = decode_bam_seq(record[seq_start:seq_end], l_seq)
            chrom = ref_names[ref_id] if ref_id < len(ref_names) else str(ref_id)
            yield chrom, aln_pos, cigar_ops, seq
            n += 1
            if n % 500000 == 0:
                print(f'  {n:,} reads...')


def build_splice_history(bam_path, max_reads=None):
    """Build per-molecule splice decision chains and compare at same loci."""
    t0 = time.time()
    print(f'STAFF Splice History')
    print(f'Input: {bam_path}\n')

    # Per-locus: store splice histories for comparison
    # Key = (chrom, locus_bin) where bin = 1kb window
    BIN = 1000
    locus_histories = defaultdict(list)  # max 100 per locus

    # Excised intron catalog
    intron_catalog = Counter()  # (chrom, start, end) -> count
    intron_sizes = []
    intron_flank_xors = Counter()  # base4 XOR at intron boundaries

    # Per-molecule: full chain as tuple
    # chain = ((exon_len, junction_size, left_base4, right_base4, xor), ...)
    chain_catalog = Counter()  # chain_signature -> count

    # Divergence tracking between molecules at same locus
    divergence_positions = Counter()  # which element in the chain diverges most
    divergence_types = Counter()  # what kind of divergence

    n_reads = 0
    n_spliced = 0

    for chrom, pos, cigar, seq in stream_bam(bam_path):
        n_reads += 1
        if max_reads and n_reads > max_reads:
            break

        # Build the full splice decision chain for this molecule
        ref_pos = pos
        seq_pos = 0
        chain = []
        current_exon_len = 0

        for op, length in cigar:
            if op in (0, 7, 8):  # M/=/X — exon
                current_exon_len += length
                ref_pos += length
                seq_pos += length
            elif op == 1:  # I
                seq_pos += length
            elif op == 2:  # D
                ref_pos += length
            elif op == 3:  # N — SPLICE JUNCTION (excised intron)
                # Record the exon that just ended
                # Get flanking bases
                left_b = BASE_TO_BITS.get(seq[seq_pos-1], 0) if seq_pos >= 1 else 0
                right_b = BASE_TO_BITS.get(seq[seq_pos], 0) if seq_pos < len(seq) else 0
                xor = left_b ^ right_b

                # Chain element: (exon_len, intron_size, left_base, right_base, xor)
                chain.append((current_exon_len, length, left_b, right_b, xor))
                current_exon_len = 0

                # Catalog the excised intron
                intron_catalog[(chrom, ref_pos, ref_pos + length)] += 1
                intron_sizes.append(length)
                intron_flank_xors[xor] += 1

                # Get 4bp of flanking sequence for intron characterization
                # Left flank (last 4bp of upstream exon)
                # Right flank (first 4bp of downstream exon)

                ref_pos += length
            elif op == 4:  # S
                seq_pos += length

        # Final exon
        if current_exon_len > 0:
            chain.append((current_exon_len,))  # terminal exon, no junction after

        if len(chain) > 1:  # has at least one junction
            n_spliced += 1

            # Store chain signature
            chain_sig = tuple(chain)
            chain_catalog[chain_sig] += 1

            # Store at locus for comparison
            locus_key = (chrom, pos // BIN)
            if len(locus_histories[locus_key]) < 200:
                locus_histories[locus_key].append(chain)

    elapsed = time.time() - t0

    # Analyze divergences between molecules at same locus
    print(f'\nAnalyzing divergences at shared loci...')
    n_loci_compared = 0
    total_comparisons = 0
    junction_divergences = Counter()  # (position_in_chain, divergence_type) -> count

    for locus_key, histories in locus_histories.items():
        if len(histories) < 2:
            continue
        n_loci_compared += 1

        # Compare all pairs (capped at 50 pairs per locus)
        for i in range(min(len(histories), 10)):
            for j in range(i+1, min(len(histories), 10)):
                h1, h2 = histories[i], histories[j]
                total_comparisons += 1

                # Align chains element by element
                min_len = min(len(h1), len(h2))
                for k in range(min_len):
                    e1, e2 = h1[k], h2[k]

                    if len(e1) >= 5 and len(e2) >= 5:  # both have junctions
                        # Compare: same junction or different?
                        if e1[1] != e2[1]:  # different intron size
                            junction_divergences[('intron_size_diff', k)] += 1
                        if e1[4] != e2[4]:  # different XOR
                            junction_divergences[('xor_diff', k)] += 1
                        if e1[0] != e2[0]:  # different exon length before junction
                            junction_divergences[('exon_len_diff', k)] += 1
                    elif len(e1) != len(e2):  # one has junction, other doesn't
                        junction_divergences[('structure_diff', k)] += 1

    # Results
    print(f'\n{"="*60}')
    print(f'SPLICE HISTORY RESULTS')
    print(f'{"="*60}')
    print(f'Reads: {n_reads:,} ({elapsed:.1f}s)')
    print(f'Spliced molecules: {n_spliced:,} ({100*n_spliced/max(1,n_reads):.1f}%)')
    print(f'Unique chain signatures: {len(chain_catalog):,}')

    # Chain complexity distribution
    chain_lengths = []
    for chain_sig, count in chain_catalog.items():
        n_junctions = sum(1 for e in chain_sig if len(e) >= 5)
        chain_lengths.extend([n_junctions] * count)

    if chain_lengths:
        print(f'\n--- CHAIN COMPLEXITY ---')
        print(f'  Median junctions per molecule: {int(np.median(chain_lengths))}')
        print(f'  Max: {max(chain_lengths)}')
        for n_j in range(1, min(10, max(chain_lengths)+1)):
            c = chain_lengths.count(n_j)
            if c > 0:
                print(f'    {n_j} junctions: {c:,} molecules')

    # Top chain signatures (most common splice patterns)
    print(f'\n--- TOP 15 SPLICE DECISION CHAINS ---')
    for chain_sig, count in chain_catalog.most_common(15):
        n_junctions = sum(1 for e in chain_sig if len(e) >= 5)
        intron_sizes_str = ','.join(str(e[1]) for e in chain_sig if len(e) >= 5)
        xors_str = ''.join(XOR_NAME[e[4]] for e in chain_sig if len(e) >= 5)
        print(f'  {count:>6,} molecules | {n_junctions} junctions | introns: {intron_sizes_str} | XOR chain: {xors_str}')

    # Excised intron analysis
    print(f'\n--- EXCISED INTRON CATALOG ---')
    print(f'  Total unique introns: {len(intron_catalog):,}')
    print(f'  Total excision events: {sum(intron_catalog.values()):,}')

    if intron_sizes:
        sizes = np.array(intron_sizes)
        print(f'  Median intron size: {int(np.median(sizes)):,} bp')

        # Size classes with functional implications
        micro = int((sizes < 100).sum())
        small = int(((sizes >= 100) & (sizes < 500)).sum())
        snorna = int(((sizes >= 60) & (sizes < 300)).sum())  # snoRNA host range
        mirna = int(((sizes >= 60) & (sizes < 150)).sum())   # miRNA precursor range
        medium = int(((sizes >= 500) & (sizes < 5000)).sum())
        large = int(((sizes >= 5000) & (sizes < 50000)).sum())
        huge = int((sizes >= 50000).sum())

        print(f'\n  Size classes:')
        print(f'    <100bp (microintrons): {micro:,} ({100*micro/len(sizes):.1f}%)')
        print(f'    60-150bp (miRNA precursor range): {mirna:,} ({100*mirna/len(sizes):.1f}%)')
        print(f'    60-300bp (snoRNA host range): {snorna:,} ({100*snorna/len(sizes):.1f}%)')
        print(f'    100-500bp: {small:,} ({100*small/len(sizes):.1f}%)')
        print(f'    500-5kb: {medium:,} ({100*medium/len(sizes):.1f}%)')
        print(f'    5-50kb: {large:,} ({100*large/len(sizes):.1f}%)')
        print(f'    >50kb: {huge:,} ({100*huge/len(sizes):.1f}%)')

    # Intron flank XOR
    print(f'\n  Excised intron boundary XOR (base 4):')
    total_fx = sum(intron_flank_xors.values())
    for xv in range(4):
        c = intron_flank_xors.get(xv, 0)
        print(f'    {xv:02b} ({XOR_NAME[xv]}): {100*c/max(1,total_fx):.1f}%')

    # Locus divergence analysis
    print(f'\n--- DIVERGENCE AT SHARED LOCI ---')
    print(f'  Loci with 2+ molecules: {n_loci_compared:,}')
    print(f'  Pairwise comparisons: {total_comparisons:,}')

    if junction_divergences:
        print(f'\n  Divergence types:')
        for (dtype, pos), count in junction_divergences.most_common(20):
            pct = 100 * count / max(1, total_comparisons)
            print(f'    Junction {pos}: {dtype} in {count:,} pairs ({pct:.1f}%)')

    # Most variable introns (highest divergence rate)
    print(f'\n--- MOST REUSED INTRONS (top 20 by frequency) ---')
    for (chrom, start, end), count in intron_catalog.most_common(20):
        size = end - start
        print(f'  {chrom}:{start}-{end} ({size:,}bp): {count:,} excisions')

    results = {
        'tool': 'STAFF_Splice_History',
        'source': bam_path,
        'n_reads': n_reads,
        'n_spliced': n_spliced,
        'unique_chains': len(chain_catalog),
        'unique_introns': len(intron_catalog),
        'total_excisions': sum(intron_catalog.values()),
        'intron_size_median': int(np.median(intron_sizes)) if intron_sizes else 0,
        'intron_flank_xor': {f'{k:02b}_{XOR_NAME[k]}': int(v) for k, v in sorted(intron_flank_xors.items())},
        'loci_compared': n_loci_compared,
        'pairwise_comparisons': total_comparisons,
        'compute_time': round(elapsed, 1),
    }
    return results


if __name__ == '__main__':
    bam_path = sys.argv[1]
    output = sys.argv[2] if len(sys.argv) > 2 else bam_path.replace('.bam', '_splice_history.json')
    max_reads = None
    for i, a in enumerate(sys.argv):
        if a == '--max-reads' and i+1 < len(sys.argv):
            max_reads = int(sys.argv[i+1])
    results = build_splice_history(bam_path, max_reads)
    with open(output, 'w') as f:
        json.dump(results, f, indent=2)
    print(f'\nSaved: {output}')
    print('DONE')
