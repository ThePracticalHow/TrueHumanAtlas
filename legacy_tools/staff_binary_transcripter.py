# coding: utf-8
"""
STAFF: Binary Transcripter
==========================
Project 27 — Wings Above Morning

Converts raw BAM alignments to binary transcripts.
No annotation. No database. No interpretation.
The molecule speaks for itself.

Encoding: A=00, T=01, C=10, G=11

For each molecule in the BAM:
  1. Extract the full ATCG sequence from the alignment record
  2. Convert to 2-bit binary
  3. Extract splice structure from CIGAR (N = junction)
  4. Record junction flanking dinucleotides FROM THE MOLECULE
  5. Compute per-locus binary diversity

Outputs:
  - Base composition (binary distribution)
  - Splice junction binary signatures (dinuc flanks + XOR)
  - Per-gene splice signature catalog
  - Locus diversity index
  - Cross-sample comparison when run on multiple BAMs

Usage:
  python staff_binary_transcripter.py <bam_path> [output.json] [--max-reads N]
  python staff_binary_transcripter.py compare <result1.json> <result2.json> ...

Dependencies: None (pure binary BAM parsing, only uses struct + zlib + numpy)
"""
import struct, zlib, json, sys, time, os
from collections import defaultdict, Counter
import numpy as np

# ============================================================
# BINARY ENCODING
# ============================================================

BASE_TO_BITS = {
    'A': 0b00, 'T': 0b01, 'C': 0b10, 'G': 0b11,
    'a': 0b00, 't': 0b01, 'c': 0b10, 'g': 0b11,
    'N': 0b00, 'n': 0b00,
}
BITS_TO_BASE = {0b00: 'A', 0b01: 'T', 0b10: 'C', 0b11: 'G'}

# BAM 4-bit nucleotide encoding
SEQ_DECODE = '=ACMGRSVTWYHKDBN'


def seq_to_binary(seq):
    """Convert ATCG string to numpy uint8 array of 2-bit values."""
    return np.array([BASE_TO_BITS.get(b, 0) for b in seq], dtype=np.uint8)


def decode_bam_seq(encoded_bytes, seq_len):
    """Decode BAM 4-bit packed sequence to ATCG string."""
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


# ============================================================
# BAM STREAMING (pure binary, no libraries)
# ============================================================

def read_bgzf_block(fh):
    """Read one BGZF block, return decompressed bytes or None."""
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
    """Stream BAM yielding alignments with sequence, CIGAR, and position."""
    with open(bam_path, 'rb') as fh:
        # Parse header
        first_block = read_bgzf_block(fh)
        if first_block is None or first_block[:4] != b'BAM\x01':
            raise ValueError(f'Not a valid BAM: {bam_path}')

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
        n_yielded = 0

        while True:
            while len(buffer) < 4:
                block = read_bgzf_block(fh)
                if block is None:
                    return
                buffer += block

            block_size = struct.unpack('<i', buffer[:4])[0]

            while len(buffer) < block_size + 4:
                block = read_bgzf_block(fh)
                if block is None:
                    return
                buffer += block

            record = buffer[4:block_size+4]
            buffer = buffer[block_size+4:]

            if len(record) < 32:
                continue

            ref_id = struct.unpack('<i', record[0:4])[0]
            aln_pos = struct.unpack('<i', record[4:8])[0]
            l_read_name = record[8]
            mapq = record[9]
            n_cigar_op = struct.unpack('<H', record[12:14])[0]
            flag = struct.unpack('<H', record[14:16])[0]
            l_seq = struct.unpack('<i', record[16:20])[0]

            if ref_id < 0 or flag & 4 or mapq < min_mapq:
                continue

            # CIGAR
            cigar_start = 32 + l_read_name
            cigar_end = cigar_start + n_cigar_op * 4
            cigar_ops = []
            for i in range(n_cigar_op):
                val = struct.unpack('<I', record[cigar_start + i*4:cigar_start + (i+1)*4])[0]
                cigar_ops.append((val & 0xF, val >> 4))

            # Sequence
            seq_start = cigar_end
            seq_bytes = (l_seq + 1) // 2
            seq_end = seq_start + seq_bytes
            if seq_end > len(record):
                continue

            seq = decode_bam_seq(record[seq_start:seq_end], l_seq)
            chrom = ref_names[ref_id] if ref_id < len(ref_names) else str(ref_id)

            yield chrom, aln_pos, mapq, cigar_ops, seq, flag

            n_yielded += 1
            if n_yielded % 200000 == 0:
                print(f'  {n_yielded:,} reads processed...')


# ============================================================
# CORE: Binary Transcript Builder
# ============================================================

def run(bam_path, max_reads=None):
    """Build binary transcriptome from BAM."""
    t0 = time.time()
    print(f'STAFF Binary Transcripter')
    print(f'Input: {bam_path}')
    print(f'Encoding: A=00 T=01 C=10 G=11')
    print()

    # Accumulators
    base_counts = np.zeros(4, dtype=np.int64)  # A T C G
    total_bases = 0
    n_reads = 0
    n_spliced = 0
    n_unspliced = 0

    # Junction flank dinucleotides
    junction_flanks = Counter()  # (left_2bp, right_2bp) -> count
    junction_sizes = []

    # Per-locus transcript diversity
    BIN_SIZE = 10000
    locus_hashes = defaultdict(set)  # (chrom, bin) -> set of sequence hashes

    # Splice signature catalog
    splice_sigs = Counter()  # (chrom, tuple_of_junction_positions) -> count

    # Per-chromosome junction sets (for cross-sample comparison)
    chrom_junctions = defaultdict(set)

    for chrom, pos, mapq, cigar, seq, flag in stream_bam(bam_path):
        n_reads += 1
        if max_reads and n_reads > max_reads:
            break

        # Binary base counts
        for b in seq:
            idx = BASE_TO_BITS.get(b, 0)
            base_counts[idx] += 1
        total_bases += len(seq)

        # Walk CIGAR for splice junctions
        ref_pos = pos
        seq_pos = 0
        junctions = []

        for op, length in cigar:
            if op in (0, 7, 8):  # M/=/X
                ref_pos += length
                seq_pos += length
            elif op == 1:  # I
                seq_pos += length
            elif op == 2:  # D
                ref_pos += length
            elif op == 3:  # N — SPLICE
                # Flanking bases from THE MOLECULE
                left = seq[max(0, seq_pos-2):seq_pos]
                right = seq[seq_pos:min(len(seq), seq_pos+2)]
                if len(left) == 2 and len(right) == 2:
                    junction_flanks[(left, right)] += 1
                junctions.append((ref_pos, ref_pos + length))
                junction_sizes.append(length)
                chrom_junctions[chrom].add((ref_pos, ref_pos + length))
                ref_pos += length
            elif op == 4:  # S
                seq_pos += length

        if junctions:
            n_spliced += 1
            splice_sigs[(chrom, tuple(junctions))] += 1
        else:
            n_unspliced += 1

        # Locus diversity (hash of binary sequence)
        locus_key = (chrom, pos // BIN_SIZE)
        locus_hashes[locus_key].add(hash(seq))

    elapsed = time.time() - t0

    # ============================================================
    # RESULTS
    # ============================================================
    print(f'\n{"="*60}')
    print(f'BINARY TRANSCRIPTOME RESULTS')
    print(f'{"="*60}')

    print(f'\nReads: {n_reads:,} in {elapsed:.1f}s')
    print(f'Spliced: {n_spliced:,} ({100*n_spliced/max(1,n_reads):.1f}%)')
    print(f'Unspliced: {n_unspliced:,} ({100*n_unspliced/max(1,n_reads):.1f}%)')
    print(f'Total bases: {total_bases:,}')
    print(f'Unique splice signatures: {len(splice_sigs):,}')
    print(f'Unique junctions: {sum(len(v) for v in chrom_junctions.values()):,}')

    # Base composition
    print(f'\n--- BASE COMPOSITION ---')
    labels = ['A(00)', 'T(01)', 'C(10)', 'G(11)']
    for i in range(4):
        pct = 100 * base_counts[i] / max(1, total_bases)
        bar = '#' * int(pct)
        print(f'  {labels[i]}: {base_counts[i]:>15,} ({pct:5.2f}%) {bar}')

    at = int(base_counts[0] + base_counts[1])
    gc = int(base_counts[2] + base_counts[3])
    gc_pct = 100 * gc / max(1, at + gc)
    print(f'  A+T (0x bits): {at:,} ({100*at/max(1,total_bases):.2f}%)')
    print(f'  G+C (1x bits): {gc:,} ({100*gc/max(1,total_bases):.2f}%)')
    print(f'  GC content: {gc_pct:.2f}%')

    # Bit-level stats
    bit0_ones = int(base_counts[2] + base_counts[3])  # C(10) + G(11) have bit0=1
    bit1_ones = int(base_counts[1] + base_counts[3])  # T(01) + G(11) have bit1=1
    print(f'\n--- BIT-LEVEL ---')
    print(f'  High bit (purine/pyrimidine): {100*bit0_ones/max(1,total_bases):.2f}% ones')
    print(f'  Low bit (amino/keto):         {100*bit1_ones/max(1,total_bases):.2f}% ones')

    # Junction dinucleotides
    print(f'\n--- SPLICE JUNCTION BINARY SIGNATURES ---')
    print(f'  (Flanking dinucleotides from the RNA molecule itself)')
    total_junc = sum(junction_flanks.values())
    print(f'  Total junctions: {total_junc:,}')
    print(f'  {"Left|Right":12s} {"Binary":16s} {"XOR":6s} {"Count":>8s} {"Pct":>6s}')
    for (left, right), count in junction_flanks.most_common(25):
        left_bin = ''.join(f'{BASE_TO_BITS[b]:02b}' for b in left)
        right_bin = ''.join(f'{BASE_TO_BITS[b]:02b}' for b in right)
        xor_val = int(left_bin, 2) ^ int(right_bin, 2)
        pct = 100 * count / max(1, total_junc)
        print(f'  {left}|{right}       {left_bin}|{right_bin}   {xor_val:04b}  {count:>8,} {pct:5.1f}%')

    # XOR distribution across all junctions
    print(f'\n--- XOR DISTRIBUTION ---')
    xor_counts = Counter()
    for (left, right), count in junction_flanks.items():
        if len(left) == 2 and len(right) == 2:
            lb = ''.join(f'{BASE_TO_BITS[b]:02b}' for b in left)
            rb = ''.join(f'{BASE_TO_BITS[b]:02b}' for b in right)
            xor_val = int(lb, 2) ^ int(rb, 2)
            xor_counts[xor_val] += count
    for xor_val in sorted(xor_counts.keys()):
        count = xor_counts[xor_val]
        pct = 100 * count / max(1, total_junc)
        print(f'  XOR={xor_val:04b} ({xor_val:>2d}): {count:>8,} ({pct:5.1f}%)')

    # Junction size distribution
    if junction_sizes:
        sizes = np.array(junction_sizes)
        print(f'\n--- JUNCTION SIZE DISTRIBUTION ---')
        print(f'  Median: {int(np.median(sizes)):,} bp')
        print(f'  Mean: {int(np.mean(sizes)):,} bp')
        bins = [(0, 100, 'micro (<100bp)'), (100, 1000, 'small (100-1kb)'),
                (1000, 10000, 'medium (1-10kb)'), (10000, 100000, 'large (10-100kb)'),
                (100000, 10000000, 'huge (>100kb)')]
        for lo, hi, label in bins:
            c = int(((sizes >= lo) & (sizes < hi)).sum())
            print(f'  {label}: {c:,} ({100*c/len(sizes):.1f}%)')

    # Locus diversity
    n_loci = len(locus_hashes)
    diversities = [len(v) for v in locus_hashes.values()]
    unique_pct = sum(1 for d in diversities if d >= 2) / max(1, n_loci) * 100
    print(f'\n--- LOCUS DIVERSITY ---')
    print(f'  Active loci (10kb bins): {n_loci:,}')
    print(f'  Mean unique sequences per locus: {np.mean(diversities):.1f}')
    print(f'  Loci with >1 unique sequence: {unique_pct:.1f}%')

    # Build output
    results = {
        'tool': 'STAFF_Binary_Transcripter',
        'version': '1.0',
        'source': bam_path,
        'n_reads': n_reads,
        'n_spliced': n_spliced,
        'n_unspliced': n_unspliced,
        'total_bases': total_bases,
        'base_composition': {
            'A_00': int(base_counts[0]), 'T_01': int(base_counts[1]),
            'C_10': int(base_counts[2]), 'G_11': int(base_counts[3]),
        },
        'gc_content': round(gc_pct, 4),
        'bit_high_pct': round(100 * bit0_ones / max(1, total_bases), 4),
        'bit_low_pct': round(100 * bit1_ones / max(1, total_bases), 4),
        'unique_splice_sigs': len(splice_sigs),
        'unique_junctions': sum(len(v) for v in chrom_junctions.values()),
        'junction_dinucs': {f'{l}|{r}': int(c) for (l, r), c in junction_flanks.most_common(50)},
        'xor_distribution': {f'{k:04b}': int(v) for k, v in sorted(xor_counts.items())},
        'junction_size_median': int(np.median(junction_sizes)) if junction_sizes else 0,
        'n_active_loci': n_loci,
        'per_chrom_junctions': {ch: len(js) for ch, js in sorted(chrom_junctions.items())},
        'compute_time': round(elapsed, 1),
    }

    return results


def compare(result_files):
    """Compare binary transcriptomes across samples."""
    print(f'STAFF Binary Transcripter — COMPARE')
    print(f'Comparing {len(result_files)} samples\n')

    results = []
    for path in result_files:
        with open(path) as f:
            results.append(json.load(f))

    # Header
    names = [os.path.basename(r['source']).replace('.bam', '') for r in results]
    print(f'{"Metric":30s}', end='')
    for name in names:
        print(f' {name[:15]:>15s}', end='')
    print()
    print('-' * (30 + 16 * len(names)))

    # Metrics
    for metric, key, fmt in [
        ('Reads', 'n_reads', ',d'),
        ('Spliced %', lambda r: 100*r['n_spliced']/max(1,r['n_reads']), '.1f'),
        ('GC content %', 'gc_content', '.2f'),
        ('High bit %', 'bit_high_pct', '.2f'),
        ('Low bit %', 'bit_low_pct', '.2f'),
        ('Unique junctions', 'unique_junctions', ',d'),
        ('Unique splice sigs', 'unique_splice_sigs', ',d'),
        ('Junction size median', 'junction_size_median', ',d'),
    ]:
        print(f'{metric:30s}', end='')
        for r in results:
            if callable(key):
                val = key(r)
            else:
                val = r.get(key, 0)
            print(f' {val:>15{fmt}}', end='')
        print()

    # XOR distribution comparison
    print(f'\n--- XOR DISTRIBUTION COMPARISON ---')
    all_xors = set()
    for r in results:
        all_xors.update(r.get('xor_distribution', {}).keys())

    print(f'{"XOR":6s}', end='')
    for name in names:
        print(f' {name[:12]:>12s}', end='')
    print()

    for xor_key in sorted(all_xors):
        print(f'{xor_key:6s}', end='')
        for r in results:
            total = sum(r.get('xor_distribution', {}).values())
            count = r.get('xor_distribution', {}).get(xor_key, 0)
            pct = 100 * count / max(1, total)
            print(f' {pct:11.1f}%', end='')
        print()


# ============================================================
# CLI
# ============================================================

if __name__ == '__main__':
    if len(sys.argv) < 2:
        print(__doc__)
        sys.exit(1)

    if sys.argv[1] == 'compare':
        compare(sys.argv[2:])
    else:
        bam_path = sys.argv[1]
        output = sys.argv[2] if len(sys.argv) > 2 else bam_path.replace('.bam', '_binary_transcript.json')

        max_reads = None
        for i, arg in enumerate(sys.argv):
            if arg == '--max-reads' and i + 1 < len(sys.argv):
                max_reads = int(sys.argv[i + 1])

        results = run(bam_path, max_reads=max_reads)

        with open(output, 'w') as f:
            json.dump(results, f, indent=2)
        print(f'\nSaved: {output}')
        print('DONE')
