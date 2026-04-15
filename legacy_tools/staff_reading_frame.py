# coding: utf-8
"""
STAFF: Reading Frame Analysis
==============================

At each splice junction, check:
1. What is the junction size (intron length) mod 3?
   - 0 = frame-preserving splice
   - 1 or 2 = frame-shifting splice
2. Does the reading frame change across the junction?
3. What is the codon context at the junction boundary?

If cancer changes reading frames differently than fetal,
the protein output is fundamentally altered — not just
which exons are included, but HOW they're read.

A=00 T=01 C=10 G=11
Codons = 6-bit numbers (3 bases x 2 bits = 64 possible codons)
"""
import struct, zlib, json, sys, time, os
from collections import defaultdict, Counter
import numpy as np

BASE_TO_BITS = {'A': 0, 'T': 1, 'C': 2, 'G': 3,
                'a': 0, 't': 1, 'c': 2, 'g': 3, 'N': 0, 'n': 0}
SEQ_DECODE = '=ACMGRSVTWYHKDBN'

CODON_TABLE = {
    'TTT': 'F', 'TTC': 'F', 'TTA': 'L', 'TTG': 'L',
    'CTT': 'L', 'CTC': 'L', 'CTA': 'L', 'CTG': 'L',
    'ATT': 'I', 'ATC': 'I', 'ATA': 'I', 'ATG': 'M',
    'GTT': 'V', 'GTC': 'V', 'GTA': 'V', 'GTG': 'V',
    'TCT': 'S', 'TCC': 'S', 'TCA': 'S', 'TCG': 'S',
    'CCT': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P',
    'ACT': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T',
    'GCT': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A',
    'TAT': 'Y', 'TAC': 'Y', 'TAA': '*', 'TAG': '*',
    'CAT': 'H', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q',
    'AAT': 'N', 'AAC': 'N', 'AAA': 'K', 'AAG': 'K',
    'GAT': 'D', 'GAC': 'D', 'GAA': 'E', 'GAG': 'E',
    'TGT': 'C', 'TGC': 'C', 'TGA': '*', 'TGG': 'W',
    'CGT': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R',
    'AGT': 'S', 'AGC': 'S', 'AGA': 'R', 'AGG': 'R',
    'GGT': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G',
}


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
            if n % 200000 == 0:
                print(f'  {n:,} reads...')


def run(bam_path, max_reads=None):
    t0 = time.time()
    print(f'STAFF Reading Frame Analysis')
    print(f'A=00 T=01 C=10 G=11')
    print(f'Input: {bam_path}\n')

    n_reads = 0
    n_spliced = 0
    n_junctions = 0

    # Junction size mod 3
    frame_shift = Counter()  # 0, 1, 2

    # Codon at junction boundary
    # Last 3 bases before cut = pre-codon
    # First 3 bases after cut = post-codon
    pre_codons = Counter()
    post_codons = Counter()

    # Cross-junction codon (last 1-2 bases before + first 1-2 after = split codon)
    split_codons = Counter()  # the 3-base codon that spans the junction

    # Amino acid at junction
    pre_aa = Counter()
    post_aa = Counter()
    split_aa = Counter()

    # Reading frame position of junction within the exonic sequence
    junction_frame_pos = Counter()  # seq_pos % 3 at junction

    # Start/stop codon proximity
    starts_near_junction = 0  # ATG within 9bp of junction
    stops_near_junction = 0   # TAA/TAG/TGA within 9bp of junction

    # Base-4 XOR at junction (single base, from v2)
    junction_xor = Counter()

    for chrom, pos, cigar, seq in stream_bam(bam_path):
        n_reads += 1
        if max_reads and n_reads > max_reads:
            break

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
            elif op == 3:  # N — SPLICE JUNCTION
                n_junctions += 1

                # Frame shift from intron size
                frame_shift[length % 3] += 1

                # Position in reading frame
                junction_frame_pos[seq_pos % 3] += 1

                # Single-base XOR
                if seq_pos >= 1 and seq_pos < len(seq):
                    lb = BASE_TO_BITS.get(seq[seq_pos-1], 0)
                    rb = BASE_TO_BITS.get(seq[seq_pos], 0)
                    junction_xor[lb ^ rb] += 1

                # Pre-junction codon (last 3 bases of upstream exon)
                if seq_pos >= 3:
                    pre = seq[seq_pos-3:seq_pos]
                    pre_codons[pre] += 1
                    aa = CODON_TABLE.get(pre, '?')
                    pre_aa[aa] += 1

                # Post-junction codon (first 3 bases of downstream exon)
                if seq_pos + 3 <= len(seq):
                    post = seq[seq_pos:seq_pos+3]
                    post_codons[post] += 1
                    aa = CODON_TABLE.get(post, '?')
                    post_aa[aa] += 1

                # Split codon spanning junction
                # Depends on frame position
                fp = seq_pos % 3
                if fp == 0:
                    # Junction at codon boundary — clean split
                    pass
                elif fp == 1:
                    # 1 base from old frame + 2 from new
                    if seq_pos >= 1 and seq_pos + 2 <= len(seq):
                        split = seq[seq_pos-1] + seq[seq_pos:seq_pos+2]
                        split_codons[split] += 1
                        aa = CODON_TABLE.get(split, '?')
                        split_aa[aa] += 1
                elif fp == 2:
                    # 2 bases from old frame + 1 from new
                    if seq_pos >= 2 and seq_pos + 1 <= len(seq):
                        split = seq[seq_pos-2:seq_pos] + seq[seq_pos]
                        split_codons[split] += 1
                        aa = CODON_TABLE.get(split, '?')
                        split_aa[aa] += 1

                # Check for start/stop codons near junction
                context = seq[max(0, seq_pos-9):min(len(seq), seq_pos+9)]
                if 'ATG' in context:
                    starts_near_junction += 1
                if any(stop in context for stop in ['TAA', 'TAG', 'TGA']):
                    stops_near_junction += 1

                junctions.append(seq_pos)
                ref_pos += length
            elif op == 4:  # S
                seq_pos += length

        if junctions:
            n_spliced += 1

    elapsed = time.time() - t0

    print(f'\n{"="*60}')
    print(f'READING FRAME ANALYSIS')
    print(f'{"="*60}')
    print(f'Reads: {n_reads:,} ({elapsed:.1f}s)')
    print(f'Spliced: {n_spliced:,}')
    print(f'Junctions: {n_junctions:,}')

    # Frame shift from intron size
    total_j = sum(frame_shift.values())
    print(f'\n--- INTRON SIZE MOD 3 (reading frame preservation) ---')
    for mod in [0, 1, 2]:
        count = frame_shift.get(mod, 0)
        pct = 100 * count / max(1, total_j)
        label = 'FRAME-PRESERVING' if mod == 0 else f'FRAME-SHIFT +{mod}'
        print(f'  mod 3 = {mod} ({label}): {count:,} ({pct:.1f}%)')

    # Junction position in reading frame
    print(f'\n--- JUNCTION POSITION IN READING FRAME ---')
    total_fp = sum(junction_frame_pos.values())
    for fp in [0, 1, 2]:
        count = junction_frame_pos.get(fp, 0)
        pct = 100 * count / max(1, total_fp)
        label = 'codon boundary' if fp == 0 else f'mid-codon pos {fp}'
        print(f'  pos % 3 = {fp} ({label}): {count:,} ({pct:.1f}%)')

    # Base-4 XOR at junction
    XOR_MEANING = {0: 'IDENTITY', 1: 'COMPLEMENT', 2: 'TRANSVERSION', 3: 'TRANSITION'}
    total_xor = sum(junction_xor.values())
    print(f'\n--- BASE-4 XOR AT JUNCTION ---')
    for xv in range(4):
        count = junction_xor.get(xv, 0)
        pct = 100 * count / max(1, total_xor)
        print(f'  {xv:02b} ({XOR_MEANING[xv]}): {pct:.2f}%')

    # Pre-junction amino acids
    print(f'\n--- AMINO ACID BEFORE JUNCTION (last codon of upstream exon) ---')
    total_pre = sum(pre_aa.values())
    for aa, count in pre_aa.most_common(10):
        pct = 100 * count / max(1, total_pre)
        name = {'*': 'STOP', 'M': 'Met(START)', 'G': 'Gly', 'A': 'Ala', 'V': 'Val',
                'L': 'Leu', 'I': 'Ile', 'P': 'Pro', 'F': 'Phe', 'W': 'Trp',
                'S': 'Ser', 'T': 'Thr', 'C': 'Cys', 'Y': 'Tyr', 'H': 'His',
                'D': 'Asp', 'E': 'Glu', 'N': 'Asn', 'Q': 'Gln', 'K': 'Lys',
                'R': 'Arg'}.get(aa, aa)
        print(f'  {aa} ({name}): {pct:.1f}%')

    # Post-junction amino acids
    print(f'\n--- AMINO ACID AFTER JUNCTION (first codon of downstream exon) ---')
    total_post = sum(post_aa.values())
    for aa, count in post_aa.most_common(10):
        pct = 100 * count / max(1, total_post)
        name = {'*': 'STOP', 'M': 'Met(START)', 'G': 'Gly', 'A': 'Ala', 'V': 'Val',
                'L': 'Leu', 'I': 'Ile', 'P': 'Pro', 'F': 'Phe', 'W': 'Trp',
                'S': 'Ser', 'T': 'Thr', 'C': 'Cys', 'Y': 'Tyr', 'H': 'His',
                'D': 'Asp', 'E': 'Glu', 'N': 'Asn', 'Q': 'Gln', 'K': 'Lys',
                'R': 'Arg'}.get(aa, aa)
        print(f'  {aa} ({name}): {pct:.1f}%')

    # Split codons (spanning junction)
    print(f'\n--- SPLIT CODON AMINO ACID (codon spanning the junction) ---')
    total_split = sum(split_aa.values())
    if total_split > 0:
        for aa, count in split_aa.most_common(10):
            pct = 100 * count / max(1, total_split)
            name = {'*': 'STOP', 'M': 'Met(START)'}.get(aa, aa)
            print(f'  {aa} ({name}): {pct:.1f}%')
        # STOP codons in split position
        stop_pct = 100 * split_aa.get('*', 0) / max(1, total_split)
        print(f'  STOP codons at split: {stop_pct:.2f}%')
    else:
        print(f'  (all junctions at codon boundaries)')

    # Start/stop proximity
    print(f'\n--- START/STOP CODON PROXIMITY ---')
    print(f'  ATG within 9bp of junction: {starts_near_junction:,} ({100*starts_near_junction/max(1,n_junctions):.1f}%)')
    print(f'  Stop within 9bp of junction: {stops_near_junction:,} ({100*stops_near_junction/max(1,n_junctions):.1f}%)')

    results = {
        'tool': 'STAFF_Reading_Frame',
        'source': bam_path,
        'n_reads': n_reads, 'n_spliced': n_spliced, 'n_junctions': n_junctions,
        'intron_mod3': {str(k): int(v) for k, v in sorted(frame_shift.items())},
        'junction_frame_pos': {str(k): int(v) for k, v in sorted(junction_frame_pos.items())},
        'junction_xor_base4': {f'{k:02b}': int(v) for k, v in sorted(junction_xor.items())},
        'pre_aa_top10': {aa: int(c) for aa, c in pre_aa.most_common(10)},
        'post_aa_top10': {aa: int(c) for aa, c in post_aa.most_common(10)},
        'split_aa_top10': {aa: int(c) for aa, c in split_aa.most_common(10)},
        'starts_near_junction': starts_near_junction,
        'stops_near_junction': stops_near_junction,
        'compute_time': round(elapsed, 1),
    }
    return results


if __name__ == '__main__':
    bam_path = sys.argv[1]
    output = sys.argv[2] if len(sys.argv) > 2 else bam_path.replace('.bam', '_frame.json')
    max_reads = None
    for i, a in enumerate(sys.argv):
        if a == '--max-reads' and i+1 < len(sys.argv):
            max_reads = int(sys.argv[i+1])
    results = run(bam_path, max_reads)
    with open(output, 'w') as f:
        json.dump(results, f, indent=2)
    print(f'\nSaved: {output}')
    print('DONE')
