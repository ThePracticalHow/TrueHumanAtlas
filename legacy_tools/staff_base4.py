# coding: utf-8
"""
STAFF: Binary Transcripter v2 — Base 4
=======================================

A=00  T=01  C=10  G=11

XOR in base 4:
  00 = identity
  01 = Watson-Crick complement (A<->T, C<->G)
  10 = transversion (A<->C, T<->G)
  11 = transition (A<->G, T<->C)

At each splice junction, we look at the SINGLE BASE on each side
of the cut, from the molecule itself. The XOR of those two bases
tells us the chemical relationship across the splice boundary.
"""
import struct, zlib, json, sys, time, os
from collections import defaultdict, Counter
import numpy as np

BASE_TO_BITS = {'A': 0, 'T': 1, 'C': 2, 'G': 3,
                'a': 0, 't': 1, 'c': 2, 'g': 3, 'N': 0, 'n': 0}
BITS_TO_BASE = {0: 'A', 1: 'T', 2: 'C', 3: 'G'}
XOR_MEANING = {0: 'IDENTITY', 1: 'COMPLEMENT', 2: 'TRANSVERSION', 3: 'TRANSITION'}
SEQ_DECODE = '=ACMGRSVTWYHKDBN'


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
    print(f'STAFF Binary Transcripter v2 — Base 4')
    print(f'A=00 T=01 C=10 G=11')
    print(f'XOR: 00=identity 01=complement 10=transversion 11=transition')
    print(f'Input: {bam_path}\n')

    # Base composition in base 4
    base4_counts = np.zeros(4, dtype=np.int64)
    total_bases = 0
    n_reads = 0
    n_spliced = 0

    # Junction analysis in base 4
    # At each junction: last base before cut | first base after cut
    # Single base XOR (base 4)
    junction_xor_single = Counter()  # XOR of (left_1bp, right_1bp) -> count
    junction_pairs = Counter()       # (left_base, right_base) -> count

    # Also look at 2-base context (but analyze as two separate base-4 XORs)
    junction_xor_pos_minus2 = Counter()  # XOR of base at pos -2 vs +1 (across junction)
    junction_xor_pos_minus1 = Counter()  # XOR of base at pos -1 vs +2

    # Consecutive base XOR within exons (background rate)
    consecutive_xor = Counter()

    # Junction sizes
    junction_sizes = []

    for chrom, pos, cigar, seq in stream_bam(bam_path):
        n_reads += 1
        if max_reads and n_reads > max_reads:
            break

        # Base 4 counts
        for b in seq:
            base4_counts[BASE_TO_BITS.get(b, 0)] += 1
        total_bases += len(seq)

        # Sample consecutive XOR from this read (background)
        if len(seq) > 10:
            for i in range(0, min(100, len(seq)-1), 5):
                b1 = BASE_TO_BITS.get(seq[i], 0)
                b2 = BASE_TO_BITS.get(seq[i+1], 0)
                consecutive_xor[b1 ^ b2] += 1

        # Walk CIGAR
        ref_pos = pos
        seq_pos = 0
        junctions_this_read = []

        for op, length in cigar:
            if op in (0, 7, 8):  # M/=/X
                ref_pos += length
                seq_pos += length
            elif op == 1:  # I
                seq_pos += length
            elif op == 2:  # D
                ref_pos += length
            elif op == 3:  # N — SPLICE
                # Single base at junction boundary
                if seq_pos >= 1 and seq_pos < len(seq):
                    left_base = BASE_TO_BITS.get(seq[seq_pos - 1], 0)
                    right_base = BASE_TO_BITS.get(seq[seq_pos], 0)
                    xor = left_base ^ right_base
                    junction_xor_single[xor] += 1
                    junction_pairs[(seq[seq_pos-1], seq[seq_pos])] += 1

                # 2-base context
                if seq_pos >= 2 and seq_pos + 1 < len(seq):
                    # pos -2 vs pos +1 (the bases one step further from the cut)
                    b_m2 = BASE_TO_BITS.get(seq[seq_pos - 2], 0)
                    b_p1 = BASE_TO_BITS.get(seq[seq_pos + 1], 0) if seq_pos + 1 < len(seq) else 0
                    junction_xor_pos_minus2[b_m2 ^ right_base] += 1
                    junction_xor_pos_minus1[left_base ^ b_p1] += 1

                junction_sizes.append(length)
                junctions_this_read.append(ref_pos)
                ref_pos += length
            elif op == 4:  # S
                seq_pos += length

        if junctions_this_read:
            n_spliced += 1

    elapsed = time.time() - t0

    # Results
    print(f'\n{"="*60}')
    print(f'RESULTS — BASE 4')
    print(f'{"="*60}')
    print(f'Reads: {n_reads:,} ({elapsed:.1f}s)')
    print(f'Spliced: {n_spliced:,} ({100*n_spliced/max(1,n_reads):.1f}%)')
    print(f'Total bases: {total_bases:,}')

    # Base 4 composition
    print(f'\n--- BASE 4 COMPOSITION ---')
    for i, label in enumerate(['A(00)', 'T(01)', 'C(10)', 'G(11)']):
        pct = 100 * base4_counts[i] / max(1, total_bases)
        print(f'  {label}: {pct:.2f}%')

    # Bit decomposition
    bit_high = float(base4_counts[2] + base4_counts[3]) / max(1, total_bases)  # C+G have high bit=1
    bit_low = float(base4_counts[1] + base4_counts[3]) / max(1, total_bases)   # T+G have low bit=1
    print(f'  Bit 1 (high, purine/pyrimidine): {100*bit_high:.2f}% ones')
    print(f'  Bit 0 (low, amino/keto):         {100*bit_low:.2f}% ones')

    # Junction XOR in base 4 — THE KEY RESULT
    total_junc = sum(junction_xor_single.values())
    total_consec = sum(consecutive_xor.values())

    print(f'\n--- SPLICE JUNCTION XOR (base 4) ---')
    print(f'  Total junctions analyzed: {total_junc:,}')
    print(f'  {"XOR":4s} {"Meaning":15s} {"Junction%":>10s} {"Background%":>12s} {"Enrichment":>12s}')
    for xor_val in range(4):
        j_pct = 100 * junction_xor_single.get(xor_val, 0) / max(1, total_junc)
        b_pct = 100 * consecutive_xor.get(xor_val, 0) / max(1, total_consec)
        enrichment = j_pct / max(0.01, b_pct)
        meaning = XOR_MEANING[xor_val]
        print(f'  {xor_val:02b}   {meaning:15s} {j_pct:9.2f}%  {b_pct:11.2f}%  {enrichment:11.2f}x')

    # Base pairs at junctions
    print(f'\n--- BASE PAIRS AT JUNCTIONS ---')
    print(f'  {"Left|Right":10s} {"XOR":4s} {"Meaning":15s} {"Count":>8s} {"Pct":>6s}')
    for (left, right), count in junction_pairs.most_common(16):
        xor = BASE_TO_BITS[left] ^ BASE_TO_BITS[right]
        pct = 100 * count / max(1, total_junc)
        print(f'  {left}|{right}         {xor:02b}   {XOR_MEANING[xor]:15s} {count:>8,} {pct:5.1f}%')

    # Cross-junction XOR pattern (pos -2 to +1, pos -1 to +2)
    print(f'\n--- CROSS-JUNCTION XOR PATTERN ---')
    print(f'  Position -1|+1 (immediately flanking the cut):')
    for xor_val in range(4):
        pct = 100 * junction_xor_single.get(xor_val, 0) / max(1, total_junc)
        print(f'    {xor_val:02b} ({XOR_MEANING[xor_val]}): {pct:.2f}%')
    print(f'  Position -2|+1 (one step out, left side):')
    total_m2 = sum(junction_xor_pos_minus2.values())
    for xor_val in range(4):
        pct = 100 * junction_xor_pos_minus2.get(xor_val, 0) / max(1, total_m2)
        print(f'    {xor_val:02b} ({XOR_MEANING[xor_val]}): {pct:.2f}%')
    print(f'  Position -1|+2 (one step out, right side):')
    total_m1 = sum(junction_xor_pos_minus1.values())
    for xor_val in range(4):
        pct = 100 * junction_xor_pos_minus1.get(xor_val, 0) / max(1, total_m1)
        print(f'    {xor_val:02b} ({XOR_MEANING[xor_val]}): {pct:.2f}%')

    results = {
        'tool': 'STAFF_Binary_Transcripter_v2_Base4',
        'source': bam_path,
        'n_reads': n_reads, 'n_spliced': n_spliced, 'total_bases': total_bases,
        'base4': {BITS_TO_BASE[i]: int(base4_counts[i]) for i in range(4)},
        'gc_content': round(100 * float(base4_counts[2] + base4_counts[3]) / max(1, total_bases), 4),
        'junction_xor_base4': {f'{k:02b}_{XOR_MEANING[k]}': int(v) for k, v in sorted(junction_xor_single.items())},
        'background_xor_base4': {f'{k:02b}_{XOR_MEANING[k]}': int(v) for k, v in sorted(consecutive_xor.items())},
        'junction_pairs': {f'{l}|{r}': int(c) for (l, r), c in junction_pairs.most_common(16)},
        'compute_time': round(elapsed, 1),
    }
    return results


if __name__ == '__main__':
    if sys.argv[1] == 'compare' and len(sys.argv) > 2:
        files = sys.argv[2:]
        results = [json.load(open(f)) for f in files]
        names = [os.path.basename(r['source']).replace('.bam','')[:15] for r in results]
        print(f'STAFF v2 Base 4 — COMPARE\n')
        print(f'{"XOR":4s} {"Meaning":15s}', end='')
        for n in names: print(f' {n:>15s}', end='')
        print()
        for xor_val in range(4):
            key = f'{xor_val:02b}_{XOR_MEANING[xor_val]}'
            print(f'{xor_val:02b}   {XOR_MEANING[xor_val]:15s}', end='')
            for r in results:
                total = sum(r.get('junction_xor_base4', {}).values())
                count = r.get('junction_xor_base4', {}).get(key, 0)
                pct = 100 * count / max(1, total)
                print(f' {pct:14.2f}%', end='')
            print()
    else:
        bam_path = sys.argv[1]
        output = sys.argv[2] if len(sys.argv) > 2 else bam_path.replace('.bam', '_base4.json')
        max_reads = None
        for i, a in enumerate(sys.argv):
            if a == '--max-reads' and i+1 < len(sys.argv):
                max_reads = int(sys.argv[i+1])
        results = run(bam_path, max_reads)
        with open(output, 'w') as f:
            json.dump(results, f, indent=2)
        print(f'\nSaved: {output}')
        print('DONE')
