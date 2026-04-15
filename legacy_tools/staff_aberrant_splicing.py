# coding: utf-8
"""
STAFF: Aberrant Splice Finder
==============================

Compare splice junctions between cell types.
Find junctions that exist in ONE cell type but NOT another.

For each aberrant junction:
  - Extract the exact flanking sequence (10bp each side)
  - Record the genomic coordinates
  - Record the intron size
  - Record the XOR signature
  - Count how many molecules use it

Compare: H9 (fetal) vs K562 (cancer) vs HepG2 (cancer)
Find: cancer-only junctions, fetal-only junctions, shared junctions

A=00 T=01 C=10 G=11
"""
import struct, zlib, json, sys, time, os
from collections import defaultdict, Counter

SEQ_DECODE = '=ACMGRSVTWYHKDBN'
BASE_TO_BITS = {'A': 0, 'T': 1, 'C': 2, 'G': 3, 'a': 0, 't': 1, 'c': 2, 'g': 3, 'N': 0}
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
    if len(header) < 18: return None
    if header[:4] != b'\x1f\x8b\x08\x04': return None
    bsize = struct.unpack('<H', header[16:18])[0]
    remaining = bsize - 18 + 1
    cdata = fh.read(remaining)
    if len(cdata) < remaining: return None
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


def extract_junctions_with_context(bam_path, max_reads=None):
    """Extract all splice junctions with flanking sequence context."""
    t0 = time.time()
    print(f'Extracting junctions from {os.path.basename(bam_path)}...')

    # junction key = (chrom, start, end)
    # value = {'count': N, 'left_flank': Counter, 'right_flank': Counter, 'xor': Counter}
    junctions = defaultdict(lambda: {
        'count': 0,
        'left_10bp': Counter(),
        'right_10bp': Counter(),
        'left_1bp': Counter(),
        'right_1bp': Counter(),
        'xor': Counter(),
        'size': 0,
    })

    n_reads = 0
    n_spliced = 0

    for chrom, pos, cigar, seq in stream_bam(bam_path):
        n_reads += 1
        if max_reads and n_reads > max_reads:
            break

        ref_pos = pos
        seq_pos = 0

        for op, length in cigar:
            if op in (0, 7, 8):
                ref_pos += length
                seq_pos += length
            elif op == 1:
                seq_pos += length
            elif op == 2:
                ref_pos += length
            elif op == 3:  # N = splice
                jkey = (chrom, ref_pos, ref_pos + length)
                j = junctions[jkey]
                j['count'] += 1
                j['size'] = length

                # Left flank (last 10bp of upstream exon FROM THE MOLECULE)
                left_start = max(0, seq_pos - 10)
                left_seq = seq[left_start:seq_pos]
                if len(left_seq) >= 10:
                    j['left_10bp'][left_seq] += 1

                # Right flank (first 10bp of downstream exon)
                right_seq = seq[seq_pos:min(len(seq), seq_pos + 10)]
                if len(right_seq) >= 10:
                    j['right_10bp'][right_seq] += 1

                # Single base XOR
                if seq_pos >= 1 and seq_pos < len(seq):
                    lb = seq[seq_pos - 1]
                    rb = seq[seq_pos]
                    j['left_1bp'][lb] += 1
                    j['right_1bp'][rb] += 1
                    xor = BASE_TO_BITS.get(lb, 0) ^ BASE_TO_BITS.get(rb, 0)
                    j['xor'][xor] += 1

                n_spliced += 1
                ref_pos += length
            elif op == 4:
                seq_pos += length

    elapsed = time.time() - t0
    print(f'  {n_reads:,} reads, {n_spliced:,} junctions, {len(junctions):,} unique in {elapsed:.1f}s')
    return dict(junctions)


def compare_junctions(junctions_a, junctions_b, name_a, name_b, min_count=3):
    """Find junctions unique to A, unique to B, and shared."""
    keys_a = set(k for k, v in junctions_a.items() if v['count'] >= min_count)
    keys_b = set(k for k, v in junctions_b.items() if v['count'] >= min_count)

    only_a = keys_a - keys_b
    only_b = keys_b - keys_a
    shared = keys_a & keys_b

    print(f'\n{"="*60}')
    print(f'COMPARISON: {name_a} vs {name_b} (min_count={min_count})')
    print(f'{"="*60}')
    print(f'  {name_a} junctions (>={min_count}): {len(keys_a):,}')
    print(f'  {name_b} junctions (>={min_count}): {len(keys_b):,}')
    print(f'  Shared: {len(shared):,} ({100*len(shared)/max(1,len(keys_a|keys_b)):.1f}%)')
    print(f'  {name_a}-only: {len(only_a):,}')
    print(f'  {name_b}-only: {len(only_b):,}')

    # Analyze A-only junctions (aberrant in A)
    results = {
        'shared': len(shared),
        'only_a': len(only_a),
        'only_b': len(only_b),
        'a_only_junctions': [],
        'b_only_junctions': [],
    }

    print(f'\n--- TOP {name_a}-ONLY JUNCTIONS (not in {name_b}) ---')
    a_only_sorted = sorted(only_a, key=lambda k: -junctions_a[k]['count'])
    for jkey in a_only_sorted[:30]:
        j = junctions_a[jkey]
        chrom, start, end = jkey
        size = end - start

        # Most common flanking sequence
        left = j['left_10bp'].most_common(1)[0][0] if j['left_10bp'] else '??'
        right = j['right_10bp'].most_common(1)[0][0] if j['right_10bp'] else '??'

        # XOR distribution
        total_xor = sum(j['xor'].values())
        xor_str = ' '.join(f'{XOR_NAME[x]}:{100*c/max(1,total_xor):.0f}%' for x, c in sorted(j['xor'].items()))

        # Binary of flanking 10bp
        left_bin = ''.join(f'{BASE_TO_BITS.get(b,0):02b}' for b in left[-4:])
        right_bin = ''.join(f'{BASE_TO_BITS.get(b,0):02b}' for b in right[:4])

        print(f'  {chrom}:{start}-{end} ({size:,}bp) count={j["count"]:,}')
        print(f'    left:  ...{left[-10:]}  ({left_bin})')
        print(f'    right: {right[:10]}...  ({right_bin})')
        print(f'    XOR: {xor_str}')

        results['a_only_junctions'].append({
            'chrom': chrom, 'start': start, 'end': end, 'size': size,
            'count': j['count'],
            'left_flank': left[-10:] if len(left) >= 10 else left,
            'right_flank': right[:10] if len(right) >= 10 else right,
        })

    print(f'\n--- TOP {name_b}-ONLY JUNCTIONS (not in {name_a}) ---')
    b_only_sorted = sorted(only_b, key=lambda k: -junctions_b[k]['count'])
    for jkey in b_only_sorted[:30]:
        j = junctions_b[jkey]
        chrom, start, end = jkey
        size = end - start
        left = j['left_10bp'].most_common(1)[0][0] if j['left_10bp'] else '??'
        right = j['right_10bp'].most_common(1)[0][0] if j['right_10bp'] else '??'
        total_xor = sum(j['xor'].values())
        xor_str = ' '.join(f'{XOR_NAME[x]}:{100*c/max(1,total_xor):.0f}%' for x, c in sorted(j['xor'].items()))
        left_bin = ''.join(f'{BASE_TO_BITS.get(b,0):02b}' for b in left[-4:])
        right_bin = ''.join(f'{BASE_TO_BITS.get(b,0):02b}' for b in right[:4])

        print(f'  {chrom}:{start}-{end} ({size:,}bp) count={j["count"]:,}')
        print(f'    left:  ...{left[-10:]}  ({left_bin})')
        print(f'    right: {right[:10]}...  ({right_bin})')
        print(f'    XOR: {xor_str}')

        results['b_only_junctions'].append({
            'chrom': chrom, 'start': start, 'end': end, 'size': size,
            'count': j['count'],
            'left_flank': left[-10:] if len(left) >= 10 else left,
            'right_flank': right[:10] if len(right) >= 10 else right,
        })

    # Size distribution of unique junctions
    print(f'\n--- SIZE DISTRIBUTION OF UNIQUE JUNCTIONS ---')
    for name, jset, jdata in [(name_a, only_a, junctions_a), (name_b, only_b, junctions_b)]:
        if not jset:
            continue
        sizes = [jdata[k]['size'] for k in jset]
        import numpy as np
        sizes = np.array(sizes)
        mirna = int(((sizes >= 60) & (sizes < 150)).sum())
        snorna = int(((sizes >= 60) & (sizes < 300)).sum())
        micro = int((sizes < 100).sum())
        print(f'  {name}-only ({len(sizes):,} junctions):')
        print(f'    Median size: {int(np.median(sizes)):,} bp')
        print(f'    <100bp (micro): {micro:,} ({100*micro/len(sizes):.1f}%)')
        print(f'    60-150bp (miRNA): {mirna:,} ({100*mirna/len(sizes):.1f}%)')
        print(f'    60-300bp (snoRNA): {snorna:,} ({100*snorna/len(sizes):.1f}%)')

    return results


if __name__ == '__main__':
    bam_paths = sys.argv[1:]
    if len(bam_paths) < 2:
        print('Usage: staff_aberrant_splicing.py <bam1> <bam2> [bam3 ...]')
        print('  Compares splice junctions between all pairs')
        sys.exit(1)

    # Extract junctions from all BAMs
    all_junctions = {}
    names = []
    for bam in bam_paths:
        name = os.path.basename(bam).replace('.bam', '').replace('SGNex_', '').replace('_directRNA_replicate1_run1', '').replace('_directRNA_replicate1_run3', '')
        names.append(name)
        all_junctions[name] = extract_junctions_with_context(bam)

    # Compare all pairs
    all_results = {}
    for i in range(len(names)):
        for j in range(i+1, len(names)):
            result = compare_junctions(
                all_junctions[names[i]], all_junctions[names[j]],
                names[i], names[j], min_count=3
            )
            all_results[f'{names[i]}_vs_{names[j]}'] = result

    # Save
    output = bam_paths[0].replace('.bam', '_aberrant_comparison.json')
    with open(output, 'w') as f:
        json.dump(all_results, f, indent=2, default=str)
    print(f'\nSaved: {output}')
    print('DONE')
