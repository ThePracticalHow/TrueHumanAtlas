"""
STAFF BAM: De Novo Splice Map
Fixed version -- uses shared stream instead of broken external import.
"""
import json, sys, time, os
from collections import defaultdict, Counter
import numpy as np

from .stream import stream_bam


def run(bam_path, max_reads=None):
    """Extract de novo splice map from BAM. Returns result dict."""
    t0 = time.time()
    print(f'STAFF De Novo Splice Map')
    print(f'Input: {bam_path}')

    junction_counts = Counter()
    molecule_patterns = Counter()
    chrom_junctions = defaultdict(Counter)
    n_reads = 0
    n_spliced = 0

    for chrom, pos, mapq, cigar, seq, flag in stream_bam(bam_path):
        n_reads += 1
        if max_reads and n_reads > max_reads:
            break

        ref_pos = pos
        junctions = []
        for op, length in cigar:
            if op in (0, 7, 8):
                ref_pos += length
            elif op == 2:
                ref_pos += length
            elif op == 3:
                junctions.append((ref_pos, ref_pos + length))
                chrom_junctions[chrom][(ref_pos, ref_pos + length)] += 1
                junction_counts[(chrom, ref_pos, ref_pos + length)] += 1
                ref_pos += length

        if junctions:
            n_spliced += 1
            key = (chrom, tuple(sorted(junctions)))
            molecule_patterns[key] += 1

    elapsed = time.time() - t0
    total_junctions = sum(junction_counts.values())

    results = {
        'tool': 'STAFF_Splice_DeNovo',
        'version': '2.0',
        'source': bam_path,
        'n_reads': n_reads,
        'n_spliced': n_spliced,
        'unique_junctions': len(junction_counts),
        'total_junction_observations': total_junctions,
        'unique_molecule_patterns': len(molecule_patterns),
        'per_chrom': {ch: len(js) for ch, js in sorted(chrom_junctions.items())},
        'top_junctions': [
            {'chrom': c, 'start': s, 'end': e, 'count': int(n), 'size': e - s}
            for (c, s, e), n in junction_counts.most_common(100)
        ],
        'compute_time': round(elapsed, 1),
    }
    return results


def run_and_save(bam_path, output_path=None, max_reads=None):
    results = run(bam_path, max_reads=max_reads)
    if output_path:
        with open(output_path, 'w') as f:
            json.dump(results, f, indent=2)
    return results
