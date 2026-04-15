# coding: utf-8
"""
De novo splice map from BAM CIGAR strings. No annotation. No assumptions.

Every molecule defines its own intron/exon structure:
  - M (match) = exon
  - N (skip) = intron/splice junction

For each gene region, collect ALL observed splice junctions.
Then: which junctions are shared? Which are unique to specific cell types?
Are there junctions the databases don't know about?
"""
import struct, zlib, json, sys, time, os
from collections import defaultdict, Counter
import numpy as np

sys.path.insert(0, '/home/jixia')
from bam_stream import stream_bam

def extract_splice_map(bam_path, min_mapq=10):
    """Extract every splice junction from every molecule in the BAM.

    Returns:
      junctions: Counter of (chrom, junction_start, junction_end) -> count
      molecules: list of per-molecule junction sets
      junction_combos: Counter of frozenset(junctions) per molecule -> count
    """
    t0 = time.time()
    print(f'Extracting de novo splice map from {bam_path}...')

    # Per-junction counts
    junction_counts = Counter()  # (chrom, start, end) -> count

    # Per-molecule junction combinations
    # Key = (chrom, sorted tuple of (start, end) pairs)
    molecule_patterns = Counter()

    # Per-chromosome junction catalog
    chrom_junctions = defaultdict(Counter)  # chrom -> {(start, end): count}

    n_reads = 0
    n_spliced = 0
    n_unspliced = 0

    for aln in stream_bam(bam_path):
        n_reads += 1
        if n_reads % 500000 == 0:
            print(f'  {n_reads:,} reads, {n_spliced:,} spliced, {n_unspliced:,} unspliced...')

        if aln['mapq'] < min_mapq:
            continue

        chrom = aln['chrom']
        pos = aln['pos']
        cigar = aln['cigar_ops']

        # Walk the CIGAR to find all N (splice) operations
        ref_pos = pos
        junctions = []

        for op, length in cigar:
            if op in (0, 7, 8):  # M, =, X
                ref_pos += length
            elif op == 1:  # I
                pass
            elif op == 2:  # D
                ref_pos += length
            elif op == 3:  # N - THIS IS A SPLICE JUNCTION
                junction_start = ref_pos
                junction_end = ref_pos + length
                junctions.append((junction_start, junction_end))
                junction_counts[(chrom, junction_start, junction_end)] += 1
                chrom_junctions[chrom][(junction_start, junction_end)] += 1
                ref_pos += length
            elif op == 4:  # S
                pass
            elif op == 5:  # H
                pass

        if junctions:
            n_spliced += 1
            # Record the combination of junctions for this molecule
            pattern = (chrom, tuple(sorted(junctions)))
            molecule_patterns[pattern] += 1
        else:
            n_unspliced += 1

    elapsed = time.time() - t0
    print(f'  Done: {n_reads:,} reads in {elapsed:.1f}s')
    print(f'  Spliced: {n_spliced:,} ({100*n_spliced/max(1,n_reads):.1f}%)')
    print(f'  Unspliced: {n_unspliced:,} ({100*n_unspliced/max(1,n_reads):.1f}%)')
    print(f'  Unique junctions: {len(junction_counts):,}')
    print(f'  Unique molecule patterns: {len(molecule_patterns):,}')

    return junction_counts, molecule_patterns, chrom_junctions


def analyze_junctions(junction_counts, molecule_patterns, chrom_junctions):
    """Analyze the de novo splice landscape."""

    results = {}

    # 1. Junction frequency distribution
    counts = list(junction_counts.values())
    results['total_junctions'] = len(junction_counts)
    results['total_junction_usage'] = sum(counts)
    results['median_junction_count'] = int(np.median(counts))
    results['mean_junction_count'] = round(float(np.mean(counts)), 1)

    # Rare vs common junctions
    rare = sum(1 for c in counts if c <= 2)
    common = sum(1 for c in counts if c >= 100)
    results['rare_junctions_le2'] = rare
    results['common_junctions_ge100'] = common
    results['pct_rare'] = round(100 * rare / max(1, len(counts)), 1)

    print(f'\n=== JUNCTION LANDSCAPE ===')
    print(f'  Total unique junctions: {len(junction_counts):,}')
    print(f'  Rare (<=2 reads): {rare:,} ({results["pct_rare"]:.1f}%)')
    print(f'  Common (>=100 reads): {common:,}')
    print(f'  Median usage: {results["median_junction_count"]}')

    # 2. Junction size distribution (intron length)
    sizes = [end - start for (chrom, start, end) in junction_counts.keys()]
    results['median_intron_size'] = int(np.median(sizes))
    results['mean_intron_size'] = int(np.mean(sizes))

    # Size bins
    tiny = sum(1 for s in sizes if s < 100)       # <100bp - microintrons
    small = sum(1 for s in sizes if 100 <= s < 1000)
    medium = sum(1 for s in sizes if 1000 <= s < 10000)
    large = sum(1 for s in sizes if 10000 <= s < 100000)
    huge = sum(1 for s in sizes if s >= 100000)

    print(f'\n=== INTRON SIZE DISTRIBUTION ===')
    print(f'  <100bp (microintrons): {tiny:,} ({100*tiny/max(1,len(sizes)):.1f}%)')
    print(f'  100-1kb: {small:,} ({100*small/max(1,len(sizes)):.1f}%)')
    print(f'  1-10kb: {medium:,} ({100*medium/max(1,len(sizes)):.1f}%)')
    print(f'  10-100kb: {large:,} ({100*large/max(1,len(sizes)):.1f}%)')
    print(f'  >100kb: {huge:,} ({100*huge/max(1,len(sizes)):.1f}%)')
    results['size_dist'] = {'micro': tiny, 'small': small, 'medium': medium, 'large': large, 'huge': huge}

    # 3. Molecule complexity (junctions per molecule)
    mol_junction_counts = []
    for (chrom, junctions), count in molecule_patterns.items():
        n_junc = len(junctions)
        mol_junction_counts.extend([n_junc] * count)

    if mol_junction_counts:
        results['median_junctions_per_molecule'] = int(np.median(mol_junction_counts))
        results['max_junctions_per_molecule'] = int(np.max(mol_junction_counts))

        print(f'\n=== MOLECULE COMPLEXITY ===')
        print(f'  Median junctions per molecule: {results["median_junctions_per_molecule"]}')
        print(f'  Max junctions per molecule: {results["max_junctions_per_molecule"]}')

        # Distribution
        for n in range(1, min(20, results['max_junctions_per_molecule']+1)):
            c = mol_junction_counts.count(n)
            if c > 0:
                print(f'    {n} junctions: {c:,} molecules')

    # 4. Top molecule patterns (most common splice combinations)
    print(f'\n=== TOP 20 MOLECULE SPLICE PATTERNS ===')
    results['top_patterns'] = []
    for (chrom, junctions), count in molecule_patterns.most_common(20):
        n_junc = len(junctions)
        # Intron sizes
        sizes_str = ','.join(f'{e-s}' for s, e in junctions)
        print(f'  {chrom}:{junctions[0][0]}-{junctions[-1][1]} | {n_junc} junctions | {count:,} molecules | sizes: {sizes_str}')
        results['top_patterns'].append({
            'chrom': chrom,
            'start': junctions[0][0],
            'end': junctions[-1][1],
            'n_junctions': n_junc,
            'count': count,
            'intron_sizes': [e - s for s, e in junctions],
        })

    # 5. Per-chromosome stats
    print(f'\n=== PER-CHROMOSOME JUNCTION COUNTS ===')
    results['per_chrom'] = {}
    for chrom in sorted(chrom_junctions.keys(), key=lambda x: -sum(chrom_junctions[x].values())):
        n_junc = len(chrom_junctions[chrom])
        total = sum(chrom_junctions[chrom].values())
        results['per_chrom'][chrom] = {'unique_junctions': n_junc, 'total_usage': total}
        if n_junc > 100:
            print(f'  {chrom}: {n_junc:,} unique junctions, {total:,} total usage')

    return results


if __name__ == '__main__':
    bam_path = sys.argv[1]
    output = sys.argv[2] if len(sys.argv) > 2 else bam_path.replace('.bam', '_denovo_splice.json')

    junction_counts, molecule_patterns, chrom_junctions = extract_splice_map(bam_path)
    results = analyze_junctions(junction_counts, molecule_patterns, chrom_junctions)
    results['source'] = bam_path

    # Save (junction_counts is too large for JSON, save summary)
    with open(output, 'w') as f:
        json.dump(results, f, indent=2)

    print(f'\nSaved: {output}')
    print('DONE')
