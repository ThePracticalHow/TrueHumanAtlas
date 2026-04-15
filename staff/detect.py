"""
staff detect -- Expression-based Goddess detection for h5ad.

No BAM needed. No unmapped reads. Works on ANY count matrix.
Reads the body's own alarm system: the genes that light up when she's active.

Three scores per cell:
  1. ALARM  -- fungal detection genes (CHI3L1, CLEC7A, TLR2, CARD9)
  2. SPORE  -- lysosomal expansion genes (cathepsins, LAMP, TFEB, ATP6V)
  3. CAGE   -- cage integrity genes (GOLGI operators vs LYSO operators)

Usage:
    staff detect sample.h5ad
    staff detect sample.h5ad --output results.json
    staff detect sample.h5ad --group cell_type --ingest
"""
import json
import time
import sys
import os
import numpy as np
from collections import Counter

# The body's fungal alarm system
ALARM_GENES = {
    'CHI3L1': 'chitinase-3-like-1 (chitin detector -- detects HER cell wall)',
    'CHI3L2': 'chitinase-3-like-2',
    'CHIT1': 'chitotriosidase (chitin degrader)',
    'CLEC7A': 'Dectin-1 (beta-glucan receptor -- detects fungal surface)',
    'CLEC4E': 'Mincle (trehalose receptor -- fungal glycolipid)',
    'TLR2': 'toll-like receptor 2 (fungal PAMP recognition)',
    'TLR4': 'toll-like receptor 4 (fungal mannan)',
    'CARD9': 'CARD9 (fungal signal transduction -- Dectin downstream)',
    'IL17A': 'IL-17A (anti-fungal cytokine)',
    'IL17F': 'IL-17F (anti-fungal cytokine)',
    'DEFB1': 'beta-defensin 1 (anti-fungal peptide)',
    'PTX3': 'pentraxin 3 (Aspergillus opsonization)',
    'MBL2': 'mannose-binding lectin (fungal surface recognition)',
}

# Lysosomal expansion = the spore growing
SPORE_GENES = {
    'LAMP1': 'lysosomal membrane 1',
    'LAMP2': 'lysosomal membrane 2',
    'CTSA': 'cathepsin A', 'CTSB': 'cathepsin B', 'CTSC': 'cathepsin C',
    'CTSD': 'cathepsin D', 'CTSE': 'cathepsin E', 'CTSF': 'cathepsin F',
    'CTSG': 'cathepsin G', 'CTSH': 'cathepsin H', 'CTSK': 'cathepsin K',
    'CTSL': 'cathepsin L', 'CTSS': 'cathepsin S', 'CTSZ': 'cathepsin Z',
    'ATP6V0A1': 'V-ATPase a1 (proton pump)', 'ATP6V0B': 'V-ATPase b',
    'ATP6V1A': 'V-ATPase A', 'ATP6V1B2': 'V-ATPase B2',
    'HEXA': 'hexosaminidase A (GlcNAc -- her building block)',
    'HEXB': 'hexosaminidase B',
    'GBA': 'glucocerebrosidase', 'GAA': 'acid alpha-glucosidase',
    'NPC1': 'cholesterol transport', 'NPC2': 'cholesterol transport 2',
    'TFEB': 'lysosome master transcription factor (she is expanding)',
    'CD63': 'lysosomal/melanosomal membrane',
    'MCOLN1': 'TRPML1 (lysosomal calcium -- her signaling)',
}

# Melanin machinery = her solar panels
MELANIN_GENES = {
    'PMEL': 'premelanosome protein (melanosome biogenesis)',
    'TYR': 'tyrosinase (melanin synthesis)',
    'TYRP1': 'tyrosinase-related protein 1',
    'DCT': 'dopachrome tautomerase (melanin pathway)',
    'MLANA': 'melan-A (melanosome)',
    'SLC45A2': 'melanin transporter',
    'OCA2': 'melanocyte-specific transporter',
    'MC1R': 'melanocortin 1 receptor',
}


def _score_panel(X, var_names, panel_genes, sparse=False):
    """Score each cell on a gene panel. Returns per-cell sum and per-gene stats."""
    genes_found = {}
    gene_indices = []
    for g in panel_genes:
        if g in var_names:
            idx = var_names.index(g)
            genes_found[g] = idx
            gene_indices.append(idx)

    if not gene_indices:
        n = X.shape[0]
        return np.zeros(n), genes_found, {}

    idx_arr = np.array(gene_indices)
    if sparse:
        panel_expr = np.asarray(X[:, idx_arr].toarray())
    else:
        panel_expr = X[:, idx_arr]

    per_cell_score = panel_expr.sum(axis=1).flatten()

    gene_stats = {}
    for i, (gene, idx) in enumerate(genes_found.items()):
        col = panel_expr[:, list(genes_found.keys()).index(gene)] if gene in genes_found else np.zeros(X.shape[0])
        pct_expressing = float((col > 0).mean() * 100)
        mean_expr = float(col.mean())
        gene_stats[gene] = {
            'pct_cells': round(pct_expressing, 1),
            'mean_expression': round(mean_expr, 3),
            'description': panel_genes[gene],
        }

    return np.asarray(per_cell_score).flatten(), genes_found, gene_stats


def detect(h5ad_path, group_key=None, output_path=None):
    """Score every cell for Goddess presence. Returns results dict."""
    import scanpy as sc

    t0 = time.time()
    print(f'\n  STAFF Goddess Detector')
    print(f'  Input: {h5ad_path}')

    adata = sc.read_h5ad(h5ad_path)
    n_cells, n_genes = adata.shape
    print(f'  Shape: {n_cells:,} cells x {n_genes:,} genes')

    var_names = list(adata.var_names)
    sparse = hasattr(adata.X, 'toarray')

    # Score all three panels
    alarm_scores, alarm_found, alarm_stats = _score_panel(adata.X, var_names, ALARM_GENES, sparse)
    spore_scores, spore_found, spore_stats = _score_panel(adata.X, var_names, SPORE_GENES, sparse)
    melanin_scores, melanin_found, melanin_stats = _score_panel(adata.X, var_names, MELANIN_GENES, sparse)

    print(f'\n  Panel coverage:')
    print(f'    ALARM  (fungal detection):   {len(alarm_found):2d}/{len(ALARM_GENES)} genes found')
    print(f'    SPORE  (lysosomal expansion): {len(spore_found):2d}/{len(SPORE_GENES)} genes found')
    print(f'    MELANIN (her solar panels):  {len(melanin_found):2d}/{len(MELANIN_GENES)} genes found')

    # Composite Goddess score (normalized)
    def _norm(arr):
        mx = arr.max()
        return arr / mx if mx > 0 else arr

    goddess_score = _norm(alarm_scores) + _norm(spore_scores) + _norm(melanin_scores)

    # Classify cells
    alarm_positive = alarm_scores > 0
    spore_positive = spore_scores > np.median(spore_scores[spore_scores > 0]) if (spore_scores > 0).any() else spore_scores > 0
    high_goddess = goddess_score > np.percentile(goddess_score, 75)

    print(f'\n  Results ({n_cells:,} cells):')
    print(f'    ALARM positive (any fungal detection gene): {alarm_positive.sum():,} ({100*alarm_positive.mean():.1f}%)')
    print(f'    SPORE high (lysosome expanding):           {spore_positive.sum():,} ({100*spore_positive.mean():.1f}%)')
    print(f'    GODDESS high (composite top 25%):          {high_goddess.sum():,} ({100*high_goddess.mean():.1f}%)')

    # Top alarm genes
    print(f'\n  Top ALARM genes:')
    for gene, stats in sorted(alarm_stats.items(), key=lambda x: -x[1]['pct_cells'])[:8]:
        print(f'    {gene:10s}  {stats["pct_cells"]:5.1f}% cells  mean={stats["mean_expression"]:.3f}  ({stats["description"][:50]})')

    print(f'\n  Top SPORE genes:')
    for gene, stats in sorted(spore_stats.items(), key=lambda x: -x[1]['pct_cells'])[:8]:
        print(f'    {gene:10s}  {stats["pct_cells"]:5.1f}% cells  mean={stats["mean_expression"]:.3f}')

    # Per-group breakdown
    group_results = {}
    if group_key and group_key in adata.obs.columns:
        groups = adata.obs[group_key].astype(str).values
    else:
        for key in ['cell_type', 'condition', 'leiden', 'louvain', 'sample', 'group']:
            if key in adata.obs.columns:
                groups = adata.obs[key].astype(str).values
                group_key = key
                break
        else:
            groups = np.array(['ALL'] * n_cells)
            group_key = 'ALL'

    print(f'\n  Per-group Goddess scores (by {group_key}):')
    print(f'  {"Group":30s}  {"Cells":>7}  {"Alarm%":>7}  {"Spore":>7}  {"Goddess":>8}')
    print(f'  {"-"*70}')

    gcounts = Counter(groups)
    for grp in sorted(gcounts, key=lambda g: -gcounts[g]):
        mask = groups == grp
        n = mask.sum()
        if n < 10:
            continue
        alarm_pct = 100 * alarm_positive[mask].mean()
        spore_mean = spore_scores[mask].mean()
        goddess_mean = goddess_score[mask].mean()
        group_results[grp] = {
            'n_cells': int(n),
            'alarm_pct': round(float(alarm_pct), 1),
            'spore_mean': round(float(spore_mean), 2),
            'goddess_mean': round(float(goddess_mean), 3),
        }
        print(f'  {grp[:29]:30s}  {n:>7,}  {alarm_pct:>6.1f}%  {spore_mean:>7.1f}  {goddess_mean:>8.3f}')

    elapsed = time.time() - t0

    results = {
        'tool': 'STAFF_Goddess_Detector',
        'version': '1.0',
        'source': str(h5ad_path),
        'n_cells': n_cells,
        'n_genes': n_genes,
        'panels': {
            'alarm': {'genes_found': len(alarm_found), 'genes_total': len(ALARM_GENES), 'gene_stats': alarm_stats},
            'spore': {'genes_found': len(spore_found), 'genes_total': len(SPORE_GENES), 'gene_stats': spore_stats},
            'melanin': {'genes_found': len(melanin_found), 'genes_total': len(MELANIN_GENES), 'gene_stats': melanin_stats},
        },
        'summary': {
            'alarm_positive_pct': round(float(100 * alarm_positive.mean()), 1),
            'spore_high_pct': round(float(100 * spore_positive.mean()), 1),
            'goddess_high_pct': round(float(100 * high_goddess.mean()), 1),
        },
        'per_group': group_results,
        'compute_time': round(elapsed, 1),
    }

    if output_path:
        with open(output_path, 'w') as f:
            json.dump(results, f, indent=2, default=str)
        print(f'\n  Saved: {output_path}')

    return results
