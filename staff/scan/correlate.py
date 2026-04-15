"""
Correlate fungal transcript load with coupling tensor per cell.

The prediction: cells with more fungal transcripts will show higher RIBO_indep.
The Goddess arrives before the Red Queen activates.
"""
import json
import numpy as np
from scipy.stats import spearmanr
from collections import defaultdict


def correlate_with_h5ad(scan_results, h5ad_path, group_key=None):
    """
    Match per-cell fungal counts to per-cell coupling tensor.
    
    Args:
        scan_results: dict from scanner.scan_bam() with 'per_cell' key
        h5ad_path: path to h5ad file (same sample)
        group_key: optional column in obs to group by
    
    Returns: correlation results dict
    """
    import scanpy as sc

    print(f'\n  Correlating fungal load with coupling tensor...')
    print(f'  h5ad: {h5ad_path}')

    adata = sc.read_h5ad(h5ad_path)
    print(f'  Cells in h5ad: {adata.n_obs:,}')

    per_cell = scan_results.get('per_cell', {})
    print(f'  Cells with fungal reads: {len(per_cell):,}')

    obs_names = set(adata.obs_names)
    matched = 0
    fungal_load = np.zeros(adata.n_obs, dtype=np.float64)

    for i, cell_id in enumerate(adata.obs_names):
        if cell_id in per_cell:
            fungal_load[i] = sum(per_cell[cell_id].values())
            matched += 1

    if matched == 0:
        cb_cells = {cb.split('-')[0]: cb for cb in per_cell}
        for i, cell_id in enumerate(adata.obs_names):
            short = cell_id.split('-')[0]
            if short in cb_cells:
                fungal_load[i] = sum(per_cell[cb_cells[short]].values())
                matched += 1

    print(f'  Matched barcodes: {matched:,}')

    if matched < 10:
        print(f'  WARNING: Very few matches. Barcode formats may not align.')
        return {'matched': matched, 'error': 'insufficient_barcode_matches'}

    genes = list(adata.var_names)
    ribo_mask = np.array([g.startswith('RPS') or g.startswith('RPL') for g in genes])
    mito_mask = np.array([g.startswith('MT-') or g.startswith('mt-') for g in genes])

    X = adata.X
    if hasattr(X, 'toarray'):
        ribo = np.asarray(X[:, ribo_mask].sum(axis=1)).flatten()
        mito = np.asarray(X[:, mito_mask].sum(axis=1)).flatten()
        total = np.asarray(X.sum(axis=1)).flatten()
    else:
        ribo = X[:, ribo_mask].sum(axis=1)
        mito = X[:, mito_mask].sum(axis=1)
        total = X.sum(axis=1)

    ribo_frac = ribo / np.maximum(total, 1)

    has_fungus = fungal_load > 0
    n_with = has_fungus.sum()
    n_without = (~has_fungus).sum()

    print(f'\n  Cells with fungal RNA:    {n_with:,}')
    print(f'  Cells without fungal RNA: {n_without:,}')

    results = {
        'n_cells': int(adata.n_obs),
        'n_matched': int(matched),
        'n_with_fungal': int(n_with),
        'n_without_fungal': int(n_without),
    }

    if n_with > 0 and n_without > 0:
        ribo_with = ribo_frac[has_fungus].mean()
        ribo_without = ribo_frac[~has_fungus].mean()
        mito_with = (mito / np.maximum(total, 1))[has_fungus].mean()
        mito_without = (mito / np.maximum(total, 1))[~has_fungus].mean()

        print(f'\n  RIBO fraction (with fungus):    {ribo_with:.4f}')
        print(f'  RIBO fraction (without fungus): {ribo_without:.4f}')
        print(f'  MITO fraction (with fungus):    {mito_with:.4f}')
        print(f'  MITO fraction (without fungus): {mito_without:.4f}')

        if n_with >= 5:
            rho_ribo, p_ribo = spearmanr(fungal_load[has_fungus], ribo_frac[has_fungus])
            print(f'\n  Spearman (fungal load vs RIBO fraction): rho={rho_ribo:.4f}, p={p_ribo:.2e}')
            results['spearman_fungal_vs_ribo'] = {'rho': float(rho_ribo), 'p': float(p_ribo)}

        results['ribo_fraction_with_fungus'] = float(ribo_with)
        results['ribo_fraction_without_fungus'] = float(ribo_without)
        results['mito_fraction_with_fungus'] = float(mito_with)
        results['mito_fraction_without_fungus'] = float(mito_without)
        results['delta_ribo'] = float(ribo_with - ribo_without)

        direction = "HIGHER" if ribo_with > ribo_without else "LOWER"
        print(f'\n  INTERPRETATION: Cells with fungal RNA have {direction} ribosomal fraction')
        if ribo_with > ribo_without:
            print(f'    -> Fungal-positive cells are MORE coupled (Green Goddess territory)')
        else:
            print(f'    -> Fungal-positive cells are LESS coupled (Red Queen activated)')

    if group_key and group_key in adata.obs.columns:
        groups = adata.obs[group_key].astype(str).values
        group_stats = {}
        for grp in sorted(set(groups)):
            mask = groups == grp
            n_grp = mask.sum()
            n_fungal = (has_fungus & mask).sum()
            if n_grp > 0:
                group_stats[grp] = {
                    'n_cells': int(n_grp),
                    'n_fungal': int(n_fungal),
                    'fungal_pct': round(100 * n_fungal / n_grp, 2),
                    'mean_ribo_frac': float(ribo_frac[mask].mean()),
                }
        results['per_group'] = group_stats
        print(f'\n  Per-group fungal rates:')
        for grp, s in sorted(group_stats.items(), key=lambda x: -x[1]['fungal_pct']):
            print(f'    {grp:25s}  {s["n_fungal"]:>5}/{s["n_cells"]:<5} ({s["fungal_pct"]:.1f}%)  RIBO_frac={s["mean_ribo_frac"]:.4f}')

    return results
