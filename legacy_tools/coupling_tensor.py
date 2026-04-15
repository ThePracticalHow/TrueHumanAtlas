# coding: utf-8
"""Universal coupling tensor calculator for h5ad files."""
import scanpy as sc
import numpy as np
import json, sys, time
from collections import Counter
from scipy.stats import spearmanr

def compute_tensor(adata, group_key=None, min_cells=100):
    genes = list(adata.var_names)
    ribo_mask = np.array([g.startswith('RPS') or g.startswith('RPL') for g in genes])
    mito_mask = np.array([g.startswith('MT-') or g.startswith('mt-') for g in genes])
    golgi_prefixes = ('GOLGA', 'GOLGB', 'SEC61', 'COPA', 'COPB', 'MAN1', 'MGAT', 'GORASP')
    golgi_mask = np.array([any(g.startswith(p) for p in golgi_prefixes) for g in genes])
    print(f'  Operator genes: RIBO={ribo_mask.sum()}, MITO={mito_mask.sum()}, GOLGI={golgi_mask.sum()}')
    X = adata.X
    if hasattr(X, 'toarray'):
        ribo = np.asarray(X[:, ribo_mask].sum(axis=1)).flatten()
        mito = np.asarray(X[:, mito_mask].sum(axis=1)).flatten()
        golgi = np.asarray(X[:, golgi_mask].sum(axis=1)).flatten()
        total = np.asarray(X.sum(axis=1)).flatten()
    else:
        ribo = X[:, ribo_mask].sum(axis=1)
        mito = X[:, mito_mask].sum(axis=1)
        golgi = X[:, golgi_mask].sum(axis=1)
        total = X.sum(axis=1)
    nuclear = total - ribo - mito - golgi
    ops = [ribo, mito, nuclear, golgi]
    if group_key and group_key in adata.obs.columns:
        groups = adata.obs[group_key].astype(str).values
    else:
        for key in ['condition', 'cell_type', 'leiden', 'louvain', 'batch', 'sample', 'group', 'stage', 'Phase', 'phase']:
            if key in adata.obs.columns:
                groups = adata.obs[key].astype(str).values
                group_key = key
                break
        else:
            groups = np.array(['ALL'] * adata.n_obs)
            group_key = 'ALL'
    print(f'  Grouping by: {group_key}')
    gcounts = Counter(groups)
    for g, n in gcounts.most_common(20):
        print(f'    {g}: {n:,}')
    results = {'group_key': group_key, 'n_cells_total': int(adata.n_obs), 'n_genes': int(adata.n_vars), 'groups': {}}
    for grp in sorted(gcounts.keys()):
        n = gcounts[grp]
        if n < min_cells:
            continue
        mask = groups == grp
        K = np.zeros((4, 4))
        for i in range(4):
            for j in range(4):
                if i == j:
                    K[i, j] = 1.0
                else:
                    vi = ops[i][mask]
                    vj = ops[j][mask]
                    if vi.std() > 0 and vj.std() > 0:
                        rho, _ = spearmanr(vi, vj)
                        K[i, j] = rho if not np.isnan(rho) else 0
        ribo_coupling = (abs(K[0,1]) + abs(K[0,2]) + abs(K[0,3])) / 3
        det_K = float(np.linalg.det(K))
        ribo_indep = float(1 - ribo_coupling)
        results['groups'][grp] = {
            'n_cells': int(n), 'det_K': det_K, 'RIBO_independence': ribo_indep,
            'K_RM': float(K[0,1]), 'K_RN': float(K[0,2]), 'K_RG': float(K[0,3]),
            'K_MN': float(K[1,2]), 'K_MG': float(K[1,3]), 'K_NG': float(K[2,3]),
        }
        print(f'  {grp} ({n:,}): det={det_K:.6f} RIBO_indep={ribo_indep:.3f} K_RG={K[0,3]:.3f}')
    return results

if __name__ == '__main__':
    h5ad_path = sys.argv[1]
    output_path = sys.argv[2] if len(sys.argv) > 2 else h5ad_path.replace('.h5ad', '_tensor.json')
    group_key = sys.argv[3] if len(sys.argv) > 3 else None
    t0 = time.time()
    print(f'Loading {h5ad_path}...')
    adata = sc.read_h5ad(h5ad_path)
    print(f'  Shape: {adata.shape[0]:,} cells x {adata.shape[1]:,} genes')
    print(f'  obs columns: {list(adata.obs.columns)[:20]}')
    results = compute_tensor(adata, group_key=group_key)
    results['source'] = h5ad_path
    results['compute_time'] = time.time() - t0
    with open(output_path, 'w') as f:
        json.dump(results, f, indent=2)
    print(f'Saved: {output_path} ({time.time()-t0:.1f}s)')
    print('DONE')
