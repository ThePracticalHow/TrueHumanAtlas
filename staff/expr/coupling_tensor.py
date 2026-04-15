"""
STAFF Expr: Coupling Tensor Calculator
Computes the NxN operator coupling tensor from h5ad single-cell data.

Default operators (5x5):
  RIBO  -- ribosomal (translation)
  MITO  -- mitochondrial (energy)
  NUC   -- nuclear (everything else)
  GOLGI -- Golgi/secretory (communication)
  LYSO  -- lysosomal (the domesticated spore, the recycler)

Backward compatible: if no LYSO genes found, falls back to 4x4.
The tensor dimension adapts to the data.
"""
import numpy as np
import json, sys, time
from collections import Counter
from scipy.stats import spearmanr


# Operator definitions -- order matters (index in the K matrix)
OPERATORS = [
    {
        'name': 'RIBO', 'label': 'R',
        'prefixes': ('RPS', 'RPL'),
        'desc': 'ribosomal (translation)',
    },
    {
        'name': 'MITO', 'label': 'M',
        'prefixes': ('MT-', 'mt-'),
        'desc': 'mitochondrial (energy)',
    },
    {
        'name': 'GOLGI', 'label': 'G',
        'prefixes': ('GOLGA', 'GOLGB', 'SEC61', 'COPA', 'COPB', 'MAN1', 'MGAT', 'GORASP'),
        'desc': 'Golgi/secretory (communication)',
    },
    {
        'name': 'LYSO', 'label': 'L',
        'prefixes': (
            'LAMP1', 'LAMP2', 'LAMP3',
            'CTSA', 'CTSB', 'CTSC', 'CTSD', 'CTSE', 'CTSF',
            'CTSG', 'CTSH', 'CTSK', 'CTSL', 'CTSS', 'CTSZ',
            'ATP6V0', 'ATP6V1',
            'HEXA', 'HEXB',
            'GBA', 'GAA', 'GLB1',
            'NPC1', 'NPC2',
            'TFEB',
        ),
        'desc': 'lysosomal (the domesticated spore)',
    },
]

AUTO_GROUP_KEYS = [
    'condition', 'cell_type', 'leiden', 'louvain', 'batch',
    'sample', 'group', 'stage', 'Phase', 'phase',
]


def _classify_genes(genes, operators):
    """Build boolean masks for each operator. Returns list of masks + NUC (residual)."""
    masks = []
    any_assigned = np.zeros(len(genes), dtype=bool)
    active_ops = []

    for op in operators:
        mask = np.array([any(g.startswith(p) for p in op['prefixes']) for g in genes])
        n = mask.sum()
        if n >= 3:
            masks.append(mask)
            any_assigned |= mask
            active_ops.append(op)

    return masks, any_assigned, active_ops


def compute_tensor(adata, group_key=None, min_cells=100):
    """Compute coupling tensor per group. Adapts to available operators."""
    genes = list(adata.var_names)
    masks, any_assigned, active_ops = _classify_genes(genes, OPERATORS)

    # NUC is always the residual
    op_names = [op['name'] for op in active_ops] + ['NUC']
    op_labels = [op['label'] for op in active_ops] + ['N']
    n_ops = len(op_names)

    gene_counts = {op['name']: int(m.sum()) for op, m in zip(active_ops, masks)}
    gene_counts['NUC'] = int((~any_assigned).sum())
    print(f'  Operators ({n_ops}): ' + ', '.join(f'{n}={gene_counts[n]}' for n in op_names))

    X = adata.X
    sparse = hasattr(X, 'toarray')

    def _sum_col(mask):
        if sparse:
            return np.asarray(X[:, mask].sum(axis=1)).flatten()
        return X[:, mask].sum(axis=1)

    op_vectors = [_sum_col(m) for m in masks]
    total = np.asarray(X.sum(axis=1)).flatten() if sparse else X.sum(axis=1)
    nuc = total - sum(op_vectors)
    op_vectors.append(nuc)

    # Grouping
    if group_key and group_key in adata.obs.columns:
        groups = adata.obs[group_key].astype(str).values
    else:
        for key in AUTO_GROUP_KEYS:
            if key in adata.obs.columns:
                groups = adata.obs[key].astype(str).values
                group_key = key
                break
        else:
            groups = np.array(['ALL'] * adata.n_obs)
            group_key = 'ALL'

    print(f'  Grouping by: {group_key}')
    gcounts = Counter(groups)

    results = {
        'tensor_version': '2.0',
        'operators': op_names,
        'operator_genes': gene_counts,
        'group_key': group_key,
        'n_cells_total': int(adata.n_obs),
        'n_genes': int(adata.n_vars),
        'groups': {},
    }

    for grp in sorted(gcounts.keys()):
        n = gcounts[grp]
        if n < min_cells:
            continue
        cell_mask = groups == grp

        K = np.zeros((n_ops, n_ops))
        for i in range(n_ops):
            for j in range(n_ops):
                if i == j:
                    K[i, j] = 1.0
                else:
                    vi = op_vectors[i][cell_mask]
                    vj = op_vectors[j][cell_mask]
                    if vi.std() > 0 and vj.std() > 0:
                        rho, _ = spearmanr(vi, vj)
                        K[i, j] = rho if not np.isnan(rho) else 0

        ribo_idx = op_names.index('RIBO') if 'RIBO' in op_names else 0
        ribo_others = [abs(K[ribo_idx, j]) for j in range(n_ops) if j != ribo_idx]
        ribo_indep = float(1 - np.mean(ribo_others)) if ribo_others else 0.0
        det_K = float(np.linalg.det(K))

        grp_result = {
            'n_cells': int(n),
            'det_K': det_K,
            'RIBO_independence': ribo_indep,
        }

        # Store all pairwise couplings with meaningful names: K_RM, K_RG, K_RL, etc.
        for i in range(n_ops):
            for j in range(i + 1, n_ops):
                key = f'K_{op_labels[i]}{op_labels[j]}'
                grp_result[key] = float(K[i, j])

        results['groups'][grp] = grp_result

        lyso_str = ''
        if 'LYSO' in op_names:
            li = op_names.index('LYSO')
            gi = op_names.index('GOLGI') if 'GOLGI' in op_names else -1
            k_gl = K[gi, li] if gi >= 0 else 0
            lyso_str = f' K_GL={k_gl:.3f}'

        print(f'  {grp} ({n:,}): det={det_K:.6f} R_ind={ribo_indep:.3f} K_RG={grp_result.get("K_RG", 0):.3f}{lyso_str}')

    return results


def run(h5ad_path, output_path=None, group_key=None):
    """Load h5ad and compute tensor. Returns result dict."""
    import scanpy as sc
    t0 = time.time()
    print(f'Loading {h5ad_path}...')
    adata = sc.read_h5ad(h5ad_path)
    print(f'  Shape: {adata.shape[0]:,} cells x {adata.shape[1]:,} genes')
    print(f'  obs columns: {list(adata.obs.columns)[:20]}')

    results = compute_tensor(adata, group_key=group_key)
    results['source'] = str(h5ad_path)
    results['compute_time'] = time.time() - t0

    if output_path:
        with open(output_path, 'w') as f:
            json.dump(results, f, indent=2)
        print(f'Saved: {output_path} ({time.time() - t0:.1f}s)')

    return results
