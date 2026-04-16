#!/usr/bin/env python3
"""
compute_umap.py — WORD protocol enrichment
============================================

Enriches true_human_atlas.json IN-PLACE.
No derived files. The working object IS the published object.

Steps:
  1. Read true_human_atlas.json
  2. Fuse cell_states + coupling_tensor into flat records[]
  3. Compute eigenvalue spectrum per record
  4. Compute UMAP embedding (8D tensor → 2D)
  5. Write back to true_human_atlas.json

Usage:
  python compute_umap.py                          # default
  python compute_umap.py --file data/true_human_atlas.json
"""
import json
import argparse
import numpy as np
from datetime import date

FEATURE_COLS = ['k_rm', 'k_rn', 'k_rg', 'k_mn', 'k_mg', 'k_ng', 'det_k', 'ribo_indep']


def load(path):
    with open(path, encoding='utf-8') as f:
        return json.load(f)


def fuse(atlas):
    """Merge cell_states + coupling_tensor into flat records.
    If records[] already exists, return as-is (already fused)."""
    if 'records' in atlas and isinstance(atlas['records'], list):
        return atlas['records']

    states_by_id = {}
    for s in atlas.get('cell_states', []):
        sid = s.get('id', '')
        if sid:
            states_by_id[sid] = dict(s)

    tensors = atlas.get('coupling_tensor', [])
    seen = set()
    records = []

    for t in tensors:
        sid = t.get('cell_state_id', '')
        key = (sid, t.get('condition', ''))
        if key in seen:
            continue
        seen.add(key)
        rec = {}
        if sid in states_by_id:
            rec.update(states_by_id[sid])
        rec.update(t)
        records.append(rec)

    matched_ids = set(r.get('cell_state_id') or r.get('id') for r in records)
    for sid, s in states_by_id.items():
        if sid not in matched_ids:
            rec = dict(s)
            rec.setdefault('cell_state_id', sid)
            records.append(rec)

    return records


def compute_eigenvalues(rec):
    K = np.eye(4)
    pairs = {
        (0, 1): rec.get('k_rm', 0) or 0,
        (0, 2): rec.get('k_rg', 0) or 0,
        (0, 3): rec.get('k_rn', 0) or 0,
        (1, 2): rec.get('k_mg', 0) or 0,
        (1, 3): rec.get('k_mn', 0) or 0,
        (2, 3): rec.get('k_ng', 0) or 0,
    }
    for (i, j), val in pairs.items():
        K[i, j] = val
        K[j, i] = val

    evals = np.sort(np.linalg.eigvalsh(K))[::-1]
    return {
        'eigenvalue_1': round(float(evals[0]), 6),
        'eigenvalue_2': round(float(evals[1]), 6),
        'eigenvalue_3': round(float(evals[2]), 6),
        'eigenvalue_4': round(float(evals[3]), 6),
        'has_confinement': bool(evals[-1] < 0),
        'confinement_ratio': round(float(evals[0] / abs(evals[-1])), 3) if abs(evals[-1]) > 1e-10 else None,
    }


def build_features(records):
    X = []
    valid = []
    for rec in records:
        row = []
        ok = True
        for col in FEATURE_COLS:
            val = rec.get(col, 0)
            if val is None:
                val = 0
            row.append(float(val))
            if val == 0 and col in ('k_rm', 'k_rg'):
                ok = False
        X.append(row)
        valid.append(ok)
    return np.array(X), valid


def run_umap(X, valid, n_neighbors=15, min_dist=0.3, seed=42):
    import umap
    idx = [i for i, v in enumerate(valid) if v]
    Xv = X[idx]
    if len(Xv) < 5:
        return np.full((len(X), 2), np.nan)

    nn = min(n_neighbors, len(Xv) - 1)
    emb_v = umap.UMAP(
        n_neighbors=nn, min_dist=min_dist,
        n_components=2, metric='euclidean', random_state=seed,
    ).fit_transform(Xv)

    emb = np.full((len(X), 2), np.nan)
    for i, ix in enumerate(idx):
        emb[ix] = emb_v[i]
    return emb


def catalog_fields(records):
    all_keys = set()
    for r in records:
        all_keys.update(r.keys())

    numeric, categorical = [], []
    for key in sorted(all_keys):
        vals = [r.get(key) for r in records if r.get(key) is not None]
        if not vals:
            continue
        if all(isinstance(v, (int, float, bool)) for v in vals[:20]):
            numeric.append(key)
        elif all(isinstance(v, str) for v in vals[:20]):
            categorical.append(key)
    return numeric, categorical


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--file', default='data/true_human_atlas.json')
    args = parser.parse_args()

    print(f'WORD: enriching {args.file} in-place')
    atlas = load(args.file)

    print('  fusing cell_states + coupling_tensor -> records[]...')
    records = fuse(atlas)
    print(f'  {len(records)} records')

    print('  computing eigenvalues...')
    for rec in records:
        has_tensor = any(rec.get(c) for c in FEATURE_COLS[:6])
        if has_tensor:
            rec.update(compute_eigenvalues(rec))

    print('  building features...')
    X, valid = build_features(records)
    n_valid = sum(valid)
    print(f'  {n_valid} with tensor, {len(records) - n_valid} without')

    print('  running UMAP...')
    emb = run_umap(X, valid)
    for i, rec in enumerate(records):
        if not np.isnan(emb[i, 0]):
            rec['umap_x'] = round(float(emb[i, 0]), 4)
            rec['umap_y'] = round(float(emb[i, 1]), 4)
        else:
            rec['umap_x'] = None
            rec['umap_y'] = None

    numeric, categorical = catalog_fields(records)

    output = {
        'version': atlas.get('version', '4.0'),
        'protocol': 'WORD',
        'created': atlas.get('created', ''),
        'updated': str(date.today()),
        'description': atlas.get('description', ''),
        'thesis': atlas.get('thesis', ''),
        'license_code': atlas.get('license_code', 'MIT'),
        'license_data': atlas.get('license_data', 'CC-BY-4.0'),
        'species': atlas.get('species', []),
        'n_records': len(records),
        'n_embedded': sum(1 for r in records if r.get('umap_x') is not None),
        'feature_cols': FEATURE_COLS,
        'numeric_fields': numeric,
        'categorical_fields': categorical,
        'records': records,
    }

    print(f'  writing {args.file}...')
    with open(args.file, 'w', encoding='utf-8') as f:
        json.dump(output, f, indent=2, ensure_ascii=False, default=str)

    kb = len(json.dumps(output, default=str)) / 1024
    print(f'  {len(records)} records, {n_valid} embedded, {kb:.0f} KB')
    print('  WORD: working object = published object. Push to deploy.')


if __name__ == '__main__':
    main()
