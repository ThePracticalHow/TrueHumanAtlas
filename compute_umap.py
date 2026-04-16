#!/usr/bin/env python3
"""
compute_umap.py — Unified Atlas UMAP Embedding Generator
=========================================================

Reads true_human_atlas.json, joins cell_states with coupling_tensor records,
computes UMAP from the tensor feature vector, and writes a unified
atlas_umap.json ready for the browser viewer.

Feature vector (8D):
  k_rm, k_rn, k_rg, k_mn, k_mg, k_ng, det_k, ribo_indep

Every field from both cell_states and coupling_tensor is preserved.
New fields added to either source automatically flow through.

Usage:
  python compute_umap.py                     # default paths
  python compute_umap.py --input data/true_human_atlas.json --output data/atlas_umap.json
"""
import json
import argparse
import numpy as np

FEATURE_COLS = ['k_rm', 'k_rn', 'k_rg', 'k_mn', 'k_mg', 'k_ng', 'det_k', 'ribo_indep']


def load_atlas(path):
    with open(path, encoding='utf-8') as f:
        return json.load(f)


def join_states_and_tensors(atlas):
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

    for sid, s in states_by_id.items():
        if not any(r.get('cell_state_id') == sid or r.get('id') == sid for r in records):
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
        'eigenvalue_min': round(float(evals[-1]), 6),
        'eigenvalue_max': round(float(evals[0]), 6),
        'has_confinement': bool(evals[-1] < 0),
        'confinement_ratio': round(float(evals[0] / abs(evals[-1])), 3) if abs(evals[-1]) > 1e-10 else None,
    }


def build_feature_matrix(records):
    X = []
    valid_mask = []
    for rec in records:
        row = []
        valid = True
        for col in FEATURE_COLS:
            val = rec.get(col, 0)
            if val is None:
                val = 0
            row.append(float(val))
            if val == 0 and col in ('k_rm', 'k_rg'):
                valid = False
        X.append(row)
        valid_mask.append(valid)
    return np.array(X), valid_mask


def run_umap(X, valid_mask, n_neighbors=15, min_dist=0.3, random_state=42):
    import umap
    valid_idx = [i for i, v in enumerate(valid_mask) if v]
    X_valid = X[valid_idx]

    if len(X_valid) < 5:
        print(f"  Only {len(X_valid)} valid records, skipping UMAP")
        return np.zeros((len(X), 2))

    n_neighbors_actual = min(n_neighbors, len(X_valid) - 1)
    reducer = umap.UMAP(
        n_neighbors=n_neighbors_actual,
        min_dist=min_dist,
        n_components=2,
        metric='euclidean',
        random_state=random_state,
    )
    embedding_valid = reducer.fit_transform(X_valid)

    embedding = np.full((len(X), 2), np.nan)
    for i, idx in enumerate(valid_idx):
        embedding[idx] = embedding_valid[i]

    return embedding


def main():
    parser = argparse.ArgumentParser(description='Compute UMAP for True Human Atlas')
    parser.add_argument('--input', default='data/true_human_atlas.json')
    parser.add_argument('--output', default='data/atlas_umap.json')
    args = parser.parse_args()

    print(f"Loading {args.input}...")
    atlas = load_atlas(args.input)

    print("Joining cell_states + coupling_tensor...")
    records = join_states_and_tensors(atlas)
    print(f"  {len(records)} unified records")

    print("Computing eigenvalues...")
    for rec in records:
        has_tensor = any(rec.get(c) for c in FEATURE_COLS[:6])
        if has_tensor:
            eig = compute_eigenvalues(rec)
            rec.update(eig)

    print("Building feature matrix...")
    X, valid_mask = build_feature_matrix(records)
    n_valid = sum(valid_mask)
    print(f"  {n_valid} records with valid tensor data, {len(records) - n_valid} without")

    print("Running UMAP...")
    embedding = run_umap(X, valid_mask)

    for i, rec in enumerate(records):
        if not np.isnan(embedding[i, 0]):
            rec['umap_x'] = round(float(embedding[i, 0]), 4)
            rec['umap_y'] = round(float(embedding[i, 1]), 4)
        else:
            rec['umap_x'] = None
            rec['umap_y'] = None

    all_keys = set()
    for rec in records:
        all_keys.update(rec.keys())

    numeric_keys = []
    categorical_keys = []
    for key in sorted(all_keys):
        vals = [rec.get(key) for rec in records if rec.get(key) is not None]
        if not vals:
            continue
        if all(isinstance(v, (int, float)) for v in vals[:20]):
            numeric_keys.append(key)
        elif all(isinstance(v, str) for v in vals[:20]):
            categorical_keys.append(key)

    output = {
        'version': '4.0',
        'created': atlas.get('created', ''),
        'description': atlas.get('description', ''),
        'thesis': atlas.get('thesis', ''),
        'license_code': atlas.get('license_code', 'MIT'),
        'license_data': atlas.get('license_data', 'CC-BY-4.0'),
        'species': atlas.get('species', []),
        'n_records': len(records),
        'n_with_umap': sum(1 for r in records if r.get('umap_x') is not None),
        'feature_cols': FEATURE_COLS,
        'numeric_fields': numeric_keys,
        'categorical_fields': categorical_keys,
        'records': records,
    }

    print(f"Writing {args.output}...")
    with open(args.output, 'w', encoding='utf-8') as f:
        json.dump(output, f, indent=2, ensure_ascii=False, default=str)

    size_kb = len(json.dumps(output, default=str)) / 1024
    print(f"  {len(records)} records, {size_kb:.0f} KB")
    print("Done.")


if __name__ == '__main__':
    main()
