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

# ── Geo-temporal provenance ──────────────────────────────────────────────
# Keys: (dataset_prefix, source) or just dataset_prefix.
# Values: (lat, lon, country, collection_year)
# Sources: GEO metadata, CellxGene Census metadata, Stefano email S14,
#          strain databases, publication dates.

GEO_PROVENANCE = {
    # Korean NSCLC cohort — Samsung Medical Center, Seoul
    'GSE131907':    (37.57, 126.98, 'South Korea', 2019),
    # GBM — MGH Boston + Weizmann Tel Aviv (primary site: Boston)
    'GSE131928':    (42.36, -71.07, 'USA', 2019),
    # SGNex — Nanopore, Genome Institute of Singapore
    'SGNex':        (1.29, 103.85, 'Singapore', 2020),
    # WI-38 senescence — NIA/NIH, Baltimore
    'GSE226225':    (39.30, -76.59, 'USA', 2023),
    'GSE226225_SRA':(39.30, -76.59, 'USA', 2023),
    # Proliferative/Senescent — USA
    'GSE250041':    (38.90, -77.04, 'USA', 2024),
    # Normal adult lung — Stanford
    'GSE150247':    (37.43, -122.17, 'USA', 2020),
    # HCV host — GEO (Japan/international)
    'GSE84346':     (35.68, 139.69, 'Japan', 2016),
    # Influenza A host — Duke University
    'GSE97672':     (35.99, -78.90, 'USA', 2017),
    # SARS-CoV-2 — Icahn School of Medicine, NYC
    'blanco_melo':  (40.79, -73.95, 'USA', 2020),
    # Cryptococcus neoformans — GSE67602 — Duke University
    'GSE67602':     (35.99, -78.90, 'USA', 2015),
    # Cryptococcus neoformans — other strains (environmental, global)
    'crypto_atlas': (35.99, -78.90, 'USA', 2015),
    # Saccharomyces cerevisiae reference — GSE144820 — Stanford/UCSF
    'GSE144820':    (37.43, -122.17, 'USA', 2020),
    # Candida albicans — GSE132030 — Brown University
    'GSE132030':    (41.83, -71.40, 'USA', 2019),
    'GEO_various':  (41.83, -71.40, 'USA', 2019),
    # Inner Engineering meditation — GSE174083 — University of Florida
    'GSE174083':    (29.65, -82.34, 'USA', 2021),
}

# Source-level defaults (when dataset field is empty or UUID-only)
SOURCE_PROVENANCE = {
    'census':         (37.77, -122.42, 'USA (CZI Census)', 2023),
    'mycelium_atlas': (35.99, -78.90, 'USA (multi-lab)', 2015),
    'hongyan_atlas':  (40.79, -73.95, 'USA (multi-lab)', 2020),
}

# Staff ingest: derive from the dataset filename
STAFF_PROVENANCE = {
    'd01_celltype_tensor.json': (37.57, 126.98, 'South Korea', 2019),
    'd01_tensor.json':          (37.57, 126.98, 'South Korea', 2019),
    'd01_lyso_results.json':    (37.57, 126.98, 'South Korea', 2019),
    'd04_tensor.json':          (37.77, -122.42, 'USA', 2023),
    'korea_vs_us_results.json': (37.57, 126.98, 'South Korea / USA', 2023),
    'gorospe_lyso_results.json':(39.30, -76.59, 'USA (NIA/NIH)', 2024),
    'macrophage_mouse_results.json': (42.36, -71.07, 'USA', 2023),
}

# Per-record overrides for specific IDs (ethnic/geographic splits)
ID_PROVENANCE = {
    'tensor_korea_vs_us_results_african_american_blood': (33.75, -84.39, 'USA (Atlanta)', 2023),
    'tensor_korea_vs_us_results_african_american_lung':  (33.75, -84.39, 'USA (Atlanta)', 2023),
    'tensor_korea_vs_us_results_japanese_blood':         (35.68, 139.69, 'Japan (Tokyo)', 2023),
    'tensor_korea_vs_us_results_korean_blood':           (37.57, 126.98, 'South Korea (Seoul)', 2023),
}


CATEGORY_FIXES = {
    'd01_celltype_tensor.json': ('immune_celltype', 'Homo sapiens', 'blood'),
    'd01_tensor.json':          (None, 'Homo sapiens', None),
    'd01_lyso_results.json':    ('coculture', 'Homo sapiens', None),
    'd04_tensor.json':          ('census_sweep', 'Homo sapiens', 'vasculature'),
    'korea_vs_us_results.json': ('population_comparison', 'Homo sapiens', None),
    'gorospe_lyso_results.json':('senolytic_experiment', 'Homo sapiens', 'fibroblast'),
    'macrophage_mouse_results.json': ('aging', 'Mus musculus', 'macrophage'),
}

CONDITION_CATEGORY = {
    'Proliferative': 'proliferative',
    'Senescent': 'senescent',
}


def fix_unknowns(records):
    """Fix records with category='unknown' using dataset filename heuristics."""
    fixed = 0
    for rec in records:
        if rec.get('category') != 'unknown':
            continue
        ds = rec.get('dataset', '') or ''
        cond = rec.get('condition', '') or ''
        name = rec.get('name', '') or ''

        if ds in CATEGORY_FIXES:
            cat, sp, tissue = CATEGORY_FIXES[ds]
            if cat:
                rec['category'] = cat
            elif cond in CONDITION_CATEGORY:
                rec['category'] = CONDITION_CATEGORY[cond]
            if sp:
                rec['species'] = sp
            if tissue:
                rec['tissue'] = tissue
            fixed += 1

        rid_lower = (rec.get('id', '') or '').lower()
        if 'korean' in rid_lower:
            rec['tissue'] = 'blood'
        elif 'japanese' in rid_lower:
            rec['tissue'] = 'blood'
        elif 'african_american' in rid_lower:
            rec['tissue'] = 'blood'

        if '_old' in name.lower() or '_Old' in cond:
            rec['age_range'] = 'old'
        elif '_young' in name.lower() or '_Young' in cond:
            rec['age_range'] = 'young'

    # Second pass: force-correct species and tissue for ALL records by dataset
    for rec in records:
        ds = rec.get('dataset', '') or ''
        if ds in CATEGORY_FIXES:
            cat, sp, tissue = CATEGORY_FIXES[ds]
            if sp:
                rec['species'] = sp
            if tissue and not rec.get('tissue'):
                rec['tissue'] = tissue

    return fixed


def stamp_geo_temporal(records):
    """Stamp lat, lon, country, collection_year onto every record."""
    tagged = 0
    for rec in records:
        rid = rec.get('id', rec.get('cell_state_id', ''))
        ds = rec.get('dataset', '') or ''
        src = rec.get('source', '') or ''

        geo = None

        # Priority 1: exact record ID override
        if rid in ID_PROVENANCE:
            geo = ID_PROVENANCE[rid]

        # Priority 2: GEO accession match (check if dataset starts with known key)
        if geo is None:
            for prefix, prov in GEO_PROVENANCE.items():
                if prefix in ds or prefix in rid:
                    geo = prov
                    break

        # Priority 3: staff ingest dataset filename
        if geo is None and ds in STAFF_PROVENANCE:
            geo = STAFF_PROVENANCE[ds]

        # Priority 4: source-level default
        if geo is None and src in SOURCE_PROVENANCE:
            geo = SOURCE_PROVENANCE[src]

        if geo:
            rec['lat'] = geo[0]
            rec['lon'] = geo[1]
            rec['country'] = geo[2]
            rec['collection_year'] = geo[3]
            tagged += 1
        else:
            rec.setdefault('lat', None)
            rec.setdefault('lon', None)
            rec.setdefault('country', None)
            rec.setdefault('collection_year', None)

    return tagged


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

    print('  fixing unknown categories...')
    n_fixed = fix_unknowns(records)
    print(f'  {n_fixed} records re-categorized')

    print('  stamping geo-temporal provenance...')
    n_geo = stamp_geo_temporal(records)
    print(f'  {n_geo}/{len(records)} records geo-tagged')

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
