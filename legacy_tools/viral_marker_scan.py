# coding: utf-8
"""Scan any h5ad for the viral/trophoblast/imprinting marker panel."""
import scanpy as sc
import numpy as np
import json, sys, time
from collections import Counter

MARKER_PANEL = {
    # Endogenous retroviruses
    'ERV': ['ERVW-1','ERVFRD-1','ERV3-1','ERVK13-1','ERVE-1','ERVK-28','ERVK9-11','ERVK3-1','ERVH48-1','ERVH-1','ERVMER34-1','ERVV-1','ERVV-2'],
    # Paternal imprinted (father's virus program)
    'PATERNAL': ['PEG10','PEG3','MEST','DLK1','SNRPN','MAGEL2','NNAT','PLAGL1','IGF2','NDN'],
    # Maternal imprinted (mother's suppression)
    'MATERNAL': ['H19','MEG3','CDKN1C','GRB10','PHLDA2','SLC22A18','UBE3A','KCNQ1OT1','IGF2R','GNAS'],
    # Trophoblast/placental
    'TROPHOBLAST': ['GCM1','GATA3','TFAP2C','HLA-G','LGALS13','LGALS1','LGALS3','PSG1','PSG2','PSG3','PSG4','PSG5','PSG6','PSG7','PSG8','PSG9'],
    # Immune evasion (shared cancer/trophoblast)
    'IMMUNE_EVASION': ['CD274','PDCD1LG2','IDO1','IDO2','HLA-E','B2M','IL10','TGFB1','FOXP3','CTLA4'],
    # EMT/invasion (trophoblast migration = GBM darting)
    'INVASION': ['SNAI1','SNAI2','ZEB1','ZEB2','TWIST1','TWIST2','VIM','CDH1','CDH2','FN1','MMP2','MMP9','MMP14','MMP7','SPARC','POSTN','LOX','LOXL2','PLAU','PLAUR','SERPINE1','ITGAV','ITGB1','ITGB3'],
    # Antenna lncRNAs
    'ANTENNA': ['LINC01235','LINC02154','LINC01705','NEAT1','MALAT1','H19','MEG3','PVT1','GAS5','HOTAIR','XIST'],
    # Viral receptors / fusogens
    'VIRAL_RECEPTOR': ['SLC1A5','SLC1A4','CXCR4','CCR5','FURIN','TMPRSS2','ACE2','NRP1','MFSD2A'],
    # Sex chromosomes
    'SEX': ['XIST','SRY','DDX3Y','EIF1AY','KDM5D','UTY','USP9Y','ZFY','DDX3X','KDM5C','KDM6A'],
    # Senescence/cancer
    'CANCER_SEN': ['CDKN1A','CDKN2A','TP53','RB1','MYC','GDF15','IL6','CXCL8'],
}

ALL_MARKERS = []
for cat, genes in MARKER_PANEL.items():
    ALL_MARKERS.extend(genes)
ALL_MARKERS = list(set(ALL_MARKERS))

def scan_h5ad(path, group_key=None):
    t0 = time.time()
    print(f'\n{"="*60}')
    print(f'SCANNING: {path}')
    print(f'{"="*60}')

    adata = sc.read_h5ad(path)
    print(f'Shape: {adata.shape[0]:,} cells x {adata.shape[1]:,} genes')
    print(f'Obs columns: {list(adata.obs.columns)[:15]}')

    # Find group key
    if group_key and group_key in adata.obs.columns:
        groups = adata.obs[group_key].astype(str).values
    else:
        for key in ['condition','cell_type','celltype_curated','leiden','louvain','sample','group','stage']:
            if key in adata.obs.columns:
                groups = adata.obs[key].astype(str).values
                group_key = key
                break
        else:
            groups = np.array(['ALL'] * adata.n_obs)
            group_key = 'ALL'

    print(f'Grouping by: {group_key}')
    gcounts = Counter(groups)
    for g, n in gcounts.most_common(10):
        print(f'  {g}: {n:,}')

    # Find present markers
    present = [g for g in ALL_MARKERS if g in adata.var_names]
    missing = [g for g in ALL_MARKERS if g not in adata.var_names]
    print(f'Markers present: {len(present)}/{len(ALL_MARKERS)}')

    # Extract expression
    X = adata[:, present].X
    if hasattr(X, 'toarray'):
        X = X.toarray()

    results = {'source': path, 'group_key': group_key, 'n_cells': int(adata.shape[0]),
               'n_genes': int(adata.shape[1]), 'markers_present': len(present), 'groups': {}}

    # Per-group marker expression
    for grp in sorted(gcounts.keys()):
        n = gcounts[grp]
        if n < 50:
            continue
        mask = groups == grp
        grp_data = {}

        for cat, cat_genes in MARKER_PANEL.items():
            cat_present = [g for g in cat_genes if g in present]
            if not cat_present:
                continue
            for gene in cat_present:
                gi = present.index(gene)
                vals = X[mask, gi]
                mean = float(vals.mean())
                pct = float((vals > 0).mean() * 100)
                if mean > 0.0001 or pct > 0.1:
                    grp_data[gene] = {'mean': mean, 'pct': pct, 'category': cat}

        results['groups'][grp] = {'n_cells': int(n), 'markers': grp_data}

    # Print summary by category
    for cat, cat_genes in MARKER_PANEL.items():
        cat_present = [g for g in cat_genes if g in present]
        if not cat_present:
            continue
        print(f'\n--- {cat} ---')
        header = f'{"Gene":12s}'
        grp_names = [g for g in sorted(gcounts.keys()) if gcounts[g] >= 50][:8]
        for gn in grp_names:
            header += f' {gn[:10]:>10s}'
        print(header)

        for gene in cat_present:
            gi = present.index(gene)
            line = f'{gene:12s}'
            any_expr = False
            for gn in grp_names:
                mask = groups == gn
                mean = X[mask, gi].mean()
                if mean > 0.001:
                    any_expr = True
                line += f' {mean:10.4f}'
            if any_expr:
                print(line)

    print(f'\nScan time: {time.time()-t0:.1f}s')
    return results

if __name__ == '__main__':
    paths = sys.argv[1:]
    all_results = []
    for path in paths:
        try:
            r = scan_h5ad(path)
            all_results.append(r)
        except Exception as e:
            print(f'ERROR on {path}: {e}')

    outpath = '/tmp/viral_marker_scan.json'
    with open(outpath, 'w') as f:
        json.dump(all_results, f, indent=2)
    print(f'\nSaved: {outpath}')
    print('ALL SCANS COMPLETE')
