"""
Atlas query — natural language interface.
Uses the same logic as ask_atlas.py but works with any DB path.
"""
import sqlite3
import re
import sys
from pathlib import Path
from ..atlas.db import open_atlas, default_db_path


FRONT_MAP = {
    'ours':     "cs.source NOT IN ('mycelium_atlas','hongyan_atlas')",
    'human':    "cs.source NOT IN ('mycelium_atlas','hongyan_atlas')",
    'soldiers': "cs.source NOT IN ('mycelium_atlas','hongyan_atlas')",
    'fungus':   "cs.source = 'mycelium_atlas'",
    'fungi':    "cs.source = 'mycelium_atlas'",
    'fungal':   "cs.source = 'mycelium_atlas'",
    'green':    "cs.source = 'mycelium_atlas'",
    'mycelium': "cs.source = 'mycelium_atlas'",
    'yeast':    "cs.source = 'mycelium_atlas' AND cs.species LIKE '%cerevisiae%'",
    'crypto':   "cs.source = 'mycelium_atlas' AND cs.species LIKE '%neoformans%'",
    'virus':    "cs.source = 'hongyan_atlas'",
    'viral':    "cs.source = 'hongyan_atlas'",
    'red':      "cs.source = 'hongyan_atlas'",
    'hongyan':  "cs.source = 'hongyan_atlas'",
}

SEX_MAP = {
    'female': "cs.sex = 'female'", 'woman': "cs.sex = 'female'", 'women': "cs.sex = 'female'",
    'male': "cs.sex = 'male'", 'man': "cs.sex = 'male'", 'men': "cs.sex = 'male'",
}

DISEASE_MAP = {
    'cancer': "cs.disease LIKE '%cancer%' OR cs.disease LIKE '%carcinoma%' OR cs.disease LIKE '%tumor%' OR cs.disease LIKE '%leukemia%' OR cs.category IN ('cancer','metastatic')",
    'alzheimer': "cs.disease LIKE '%Alzheimer%'", 'parkinson': "cs.disease LIKE '%Parkinson%'",
    'covid': "cs.disease LIKE '%COVID%' OR cs.disease LIKE '%SARS%'",
    'hcv': "cs.disease LIKE '%HCV%' OR cs.disease LIKE '%hepatitis%'",
    'gbm': "cs.disease LIKE '%glioblastoma%' OR cs.name LIKE '%GBM%'",
    'lung': "cs.disease LIKE '%lung%' OR cs.tissue LIKE '%lung%'",
    'breast': "cs.disease LIKE '%breast%'",
    'normal': "cs.disease = 'normal'", 'healthy': "cs.disease = 'normal'",
    'nsclc': "cs.name LIKE '%NSCLC%' OR cs.name LIKE '%nsclc%'",
}

TISSUE_MAP = {
    'brain': "cs.tissue LIKE '%brain%' OR cs.tissue LIKE '%cortex%'",
    'lung': "cs.tissue LIKE '%lung%' OR cs.tissue LIKE '%respiratory%'",
    'blood': "cs.tissue LIKE '%blood%' OR cs.tissue LIKE '%bone marrow%'",
    'liver': "cs.tissue LIKE '%liver%'",
}


def _tokens_to_where(tokens):
    clauses = []
    for tok in tokens:
        if tok in FRONT_MAP: clauses.append(FRONT_MAP[tok])
        if tok in DISEASE_MAP: clauses.append(DISEASE_MAP[tok])
        if tok in TISSUE_MAP: clauses.append(TISSUE_MAP[tok])
        if tok in SEX_MAP: clauses.append(SEX_MAP[tok])
    if not clauses:
        terms = [t for t in tokens if len(t) > 2 and t not in ('the','for','and','all','show','get','what','how','me')]
        for term in terms:
            clauses.append(f"(cs.name LIKE '%{term}%' OR cs.disease LIKE '%{term}%' OR cs.tissue LIKE '%{term}%')")
    return clauses


def _front_bar(ribo, source):
    n = int(ribo * 35)
    if source == 'mycelium_atlas': return 'x' * n
    if source == 'hongyan_atlas': return 'V' * n
    return '#' * n


def _front_label(source):
    if source == 'mycelium_atlas': return 'GREEN'
    if source == 'hongyan_atlas': return 'RED'
    return 'OURS'


def _run_query(conn, where_clauses, sort='ct.ribo_indep ASC', limit=15):
    where = " AND ".join(where_clauses) if where_clauses else "1=1"
    sql = f"""SELECT cs.name, cs.species, cs.source, cs.disease, cs.sex,
                     ct.ribo_indep, ct.k_rm, ct.k_rg, ct.det_k, ct.n_cells
              FROM coupling_tensor ct
              JOIN cell_states cs ON ct.cell_state_id = cs.id
              WHERE ct.ribo_indep IS NOT NULL AND ({where})
              ORDER BY {sort} LIMIT {limit}"""
    return conn.execute(sql).fetchall()


def _print_rows(rows, label):
    if not rows:
        print(f"\n  No results for: {label}")
        return
    print(f"\n  {label} ({len(rows)} entries)")
    print(f"  {'R_ind':>6}  {'Front':5}  {'Sex':6}  {'Cells':>9}  {'K_RM':>6}  {'K_RG':>6}  {'Condition'}")
    print(f"  {'-'*95}")
    for name, species, source, disease, sex, ribo, k_rm, k_rg, det_k, n_cells in rows:
        front = _front_label(source)
        bar = _front_bar(ribo, source)
        n_str = f"{n_cells:>9,}" if n_cells else "       --"
        k_rm_v = f"{k_rm:>6.3f}" if k_rm else "    --"
        k_rg_v = f"{k_rg:>6.3f}" if k_rg else "    --"
        nm = (name or '?')[:50]
        print(f"  {ribo:>6.3f}  {front:5}  {(sex or '?')[:6]:6}  {n_str}  {k_rm_v}  {k_rg_v}  {nm}  {bar}")


def ask(query_string, db_path=None):
    """Run a natural language query against the atlas."""
    conn = open_atlas(db_path)
    tokens = re.split(r'[\s,]+', query_string.lower().strip())
    is_compare = 'vs' in tokens or 'compare' in tokens

    # Extract global sex filter
    global_sex = None
    for tok in tokens:
        if tok in SEX_MAP:
            global_sex = SEX_MAP[tok]

    noise = {'compare','show','me','the','a','an','for','in','of','with','get','what','how'}
    sex_words = set(SEX_MAP.keys())

    print(f"\n  Query: \"{query_string}\"")

    if is_compare:
        text = query_string.lower()
        for w in noise | sex_words:
            text = re.sub(r'\b' + re.escape(w) + r'\b', ' ', text)
        text = re.sub(r'\s+', ' ', text).strip()

        groups = []
        for sep in [' vs ', ' versus ', ' compared to ']:
            if sep in text:
                groups = [p.strip() for p in text.split(sep) if p.strip()]
                break

        if len(groups) >= 2:
            print(f"\n{'='*100}")
            print(f"  COMPARISON: {' vs '.join(g.upper() for g in groups)}")
            print(f"{'='*100}")

            all_avgs = {}
            for group in groups:
                toks = re.split(r'[\s,]+', group)
                clauses = _tokens_to_where(toks)
                is_nonhuman = any(t in FRONT_MAP and ('mycelium' in FRONT_MAP.get(t,'') or 'hongyan' in FRONT_MAP.get(t,'')) for t in toks)
                if global_sex and not is_nonhuman:
                    clauses.append(global_sex)
                rows = _run_query(conn, clauses)
                _print_rows(rows, group.upper())
                ribos = [r[5] for r in rows if r[5] is not None]
                if ribos:
                    all_avgs[group] = sum(ribos) / len(ribos)

            if len(all_avgs) >= 2:
                print(f"\n  {'-'*80}")
                print(f"  SUMMARY:")
                for g, avg in sorted(all_avgs.items(), key=lambda x: x[1]):
                    print(f"    {g.upper():20s}  avg_R={avg:.3f}")
            conn.close()
            return

    # Single query
    clauses = _tokens_to_where(tokens)
    if global_sex:
        clauses.append(global_sex)
    rows = _run_query(conn, clauses, limit=30)
    label = " ".join(t for t in tokens if t not in noise)
    _print_rows(rows, label.upper())
    if rows:
        ribos = [r[5] for r in rows if r[5] is not None]
        if ribos:
            print(f"\n  Avg RIBO_indep: {sum(ribos)/len(ribos):.3f}  |  Range: [{min(ribos):.3f} - {max(ribos):.3f}]  |  {len(rows)} entries")
    conn.close()
