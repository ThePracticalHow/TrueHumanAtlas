"""
Atlas database layer — schema, init, insert, query primitives.
Extracted from build_atlas.py. This is the canonical schema.
"""
import sqlite3
import os
from pathlib import Path

SCHEMA = """
CREATE TABLE IF NOT EXISTS cell_states (
    id TEXT PRIMARY KEY,
    name TEXT NOT NULL,
    tissue TEXT,
    disease TEXT DEFAULT 'normal',
    sex TEXT DEFAULT 'mixed',
    age_range TEXT,
    category TEXT,
    species TEXT DEFAULT 'Homo sapiens',
    dataset TEXT,
    source TEXT,
    n_cells INTEGER DEFAULT 0,
    n_reads INTEGER DEFAULT 0,
    platform TEXT
);

CREATE TABLE IF NOT EXISTS coupling_tensor (
    cell_state_id TEXT NOT NULL,
    condition TEXT DEFAULT 'ALL',
    n_cells INTEGER,
    ribo_indep REAL,
    det_k REAL,
    k_rm REAL, k_rn REAL, k_rg REAL,
    k_mn REAL, k_mg REAL, k_ng REAL,
    k_rl REAL, k_ml REAL, k_nl REAL, k_gl REAL,
    k_json TEXT,
    eigenvalue_min REAL,
    eigenvalue_max REAL,
    compute_host TEXT,
    compute_time REAL,
    timestamp TEXT,
    FOREIGN KEY (cell_state_id) REFERENCES cell_states(id)
);

CREATE TABLE IF NOT EXISTS base4_xor (
    cell_state_id TEXT NOT NULL,
    xor_identity REAL,
    xor_complement REAL,
    xor_transversion REAL,
    xor_transition REAL,
    gc_content REAL,
    n_junctions INTEGER,
    platform TEXT,
    FOREIGN KEY (cell_state_id) REFERENCES cell_states(id)
);

CREATE TABLE IF NOT EXISTS reading_frame (
    cell_state_id TEXT NOT NULL,
    frame_preserving_pct REAL,
    stop_at_split_pct REAL,
    atg_near_junction_pct REAL,
    stop_near_junction_pct REAL,
    lys_before_junction_pct REAL,
    leu_after_junction_pct REAL,
    FOREIGN KEY (cell_state_id) REFERENCES cell_states(id)
);

CREATE TABLE IF NOT EXISTS splice_history (
    cell_state_id TEXT NOT NULL,
    n_spliced INTEGER,
    unique_chains INTEGER,
    chain_uniqueness_pct REAL,
    unique_introns INTEGER,
    total_excisions INTEGER,
    median_intron_size INTEGER,
    mirna_range_pct REAL,
    snorna_range_pct REAL,
    FOREIGN KEY (cell_state_id) REFERENCES cell_states(id)
);

CREATE TABLE IF NOT EXISTS structural_index (
    cell_state_id TEXT NOT NULL,
    intron_name TEXT,
    size_bp INTEGER,
    gc_pct REAL,
    mfe_per_nt REAL,
    paired_pct REAL,
    n_stems INTEGER,
    freq_mhz REAL,
    stiffness_class TEXT
);

CREATE TABLE IF NOT EXISTS meditation (
    timepoint TEXT,
    n_samples INTEGER,
    ribo_indep REAL,
    k_rm REAL,
    k_rg REAL,
    dataset TEXT
);

CREATE TABLE IF NOT EXISTS md_simulations (
    id TEXT PRIMARY KEY,
    description TEXT,
    system TEXT,
    n_steps INTEGER,
    ion_type TEXT,
    result_json TEXT,
    compute_host TEXT,
    timestamp TEXT
);

CREATE TABLE IF NOT EXISTS disease_index (
    canonical_name TEXT NOT NULL,
    cell_state_id TEXT NOT NULL,
    tissue TEXT,
    sex TEXT,
    n_cells INTEGER,
    ribo_indep REAL,
    k_rm REAL,
    k_rg REAL,
    failure_mode TEXT,
    FOREIGN KEY (cell_state_id) REFERENCES cell_states(id)
);

CREATE TABLE IF NOT EXISTS fungal_load (
    cell_state_id TEXT NOT NULL,
    bam_source TEXT,
    n_unmapped_scanned INTEGER,
    n_classified_fungal INTEGER,
    fungal_pct REAL,
    candida_hits INTEGER DEFAULT 0,
    aspergillus_hits INTEGER DEFAULT 0,
    crypto_hits INTEGER DEFAULT 0,
    yeast_hits INTEGER DEFAULT 0,
    n_cells_with_fungus INTEGER DEFAULT 0,
    scan_timestamp TEXT,
    FOREIGN KEY (cell_state_id) REFERENCES cell_states(id)
);

CREATE INDEX IF NOT EXISTS idx_cs_tissue ON cell_states(tissue);
CREATE INDEX IF NOT EXISTS idx_cs_disease ON cell_states(disease);
CREATE INDEX IF NOT EXISTS idx_cs_sex ON cell_states(sex);
CREATE INDEX IF NOT EXISTS idx_ct_ribo ON coupling_tensor(ribo_indep);
"""

# Column lists for each table (for insert)
TABLE_COLS = {
    'cell_states': ['id', 'name', 'tissue', 'disease', 'sex', 'age_range', 'category',
                     'species', 'dataset', 'source', 'n_cells', 'n_reads', 'platform'],
    'coupling_tensor': ['cell_state_id', 'condition', 'n_cells', 'ribo_indep', 'det_k',
                         'k_rm', 'k_rn', 'k_rg', 'k_mn', 'k_mg', 'k_ng',
                         'k_rl', 'k_ml', 'k_nl', 'k_gl',
                         'k_json', 'eigenvalue_min', 'eigenvalue_max',
                         'compute_host', 'compute_time', 'timestamp'],
    'base4_xor': ['cell_state_id', 'xor_identity', 'xor_complement',
                   'xor_transversion', 'xor_transition', 'gc_content', 'n_junctions', 'platform'],
    'reading_frame': ['cell_state_id', 'frame_preserving_pct', 'stop_at_split_pct',
                       'atg_near_junction_pct', 'stop_near_junction_pct',
                       'lys_before_junction_pct', 'leu_after_junction_pct'],
    'splice_history': ['cell_state_id', 'n_spliced', 'unique_chains', 'chain_uniqueness_pct',
                        'unique_introns', 'total_excisions', 'median_intron_size',
                        'mirna_range_pct', 'snorna_range_pct'],
}


def default_db_path():
    """Return the default atlas database path."""
    return Path(__file__).resolve().parent.parent.parent / "true_human_atlas.db"


def open_atlas(db_path=None):
    """Open or create the atlas database."""
    path = str(db_path or default_db_path())
    parent = os.path.dirname(path)
    if parent and not os.path.exists(parent):
        os.makedirs(parent, exist_ok=True)
    conn = sqlite3.connect(path)
    conn.executescript(SCHEMA)
    conn.commit()
    return conn


def insert_row(conn, table, row_dict):
    """Insert a row into a table from a dict. Uses INSERT OR REPLACE for cell_states."""
    cols = TABLE_COLS.get(table)
    if not cols:
        raise ValueError(f"Unknown table: {table}")
    values = [row_dict.get(c) for c in cols]
    placeholders = ','.join(['?'] * len(cols))
    col_names = ','.join(cols)
    verb = "INSERT OR REPLACE" if table == 'cell_states' else "INSERT"
    conn.execute(f"{verb} INTO {col_names.split(',')[0].strip() and table} ({col_names}) VALUES ({placeholders})", values)


def insert_adapted_rows(conn, adapted_rows):
    """Insert a list of (table, row_dict) tuples from an adapter."""
    for table, row_dict in adapted_rows:
        insert_row(conn, table, row_dict)
    conn.commit()


def count_states(conn):
    return conn.execute("SELECT COUNT(*) FROM cell_states").fetchone()[0]


def count_tensors(conn):
    return conn.execute("SELECT COUNT(*) FROM coupling_tensor").fetchone()[0]


def summary(conn):
    """Print a compact summary of the atlas."""
    n_states = count_states(conn)
    n_tensor = count_tensors(conn)
    n_cells = conn.execute("SELECT COALESCE(SUM(n_cells),0) FROM coupling_tensor").fetchone()[0]

    print(f"\n  Atlas: {n_states} cell states, {n_tensor} tensors, {int(n_cells):,} cells")

    for row in conn.execute("SELECT source, COUNT(*) FROM cell_states GROUP BY source ORDER BY COUNT(*) DESC"):
        print(f"    {(row[0] or '?'):20s} {row[1]:4d} states")
