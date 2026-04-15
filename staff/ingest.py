"""
Ingest router — detect file type, run the right tools, adapt, insert into atlas.

    staff ingest sample.bam                     # Track A: all BAM tools
    staff ingest sample.h5ad                    # Track B: coupling tensor
    staff ingest results.json                   # Direct JSON ingest
    staff ingest folder/                        # Everything in the folder
"""
import os
import json
import glob
from pathlib import Path

from .atlas.db import open_atlas, insert_adapted_rows, summary, default_db_path
from . import adapters


def ingest_h5ad(h5ad_path, db_path=None, meta=None):
    """Track B: compute coupling tensor from h5ad, adapt, insert."""
    from .expr.coupling_tensor import run as run_tensor

    print(f"\n  [TRACK B] Expression tensor: {h5ad_path}")
    meta = meta or {}

    result = run_tensor(h5ad_path, group_key=meta.get('group_key'))
    adapted = adapters.adapt_tensor_json(result, h5ad_path, meta)

    conn = open_atlas(db_path)
    insert_adapted_rows(conn, adapted)
    n_groups = len(result.get('groups', {}))
    print(f"\n  Ingested {n_groups} groups into atlas")
    summary(conn)
    conn.close()
    return result


def ingest_bam(bam_path, db_path=None, meta=None):
    """Track A: run all BAM tools, adapt, insert."""
    from .bam.binary_transcripter import run_and_save as run_binary
    from .bam.base4 import run_and_save as run_base4
    from .bam.reading_frame import run_and_save as run_frame
    from .bam.splice_history import run_and_save as run_splice

    print(f"\n  [TRACK A] BAM analysis: {bam_path}")
    meta = meta or {}
    out_dir = os.path.dirname(bam_path) or '.'
    basename = os.path.basename(bam_path).replace('.bam', '')

    print(f"\n  [1/4] Binary transcripter...")
    binary_result = run_binary(bam_path, os.path.join(out_dir, f'{basename}_binary.json'))

    print(f"\n  [2/4] Base-4 XOR...")
    base4_result = run_base4(bam_path, os.path.join(out_dir, f'{basename}_base4.json'))

    print(f"\n  [3/4] Reading frame...")
    frame_result = run_frame(bam_path, os.path.join(out_dir, f'{basename}_frame.json'))

    print(f"\n  [4/4] Splice history...")
    splice_result = run_splice(bam_path, os.path.join(out_dir, f'{basename}_splice.json'))

    # Adapt and insert
    conn = open_atlas(db_path)
    cell_rows = adapters.adapt_bam_binary(binary_result, bam_path, meta)
    insert_adapted_rows(conn, cell_rows)

    cid = cell_rows[0][1]['id']  # cell_state_id from the first row
    insert_adapted_rows(conn, adapters.adapt_bam_base4(base4_result, cid))
    insert_adapted_rows(conn, adapters.adapt_bam_frame(frame_result, cid))
    insert_adapted_rows(conn, adapters.adapt_bam_splice_history(splice_result, cid))

    print(f"\n  Ingested BAM data for {basename} into atlas")
    summary(conn)
    conn.close()

    return {
        'binary': binary_result,
        'base4': base4_result,
        'frame': frame_result,
        'splice': splice_result,
    }


def ingest_json(json_path, db_path=None, meta=None):
    """Ingest a pre-computed JSON result (tensor or BAM output)."""
    print(f"\n  [JSON] Direct ingest: {json_path}")
    with open(json_path) as f:
        data = json.load(f)

    conn = open_atlas(db_path)

    # Detect JSON type and adapt
    if 'groups' in data and 'group_key' in data:
        # Coupling tensor output
        adapted = adapters.adapt_tensor_json(data, json_path, meta)
        insert_adapted_rows(conn, adapted)
        print(f"  Ingested tensor JSON ({len(data['groups'])} groups)")
    elif isinstance(data, list):
        # Census-style flat records
        for record in data:
            adapted = adapters.adapt_census_record(record)
            insert_adapted_rows(conn, adapted)
        print(f"  Ingested {len(data)} Census records")
    elif 'conditions' in data:
        # Mycelium/Hongyan atlas format -- detect which atlas by keys
        atlas_name = data.get('atlas', data.get('name', ''))
        is_mycelium = 'mycelium' in atlas_name.lower() or 'mycelium' in json_path.lower()
        is_hongyan = 'hongyan' in atlas_name.lower() or 'hongyan' in json_path.lower()
        if is_mycelium:
            source_tag = 'mycelium_atlas'
        elif is_hongyan:
            source_tag = 'hongyan_atlas'
        else:
            source_tag = 'json_ingest'

        conditions = data['conditions']
        for name, vals in conditions.items():
            cid = f"{source_tag}_{name}".replace(' ', '_').lower()

            species = 'unknown'
            if is_mycelium:
                if 'yeast' in name: species = 'Saccharomyces cerevisiae'
                elif 'crypto' in name: species = 'Cryptococcus neoformans'
                elif 'candida' in name: species = 'Candida albicans'
                elif 'aspergillus' in name: species = 'Aspergillus fumigatus'
            elif is_hongyan:
                species = 'Homo sapiens (viral infection)'

            conn.execute("""INSERT OR REPLACE INTO cell_states
                (id, name, tissue, disease, source, n_cells, species, category)
                VALUES (?,?,?,?,?,?,?,?)""",
                (cid, name, 'whole_organism' if is_mycelium else 'host_cells',
                 vals.get('virus', 'active_growth') if is_hongyan else 'active_growth',
                 source_tag,
                 vals.get('n_samples', 0), species,
                 'fungal_pathogen' if is_mycelium else ('viral_host_infection' if is_hongyan else 'unknown')))
            conn.execute("""INSERT INTO coupling_tensor
                (cell_state_id, condition, n_cells, ribo_indep, det_k, k_rm, k_rg, k_mn)
                VALUES (?,?,?,?,?,?,?,?)""",
                (cid, name, vals.get('n_samples', 0),
                 vals.get('ribo_indep'), vals.get('det_k'),
                 vals.get('K_RM', 0), vals.get('K_RG', 0), vals.get('K_MN', 0)))
        conn.commit()
        print(f"  Ingested {len(conditions)} conditions from {source_tag}")
    elif isinstance(data, dict) and all(isinstance(v, dict) for v in data.values()):
        # Legacy tensor format: {condition_name: {tensor_values}}
        # or any dict-of-dicts where values have tensor-like keys
        wrapped = {'groups': data, 'group_key': 'condition'}
        adapted = adapters.adapt_tensor_json(wrapped, json_path, meta)
        insert_adapted_rows(conn, adapted)
        print(f"  Ingested {len(data)} conditions from dict-of-dicts JSON")
    else:
        print(f"  WARNING: Unrecognized JSON format in {json_path}")

    summary(conn)
    conn.close()


def ingest_folder(folder_path, db_path=None, meta=None):
    """Ingest everything in a folder."""
    folder = Path(folder_path)
    print(f"\n  [FOLDER] Batch ingest: {folder}")

    files = sorted(folder.iterdir())
    bams = [f for f in files if f.suffix == '.bam']
    h5ads = [f for f in files if f.suffix == '.h5ad']
    jsons = [f for f in files if f.suffix == '.json' and not f.name.startswith('.')]

    print(f"  Found: {len(bams)} BAM, {len(h5ads)} h5ad, {len(jsons)} JSON")

    for f in h5ads:
        ingest_h5ad(str(f), db_path, meta)
    for f in bams:
        ingest_bam(str(f), db_path, meta)
    for f in jsons:
        ingest_json(str(f), db_path, meta)


def ingest(path, db_path=None, meta=None):
    """Auto-detect and ingest. The main entry point."""
    p = Path(path)

    if p.is_dir():
        return ingest_folder(str(p), db_path, meta)
    elif p.suffix == '.h5ad':
        return ingest_h5ad(str(p), db_path, meta)
    elif p.suffix == '.bam':
        return ingest_bam(str(p), db_path, meta)
    elif p.suffix == '.json':
        return ingest_json(str(p), db_path, meta)
    elif str(p).upper().startswith('GSE'):
        print(f"\n  [TRACK C] GEO download not yet implemented: {p}")
        print(f"  Download the h5ad manually and run: staff ingest <file>.h5ad")
        return None
    else:
        print(f"\n  ERROR: Don't know how to ingest: {p}")
        print(f"  Supported: .bam, .h5ad, .json, folder/, GSE*")
        return None
