"""
Adapters — normalize ANY tool output JSON into the unified atlas schema.

The atlas has these tables:
  cell_states       (id, name, tissue, disease, sex, category, species, dataset, source, n_cells, platform)
  coupling_tensor   (cell_state_id, condition, n_cells, ribo_indep, det_k, k_rm, k_rn, k_rg, ...)
  base4_xor         (cell_state_id, xor_identity, xor_complement, xor_transversion, xor_transition, ...)
  reading_frame     (cell_state_id, frame_preserving_pct, stop_at_split_pct, ...)
  splice_history    (cell_state_id, n_spliced, unique_chains, chain_uniqueness_pct, ...)

Every adapter returns a list of (table_name, row_dict) tuples ready for insertion.
"""
import os
import time
import socket


def _make_id(source_path, prefix="ingest"):
    """Create a stable cell_state id from a file path."""
    basename = os.path.basename(source_path).replace('.h5ad', '').replace('.bam', '')
    basename = basename.replace('.json', '').replace(' ', '_').lower()
    return f"{prefix}_{basename}"


def _meta_defaults(meta):
    """Fill in metadata defaults."""
    return {
        'tissue': meta.get('tissue', 'unknown'),
        'disease': meta.get('disease', 'normal'),
        'sex': meta.get('sex', 'mixed'),
        'species': meta.get('species', 'Homo sapiens'),
        'category': meta.get('category', 'unknown'),
        'platform': meta.get('platform', 'unknown'),
    }


def adapt_tensor_json(tensor_result, source_path, meta=None):
    """
    Adapt coupling_tensor.py output (nested groups dict) to atlas rows.
    Input: {'group_key': ..., 'groups': {'grp1': {'n_cells': N, 'RIBO_independence': ..., ...}}}
    Output: list of (table, row_dict) tuples.
    """
    meta = _meta_defaults(meta or {})
    rows = []
    groups = tensor_result.get('groups', {})
    base_id = _make_id(source_path, 'tensor')
    host = socket.gethostname()
    ts = time.strftime('%Y-%m-%d %H:%M')

    for grp_name, grp_data in groups.items():
        cid = f"{base_id}_{grp_name}".replace(' ', '_').lower()
        n_cells = grp_data.get('n_cells', 0)

        rows.append(('cell_states', {
            'id': cid,
            'name': f"{os.path.basename(source_path)} / {grp_name}",
            'tissue': meta['tissue'],
            'disease': meta['disease'],
            'sex': meta['sex'],
            'category': meta['category'],
            'species': meta['species'],
            'dataset': source_path,
            'source': 'staff_ingest',
            'n_cells': n_cells,
            'platform': meta['platform'],
        }))

        rows.append(('coupling_tensor', {
            'cell_state_id': cid,
            'condition': grp_name,
            'n_cells': n_cells,
            'ribo_indep': grp_data.get('RIBO_independence', grp_data.get('ribo_indep')),
            'det_k': grp_data.get('det_K', grp_data.get('det_k')),
            'k_rm': grp_data.get('K_RM', grp_data.get('k_rm', 0)),
            'k_rn': grp_data.get('K_RN', grp_data.get('k_rn', 0)),
            'k_rg': grp_data.get('K_RG', grp_data.get('k_rg', 0)),
            'k_mn': grp_data.get('K_MN', grp_data.get('k_mn', 0)),
            'k_mg': grp_data.get('K_MG', grp_data.get('k_mg', 0)),
            'k_ng': grp_data.get('K_NG', grp_data.get('k_ng', 0)),
            'compute_host': host,
            'timestamp': ts,
        }))

    return rows


def adapt_bam_binary(binary_result, source_path, meta=None):
    """Adapt staff_binary_transcripter output to atlas rows."""
    meta = _meta_defaults(meta or {})
    cid = _make_id(source_path, 'bam')
    rows = []

    rows.append(('cell_states', {
        'id': cid,
        'name': f"BAM: {os.path.basename(source_path)}",
        'tissue': meta['tissue'],
        'disease': meta['disease'],
        'sex': meta['sex'],
        'category': meta['category'],
        'species': meta['species'],
        'dataset': source_path,
        'source': 'staff_ingest',
        'n_cells': binary_result.get('n_reads', 0),
        'platform': meta['platform'],
    }))

    return rows


def adapt_bam_base4(base4_result, cell_state_id):
    """Adapt staff_base4 output to atlas base4_xor row."""
    xor_dist = base4_result.get('xor_distribution', {})
    total = sum(xor_dist.values()) or 1
    return [('base4_xor', {
        'cell_state_id': cell_state_id,
        'xor_identity': 100 * xor_dist.get('0000', xor_dist.get(0, 0)) / total,
        'xor_complement': 100 * xor_dist.get('0001', xor_dist.get(1, 0)) / total,
        'xor_transversion': 100 * xor_dist.get('0010', xor_dist.get(2, 0)) / total,
        'xor_transition': 100 * xor_dist.get('0011', xor_dist.get(3, 0)) / total,
        'gc_content': base4_result.get('gc_content', 0),
        'n_junctions': base4_result.get('total_junctions', 0),
        'platform': 'long_read',
    })]


def adapt_bam_frame(frame_result, cell_state_id):
    """Adapt staff_reading_frame output to atlas reading_frame row."""
    return [('reading_frame', {
        'cell_state_id': cell_state_id,
        'frame_preserving_pct': frame_result.get('frame_preserving_pct', 0),
        'stop_at_split_pct': frame_result.get('stop_at_split_pct', 0),
        'atg_near_junction_pct': frame_result.get('atg_near_junction_pct', 0),
        'stop_near_junction_pct': frame_result.get('stop_near_junction_pct', 0),
        'lys_before_junction_pct': frame_result.get('lys_before_junction_pct', 0),
        'leu_after_junction_pct': frame_result.get('leu_after_junction_pct', 0),
    })]


def adapt_bam_splice_history(splice_result, cell_state_id):
    """Adapt staff_splice_history output to atlas splice_history row."""
    return [('splice_history', {
        'cell_state_id': cell_state_id,
        'n_spliced': splice_result.get('n_spliced', 0),
        'unique_chains': splice_result.get('unique_chains', 0),
        'chain_uniqueness_pct': splice_result.get('chain_uniqueness_pct', 0),
        'unique_introns': splice_result.get('unique_introns', 0),
        'total_excisions': splice_result.get('total_excisions', 0),
        'median_intron_size': splice_result.get('junction_size_median', splice_result.get('median_intron_size', 0)),
        'mirna_range_pct': splice_result.get('mirna_range_pct', 0),
        'snorna_range_pct': splice_result.get('snorna_range_pct', 0),
    })]


def adapt_census_record(record):
    """Adapt a CellxGene Census-style flat record (as used by atlas_results.json)."""
    import json as _json
    tissue = record.get('tissue', 'unknown')
    disease = record.get('disease', 'normal')
    sex = record.get('sex', 'mixed')
    cid = f"census_{tissue}_{disease}_{sex}".replace(' ', '_').lower()

    rows = []
    rows.append(('cell_states', {
        'id': cid,
        'name': f"{tissue} / {disease} / {sex}",
        'tissue': tissue,
        'disease': disease,
        'sex': sex,
        'category': 'census_sweep',
        'species': record.get('species', 'Homo sapiens'),
        'dataset': record.get('dataset_id', ''),
        'source': 'census',
        'n_cells': record.get('n_cells', 0),
        'platform': 'cellxgene_census',
    }))

    k_mn = k_mg = k_ng = 0.0
    K = record.get('K')
    k_json_str = None
    if isinstance(K, list) and len(K) >= 4:
        k_json_str = _json.dumps(K)
        k_mn = K[1][2] if len(K[1]) > 2 else 0
        k_mg = K[1][3] if len(K[1]) > 3 else 0
        k_ng = K[2][3] if len(K[2]) > 3 else 0

    rows.append(('coupling_tensor', {
        'cell_state_id': cid,
        'condition': disease,
        'n_cells': record.get('n_cells', 0),
        'ribo_indep': record.get('RIBO_indep', 0),
        'det_k': record.get('det_K', 0),
        'k_rm': record.get('K_RM', 0),
        'k_rn': record.get('K_RN', 0),
        'k_rg': record.get('K_RG', 0),
        'k_mn': k_mn,
        'k_mg': k_mg,
        'k_ng': k_ng,
        'k_json': k_json_str,
    }))

    return rows
