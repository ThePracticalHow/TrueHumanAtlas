"""
Fungal transcript scanner -- extract unmapped reads from BAM, classify against panel.

For every unmapped read:
  1. Decode the sequence
  2. Extract cell barcode (CB tag for 10x, or read name for nanopore)
  3. Classify against fungal k-mer panel
  4. Count hits per species per cell

Pure Python BAM parsing. No pysam.
"""
import struct
import zlib
import json
import time
import os
from collections import defaultdict, Counter
from pathlib import Path

from ..bam.stream import read_bgzf_block, decode_bam_seq


def _parse_aux_tags(record, aux_start):
    """Parse BAM auxiliary tags. Returns dict of {tag: value}."""
    tags = {}
    pos = aux_start
    while pos < len(record) - 2:
        tag = record[pos:pos+2].decode('ascii', errors='replace')
        val_type = chr(record[pos+2])
        pos += 3
        if val_type == 'Z':
            end = record.index(b'\x00', pos)
            tags[tag] = record[pos:end].decode('ascii', errors='replace')
            pos = end + 1
        elif val_type == 'A':
            tags[tag] = chr(record[pos])
            pos += 1
        elif val_type in ('c', 'C'):
            tags[tag] = struct.unpack('<b' if val_type == 'c' else '<B', record[pos:pos+1])[0]
            pos += 1
        elif val_type in ('s', 'S'):
            tags[tag] = struct.unpack('<h' if val_type == 's' else '<H', record[pos:pos+2])[0]
            pos += 2
        elif val_type in ('i', 'I'):
            tags[tag] = struct.unpack('<i' if val_type == 'i' else '<I', record[pos:pos+4])[0]
            pos += 4
        elif val_type == 'f':
            tags[tag] = struct.unpack('<f', record[pos:pos+4])[0]
            pos += 4
        elif val_type == 'H':
            end = record.index(b'\x00', pos)
            tags[tag] = record[pos:end].hex()
            pos = end + 1
        elif val_type == 'B':
            sub = chr(record[pos])
            pos += 1
            count = struct.unpack('<I', record[pos:pos+4])[0]
            pos += 4
            sz = {'c': 1, 'C': 1, 's': 2, 'S': 2, 'i': 4, 'I': 4, 'f': 4}.get(sub, 1)
            pos += count * sz
        else:
            break
    return tags


def stream_unmapped(bam_path):
    """Stream ONLY unmapped reads from a BAM. Yields (read_name, seq, cell_barcode)."""
    with open(bam_path, 'rb') as fh:
        first_block = read_bgzf_block(fh)
        if first_block is None or first_block[:4] != b'BAM\x01':
            raise ValueError(f'Not a valid BAM: {bam_path}')

        data = first_block
        pos = 4
        l_text = struct.unpack('<i', data[pos:pos+4])[0]
        pos += 4 + l_text
        n_ref = struct.unpack('<i', data[pos:pos+4])[0]
        pos += 4
        for _ in range(n_ref):
            l_name = struct.unpack('<i', data[pos:pos+4])[0]
            pos += 4 + l_name + 4

        buffer = data[pos:]
        n_total = 0
        n_unmapped = 0

        while True:
            while len(buffer) < 4:
                block = read_bgzf_block(fh)
                if block is None:
                    return
                buffer += block

            block_size = struct.unpack('<i', buffer[:4])[0]
            while len(buffer) < block_size + 4:
                block = read_bgzf_block(fh)
                if block is None:
                    return
                buffer += block

            record = buffer[4:block_size+4]
            buffer = buffer[block_size+4:]
            n_total += 1

            if len(record) < 32:
                continue

            flag = struct.unpack('<H', record[14:16])[0]

            if not (flag & 4):
                continue

            n_unmapped += 1
            l_read_name = record[8]
            n_cigar_op = struct.unpack('<H', record[12:14])[0]
            l_seq = struct.unpack('<i', record[16:20])[0]

            read_name = record[32:32+l_read_name-1].decode('ascii', errors='replace')

            cigar_end = 32 + l_read_name + n_cigar_op * 4
            seq_start = cigar_end
            seq_bytes = (l_seq + 1) // 2
            seq_end = seq_start + seq_bytes
            if seq_end > len(record):
                continue

            seq = decode_bam_seq(record[seq_start:seq_end], l_seq)

            qual_end = seq_end + l_seq
            aux_start = qual_end

            cell_barcode = None
            try:
                tags = _parse_aux_tags(record, aux_start)
                cell_barcode = tags.get('CB') or tags.get('CR')
            except Exception:
                pass

            if cell_barcode is None:
                parts = read_name.split('_')
                if len(parts) >= 2 and len(parts[0]) >= 12:
                    cell_barcode = parts[0]

            yield read_name, seq, cell_barcode

            if n_unmapped % 50000 == 0:
                pct = 100 * n_unmapped / max(1, n_total)
                print(f'  {n_total:,} total, {n_unmapped:,} unmapped ({pct:.1f}%)...')


def scan_bam(bam_path, panel, min_hits=2, max_reads=None):
    """Scan a BAM for fungal transcripts. Returns result dict."""
    t0 = time.time()
    print(f'\n  STAFF Fungal Scanner')
    print(f'  Input: {bam_path}')
    print(f'  Panel: {len(panel.species_kmers)} species, k={panel.k}')

    per_cell = defaultdict(lambda: defaultdict(int))
    bulk_hits = defaultdict(int)
    n_unmapped = 0
    n_classified = 0
    n_total_reads = 0

    for read_name, seq, cell_barcode in stream_unmapped(bam_path):
        n_unmapped += 1
        if max_reads and n_unmapped > max_reads:
            break

        species = panel.best_species(seq, min_hits=min_hits)
        if species:
            n_classified += 1
            bulk_hits[species] += 1
            cb = cell_barcode or 'BULK'
            per_cell[cb][species] += 1

    elapsed = time.time() - t0

    print(f'\n  Results:')
    print(f'    Unmapped reads scanned: {n_unmapped:,}')
    print(f'    Classified as fungal:   {n_classified:,} ({100*n_classified/max(1,n_unmapped):.2f}%)')
    print(f'    Time: {elapsed:.1f}s')

    if bulk_hits:
        print(f'\n  Species breakdown:')
        for sp, count in sorted(bulk_hits.items(), key=lambda x: -x[1]):
            info = panel.species_info.get(sp, {})
            pct = 100 * count / max(1, n_classified)
            print(f'    {info.get("name","?"):30s}  {count:>8,} reads ({pct:.1f}%)  [{info.get("role","")}]')

    n_cells_with_fungus = sum(1 for cb in per_cell if cb != 'BULK')
    if n_cells_with_fungus:
        print(f'\n  Cells with fungal reads: {n_cells_with_fungus:,}')

    results = {
        'tool': 'STAFF_Fungal_Scanner',
        'version': '1.0',
        'source': str(bam_path),
        'n_unmapped_scanned': n_unmapped,
        'n_classified_fungal': n_classified,
        'fungal_pct_of_unmapped': round(100 * n_classified / max(1, n_unmapped), 4),
        'species_hits': {sp: int(c) for sp, c in bulk_hits.items()},
        'n_cells_with_fungus': n_cells_with_fungus,
        'per_cell': {cb: dict(sp_counts) for cb, sp_counts in per_cell.items()},
        'compute_time': round(elapsed, 1),
    }

    return results


def run(bam_path, output_path=None, max_reads=None, min_hits=2):
    """Full scan pipeline: load panel, scan BAM, save results."""
    from .panel import load_panel

    panel = load_panel()
    results = scan_bam(bam_path, panel, min_hits=min_hits, max_reads=max_reads)

    if output_path is None:
        output_path = str(bam_path).replace('.bam', '_fungal_scan.json')

    with open(output_path, 'w') as f:
        json.dump(results, f, indent=2, default=str)
    print(f'\n  Saved: {output_path}')

    return results
