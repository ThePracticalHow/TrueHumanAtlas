"""
STAFF CLI — unified entry point.

Usage:
    staff ingest <path> [--tissue T] [--disease D] [--sex S] [--species SP] [--db DB]
    staff ask "<natural language query>"
    staff tensor <h5ad> [--group KEY] [--output PATH]
    staff binary <bam> [--output PATH] [--max-reads N]
    staff base4 <bam> [--output PATH]
    staff frame <bam> [--output PATH]
    staff splice <bam> [--output PATH]
    staff summary
    staff fronts
    staff axis
"""
import sys
import argparse


def cmd_ingest(args):
    from .ingest import ingest
    meta = {}
    if args.tissue:
        meta['tissue'] = args.tissue
    if args.disease:
        meta['disease'] = args.disease
    if args.sex:
        meta['sex'] = args.sex
    if args.species:
        meta['species'] = args.species
    if args.group:
        meta['group_key'] = args.group
    if args.category:
        meta['category'] = args.category

    ingest(args.path, db_path=args.db, meta=meta or None)


def cmd_ask(args):
    query_string = ' '.join(args.query)
    from .atlas.query import ask
    ask(query_string, db_path=args.db)


def cmd_tensor(args):
    from .expr.coupling_tensor import run
    out = args.output or args.h5ad.replace('.h5ad', '_tensor.json')
    run(args.h5ad, output_path=out, group_key=args.group)


def cmd_binary(args):
    from .bam.binary_transcripter import run_and_save
    out = args.output or args.bam.replace('.bam', '_binary.json')
    run_and_save(args.bam, output_path=out, max_reads=args.max_reads)


def cmd_base4(args):
    from .bam.base4 import run_and_save
    out = args.output or args.bam.replace('.bam', '_base4.json')
    run_and_save(args.bam, output_path=out)


def cmd_frame(args):
    from .bam.reading_frame import run_and_save
    out = args.output or args.bam.replace('.bam', '_frame.json')
    run_and_save(args.bam, output_path=out)


def cmd_splice(args):
    from .bam.splice_history import run_and_save
    out = args.output or args.bam.replace('.bam', '_splice.json')
    run_and_save(args.bam, output_path=out)


def cmd_scan(args):
    from .scan.scanner import run as run_scan
    results = run_scan(args.bam, output_path=args.output,
                       max_reads=args.max_reads, min_hits=args.min_hits)
    if args.with_tensor:
        from .scan.correlate import correlate_with_h5ad
        corr = correlate_with_h5ad(results, args.with_tensor, group_key=args.group)
        out_corr = (args.output or args.bam.replace('.bam', '_fungal_scan.json')).replace('.json', '_correlation.json')
        import json
        with open(out_corr, 'w') as f:
            json.dump(corr, f, indent=2, default=str)
        print(f'  Correlation saved: {out_corr}')


def cmd_detect(args):
    from .detect import detect
    out = args.output or args.h5ad.replace('.h5ad', '_goddess_detect.json')
    detect(args.h5ad, group_key=args.group, output_path=out)


def cmd_scan_panel(args):
    from .scan.panel import build_panel
    build_panel(force=args.force)


def cmd_view(args):
    from .view import view, parse_view_args, list_columns
    query_parts = args.query
    if not query_parts or query_parts == ['columns']:
        list_columns(db_path=args.db)
        return
    x, y, z, color, size, where = parse_view_args(query_parts)
    if not x:
        list_columns(db_path=args.db)
        return
    view(x, y=y, z=z, color=color, size=size, where=where,
         db_path=args.db, output=args.output, show=not args.output)


def cmd_summary(args):
    from .atlas.db import open_atlas, summary
    conn = open_atlas(args.db)
    summary(conn)
    conn.close()


def cmd_fronts(args):
    """Delegate to build_atlas.py --front for the full three-front map."""
    import subprocess
    from pathlib import Path
    tools = Path(__file__).resolve().parent.parent / "tools"
    subprocess.run([sys.executable, str(tools / 'build_atlas.py'), '--front'])


def cmd_axis(args):
    import subprocess
    from pathlib import Path
    tools = Path(__file__).resolve().parent.parent / "tools"
    subprocess.run([sys.executable, str(tools / 'build_atlas.py'), '--axis'])


def main():
    parser = argparse.ArgumentParser(
        prog='staff',
        description='STAFF -- Spliceosome Tensor Atlas of Functional Features',
    )
    parser.add_argument('--db', default=None, help='Path to atlas database')
    sub = parser.add_subparsers(dest='command')

    # ingest
    p_ingest = sub.add_parser('ingest', help='Ingest data into the atlas')
    p_ingest.add_argument('path', help='File or folder to ingest (.bam, .h5ad, .json, folder/)')
    p_ingest.add_argument('--tissue', help='Tissue type')
    p_ingest.add_argument('--disease', help='Disease label')
    p_ingest.add_argument('--sex', help='Sex (male/female/mixed)')
    p_ingest.add_argument('--species', help='Species')
    p_ingest.add_argument('--category', help='Category (cancer, normal, fetal, etc.)')
    p_ingest.add_argument('--group', help='Group column for h5ad')
    p_ingest.set_defaults(func=cmd_ingest)

    # ask
    p_ask = sub.add_parser('ask', help='Query the atlas in natural language')
    p_ask.add_argument('query', nargs='+', help='Your question')
    p_ask.set_defaults(func=cmd_ask)

    # tensor
    p_tensor = sub.add_parser('tensor', help='Compute coupling tensor from h5ad')
    p_tensor.add_argument('h5ad', help='Path to h5ad file')
    p_tensor.add_argument('--group', help='Group column name')
    p_tensor.add_argument('--output', help='Output JSON path')
    p_tensor.set_defaults(func=cmd_tensor)

    # BAM tools
    for name, helptext in [('binary', 'Binary transcripter'), ('base4', 'Base-4 XOR'),
                            ('frame', 'Reading frame'), ('splice', 'Splice history')]:
        p = sub.add_parser(name, help=helptext)
        p.add_argument('bam', help='Path to BAM file')
        p.add_argument('--output', help='Output JSON path')
        if name == 'binary':
            p.add_argument('--max-reads', type=int, help='Max reads to process')
        p.set_defaults(func=locals().get(f'cmd_{name}', lambda a: print(f'Not implemented: {name}')))

    p_binary = sub._name_parser_map.get('binary')
    if p_binary:
        p_binary.set_defaults(func=cmd_binary)
    p_base4 = sub._name_parser_map.get('base4')
    if p_base4:
        p_base4.set_defaults(func=cmd_base4)
    p_frame = sub._name_parser_map.get('frame')
    if p_frame:
        p_frame.set_defaults(func=cmd_frame)
    p_splice = sub._name_parser_map.get('splice')
    if p_splice:
        p_splice.set_defaults(func=cmd_splice)

    # view -- adaptive projection engine
    p_view = sub.add_parser('view', help='Project atlas onto any axes. The UMAP that morphs.')
    p_view.add_argument('query', nargs='*', help='e.g. "ribo_indep vs k_rg color=source" or "columns"')
    p_view.add_argument('--output', help='Save plot to file instead of displaying')
    p_view.set_defaults(func=cmd_view)

    # detect -- expression-based Goddess detection (works on any h5ad)
    p_detect = sub.add_parser('detect', help='Detect the Goddess in expression data (h5ad). No BAM needed.')
    p_detect.add_argument('h5ad', help='Path to h5ad file')
    p_detect.add_argument('--output', help='Output JSON path')
    p_detect.add_argument('--group', help='Group column for per-group breakdown')
    p_detect.set_defaults(func=cmd_detect)

    # scan -- fungal transcript scanner
    p_scan = sub.add_parser('scan', help='Scan BAM for fungal transcripts')
    p_scan.add_argument('bam', help='Path to BAM file')
    p_scan.add_argument('--output', help='Output JSON path')
    p_scan.add_argument('--max-reads', type=int, help='Max unmapped reads to scan')
    p_scan.add_argument('--min-hits', type=int, default=2, help='Min k-mer hits to classify (default 2)')
    p_scan.add_argument('--with-tensor', help='h5ad file to correlate fungal load with coupling')
    p_scan.add_argument('--group', help='Group column for per-group correlation')
    p_scan.set_defaults(func=cmd_scan)

    p_scan_panel = sub.add_parser('scan-panel', help='Download/build fungal reference panel')
    p_scan_panel.add_argument('--force', action='store_true', help='Force re-download')
    p_scan_panel.set_defaults(func=cmd_scan_panel)

    # summary / fronts / axis
    p_summary = sub.add_parser('summary', help='Atlas summary')
    p_summary.set_defaults(func=cmd_summary)
    p_fronts = sub.add_parser('fronts', help='Three-front map')
    p_fronts.set_defaults(func=cmd_fronts)
    p_axis = sub.add_parser('axis', help='Full independence axis')
    p_axis.set_defaults(func=cmd_axis)

    args = parser.parse_args()
    if not args.command:
        parser.print_help()
        return

    args.func(args)


if __name__ == '__main__':
    main()
