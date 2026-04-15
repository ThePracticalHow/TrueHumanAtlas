"""
Adaptive view engine -- project the atlas onto any axes chosen by the question.

The database IS the dimension space. Any column is an axis.
The view engine discovers available columns and never needs updating.

Usage:
    staff view ribo_indep vs k_rg color=source
    staff view k_gl vs k_rl color=disease
    staff view fungal_load vs ribo_indep color=tissue
    staff view ribo_indep vs k_rg vs det_k color=source    # 3D
    staff view ribo_indep                                   # histogram
    staff view ribo_indep vs k_rg where="tissue='brain'"
"""
import sqlite3
import re
import sys
import os
import numpy as np

from .atlas.db import open_atlas, default_db_path

# Column name aliases for convenience
ALIASES = {
    'r_ind': 'ribo_indep', 'rind': 'ribo_indep', 'independence': 'ribo_indep',
    'det': 'det_k', 'determinant': 'det_k',
    'cells': 'n_cells', 'ncells': 'n_cells',
    'krm': 'k_rm', 'krn': 'k_rn', 'krg': 'k_rg',
    'kmn': 'k_mn', 'kmg': 'k_mg', 'kng': 'k_ng',
    'krl': 'k_rl', 'kml': 'k_ml', 'knl': 'k_nl', 'kgl': 'k_gl',
}


def _resolve_alias(name):
    return ALIASES.get(name.lower().strip(), name.lower().strip())


def discover_columns(conn):
    """Discover all plottable columns from the joined cell_states + coupling_tensor view."""
    row = conn.execute("""
        SELECT * FROM cell_states cs
        JOIN coupling_tensor ct ON ct.cell_state_id = cs.id
        LIMIT 1
    """).fetchone()
    if not row:
        return []
    desc = conn.execute("""
        SELECT * FROM cell_states cs
        JOIN coupling_tensor ct ON ct.cell_state_id = cs.id
        LIMIT 0
    """).description
    return [d[0] for d in desc]


def _categorize_columns(columns):
    """Split into numeric-likely and categorical-likely columns."""
    numeric_hints = {'ribo_indep', 'det_k', 'n_cells', 'n_reads',
                     'k_rm', 'k_rn', 'k_rg', 'k_mn', 'k_mg', 'k_ng',
                     'k_rl', 'k_ml', 'k_nl', 'k_gl',
                     'eigenvalue_min', 'eigenvalue_max', 'compute_time',
                     'fungal_pct', 'n_classified_fungal', 'n_unmapped_scanned',
                     'candida_hits', 'aspergillus_hits', 'crypto_hits', 'yeast_hits'}
    categorical = {'tissue', 'disease', 'sex', 'category', 'species', 'source',
                   'platform', 'dataset', 'name', 'id', 'cell_state_id', 'condition',
                   'age_range', 'failure_mode'}
    num = [c for c in columns if c in numeric_hints or c.startswith('k_')]
    cat = [c for c in columns if c in categorical]
    return num, cat


def fetch_data(conn, x_col, y_col=None, z_col=None, color_col=None, size_col=None, where=None):
    """Fetch data from the atlas for the requested columns."""
    cols_needed = [x_col]
    if y_col: cols_needed.append(y_col)
    if z_col: cols_needed.append(z_col)
    if color_col: cols_needed.append(color_col)
    if size_col: cols_needed.append(size_col)
    cols_needed.append('cs.name')
    cols_needed.append('cs.source')

    select = ', '.join(set(cols_needed))
    sql = f"""
        SELECT {select}
        FROM cell_states cs
        JOIN coupling_tensor ct ON ct.cell_state_id = cs.id
        WHERE ct.ribo_indep IS NOT NULL
    """
    if where:
        sql += f" AND ({where})"
    sql += f" ORDER BY {x_col} ASC"

    return conn.execute(sql).fetchall(), [d.split('.')[-1] for d in set(cols_needed)]


def parse_view_args(args):
    """Parse free-form view arguments into structured query.

    Formats:
        "ribo_indep vs k_rg color=source"
        "ribo_indep vs k_rg vs det_k color=source"
        "ribo_indep"
        "k_rm vs k_rg where=tissue='brain' color=disease size=n_cells"
    """
    text = ' '.join(args) if isinstance(args, list) else args
    text = text.strip()

    x = y = z = color = size = where = None

    # Extract key=value pairs
    for pattern, setter in [
        (r'color\s*=\s*(\S+)', 'color'),
        (r'size\s*=\s*(\S+)', 'size'),
        (r'where\s*=\s*["\']?(.+?)["\']?\s*(?:color|size|$)', 'where'),
    ]:
        m = re.search(pattern, text, re.IGNORECASE)
        if m:
            if setter == 'color': color = _resolve_alias(m.group(1))
            elif setter == 'size': size = _resolve_alias(m.group(1))
            elif setter == 'where': where = m.group(1).strip()
            text = text[:m.start()] + text[m.end():]

    # Clean up remaining text for axis names
    text = text.strip()
    parts = re.split(r'\s+vs\.?\s+|\s+versus\s+', text)
    parts = [_resolve_alias(p.strip()) for p in parts if p.strip()]

    if len(parts) >= 1: x = parts[0]
    if len(parts) >= 2: y = parts[1]
    if len(parts) >= 3: z = parts[2]

    if not color:
        color = 'source'

    return x, y, z, color, size, where


def view(x, y=None, z=None, color='source', size=None, where=None, db_path=None,
         output=None, show=True):
    """
    Project the atlas onto question-defined axes.
    Any column from cell_states or coupling_tensor can be an axis.
    """
    conn = open_atlas(db_path)
    all_cols = discover_columns(conn)
    num_cols, cat_cols = _categorize_columns(all_cols)

    # Qualify column names for the SQL join
    def q(col):
        if col in ('name', 'id', 'tissue', 'disease', 'sex', 'category',
                    'species', 'dataset', 'source', 'n_cells', 'n_reads',
                    'platform', 'age_range'):
            return f'cs.{col}'
        return f'ct.{col}'

    # Build SELECT
    select_cols = [q(x)]
    col_names = [x]
    if y:
        select_cols.append(q(y)); col_names.append(y)
    if z:
        select_cols.append(q(z)); col_names.append(z)
    if color:
        select_cols.append(q(color)); col_names.append(color)
    if size:
        select_cols.append(q(size)); col_names.append(size)
    select_cols.append('cs.name as _name')
    select_cols.append('cs.source as _source')
    col_names.extend(['_name', '_source'])

    sql = f"""SELECT {', '.join(select_cols)}
              FROM cell_states cs
              JOIN coupling_tensor ct ON ct.cell_state_id = cs.id
              WHERE ct.ribo_indep IS NOT NULL"""
    if where:
        sql += f" AND ({where})"

    rows = conn.execute(sql).fetchall()
    conn.close()

    if not rows:
        print(f"  No data found for: x={x}, y={y}, where={where}")
        print(f"  Available numeric columns: {', '.join(num_cols)}")
        print(f"  Available categorical columns: {', '.join(cat_cols)}")
        return

    print(f"\n  View: {x}" + (f" vs {y}" if y else "") + (f" vs {z}" if z else ""))
    print(f"  Color: {color}  |  Filter: {where or 'none'}  |  {len(rows)} points")

    # Extract columns by position
    data = {}
    for i, name in enumerate(col_names):
        vals = [row[i] for row in rows]
        data[name] = vals

    # Detect if we need matplotlib
    try:
        import matplotlib
        matplotlib.use('Agg') if output else None
        import matplotlib.pyplot as plt

        fig = plt.figure(figsize=(12, 8))

        # Color mapping
        color_vals = data.get(color, data.get('_source'))
        unique_colors = sorted(set(str(v) for v in color_vals if v is not None))
        cmap = plt.cm.get_cmap('tab10', max(len(unique_colors), 1))
        color_map = {v: cmap(i) for i, v in enumerate(unique_colors)}

        # Front-aware markers
        source_vals = data.get('_source', [])
        def _marker(src):
            if src == 'mycelium_atlas': return 'X'
            if src == 'hongyan_atlas': return 'v'
            return 'o'

        x_vals = np.array([float(v) if v is not None else np.nan for v in data[x]])

        if y and not z:
            # 2D scatter
            ax = fig.add_subplot(111)
            y_vals = np.array([float(v) if v is not None else np.nan for v in data[y]])
            for cv in unique_colors:
                mask = np.array([str(c) == cv for c in color_vals])
                for src_type in set(source_vals):
                    smask = mask & np.array([s == src_type for s in source_vals])
                    if smask.any():
                        ax.scatter(x_vals[smask], y_vals[smask],
                                   c=[color_map[cv]], marker=_marker(src_type),
                                   s=60, alpha=0.7, edgecolors='white', linewidths=0.5,
                                   label=f'{cv} [{src_type[:3]}]' if smask.sum() > 0 else None)
            ax.set_xlabel(x, fontsize=12, fontweight='bold')
            ax.set_ylabel(y, fontsize=12, fontweight='bold')
            ax.set_title(f'Atlas View: {x} vs {y}', fontsize=14)
            ax.grid(True, alpha=0.3)

        elif y and z:
            # 3D scatter
            ax = fig.add_subplot(111, projection='3d')
            y_vals = np.array([float(v) if v is not None else np.nan for v in data[y]])
            z_vals = np.array([float(v) if v is not None else np.nan for v in data[z]])
            for cv in unique_colors:
                mask = np.array([str(c) == cv for c in color_vals])
                if mask.any():
                    ax.scatter(x_vals[mask], y_vals[mask], z_vals[mask],
                               c=[color_map[cv]], s=40, alpha=0.7, label=cv)
            ax.set_xlabel(x)
            ax.set_ylabel(y)
            ax.set_zlabel(z)
            ax.set_title(f'Atlas View: {x} vs {y} vs {z}')

        else:
            # 1D histogram
            ax = fig.add_subplot(111)
            valid = x_vals[~np.isnan(x_vals)]
            ax.hist(valid, bins=30, color='#00ff00', alpha=0.7, edgecolor='black')
            ax.set_xlabel(x, fontsize=12, fontweight='bold')
            ax.set_ylabel('Count', fontsize=12)
            ax.set_title(f'Atlas Distribution: {x}', fontsize=14)
            ax.grid(True, alpha=0.3)

        # Legend (keep manageable)
        handles, labels = ax.get_legend_handles_labels() if hasattr(ax, 'get_legend_handles_labels') else ([], [])
        if 0 < len(labels) <= 20:
            ax.legend(fontsize=8, loc='best', framealpha=0.8)

        plt.tight_layout()
        if output:
            plt.savefig(output, dpi=150, bbox_inches='tight', facecolor='black')
            print(f"  Saved: {output}")
        if show and not output:
            plt.show()
        plt.close()

    except ImportError:
        # No matplotlib -- text output
        print(f"\n  (matplotlib not available -- text output)")
        print(f"  {'x':>10}  {'y':>10}  {'color':>15}  name")
        print(f"  {'-'*60}")
        for i in range(min(30, len(rows))):
            xv = data[x][i]
            yv = data[y][i] if y else ''
            cv = data.get(color, [''])[i] if color else ''
            nm = data['_name'][i][:35] if '_name' in data else ''
            print(f"  {str(xv):>10}  {str(yv):>10}  {str(cv):>15}  {nm}")


def list_columns(db_path=None):
    """Print all available columns that can be used as axes."""
    conn = open_atlas(db_path)
    all_cols = discover_columns(conn)
    num_cols, cat_cols = _categorize_columns(all_cols)
    conn.close()

    print(f"\n  Available axes (numeric): {', '.join(sorted(num_cols))}")
    print(f"  Available colors (categorical): {', '.join(sorted(cat_cols))}")
    print(f"\n  Usage: staff view <x> vs <y> color=<cat> where=\"<filter>\"")
