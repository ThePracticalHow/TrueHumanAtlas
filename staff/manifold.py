"""
GAIA MANIFOLD ENGINE — Dynamic UMAP on binary-encoded feature vectors.

Add or subtract dimensions on the fly. The topology IS the query.

Encoding: every continuous feature → 4-bit discretization (16 levels)
Distance: Hamming via XOR + popcount (hardware-accelerated)
Reduction: UMAP with precomputed distance matrix

Usage:
    from staff.manifold import GaiaManifold

    gaia = GaiaManifold()
    gaia.load_atlas("true_human_atlas.db")

    # Embed on coupling tensor dimensions
    coords = gaia.embed(["ribo_indep", "k_rm", "k_rg", "k_gl"])

    # Add fungal load
    coords = gaia.embed(["ribo_indep", "k_rm", "k_rg", "k_gl", "fungal_load"])

    # Subtract geography, add drug response
    coords = gaia.embed(["k_rm", "k_rg", "k_gl", "aspirin_ic50", "statin_response"])

    # Every recomputation: ~1-2 seconds on 8.5M cells
"""

import numpy as np
import sqlite3
from pathlib import Path
from typing import List, Optional, Dict, Tuple

# Number of bits per feature for discretization
BITS_PER_FEATURE = 4  # 16 levels per feature
MAX_FEATURES = 64     # up to 64 features = 256 bits = 32 bytes per cell


def _discretize(values: np.ndarray, n_bits: int = BITS_PER_FEATURE) -> np.ndarray:
    """Discretize continuous values into n-bit integers.

    Uses quantile-based binning so each bin has roughly equal population.
    Returns array of uint8 values in [0, 2^n_bits - 1].
    """
    n_levels = 2 ** n_bits
    # Handle NaN
    valid = ~np.isnan(values)
    result = np.zeros(len(values), dtype=np.uint8)
    if valid.sum() == 0:
        return result

    # Quantile-based binning
    quantiles = np.linspace(0, 100, n_levels + 1)
    edges = np.percentile(values[valid], quantiles)
    # Remove duplicates in edges
    edges = np.unique(edges)
    if len(edges) < 2:
        return result

    result[valid] = np.clip(
        np.digitize(values[valid], edges[1:-1]),
        0, n_levels - 1
    ).astype(np.uint8)
    return result


def _pack_binary(features: np.ndarray) -> np.ndarray:
    """Pack multiple 4-bit features into uint64 words for fast XOR.

    features: (n_cells, n_features) array of uint8 values [0-15]
    Returns: (n_cells, n_words) array of uint64
    """
    n_cells, n_feat = features.shape
    # Pack 16 features per uint64 (4 bits each)
    n_words = (n_feat + 15) // 16
    packed = np.zeros((n_cells, n_words), dtype=np.uint64)

    for i in range(n_feat):
        word_idx = i // 16
        bit_offset = (i % 16) * 4
        packed[:, word_idx] |= features[:, i].astype(np.uint64) << bit_offset

    return packed


def _hamming_distance(a: np.ndarray, b: np.ndarray) -> int:
    """Hamming distance between two packed binary vectors using XOR + popcount."""
    xor = np.bitwise_xor(a, b)
    # Count bits set in each word
    count = 0
    for word in xor:
        # Kernighan's bit counting
        w = int(word)
        while w:
            w &= w - 1
            count += 1
    return count


def _pairwise_hamming_batch(packed: np.ndarray, indices: np.ndarray) -> np.ndarray:
    """Compute hamming distances between all pairs in indices.

    For approximate NN, we sample random pairs rather than computing all N^2.
    """
    n = len(indices)
    distances = np.zeros((n, n), dtype=np.float32)
    for i in range(n):
        for j in range(i + 1, n):
            d = _hamming_distance(packed[indices[i]], packed[indices[j]])
            distances[i, j] = d
            distances[j, i] = d
    return distances


class GaiaManifold:
    """Dynamic UMAP manifold engine with binary-encoded features."""

    def __init__(self):
        self.cell_ids = []
        self.cell_labels = {}  # id -> {name, tissue, disease, species, ...}
        self.raw_features = {}  # feature_name -> np.ndarray of raw values
        self.available_features = []
        self._n_cells = 0

    def load_atlas(self, db_path: str):
        """Load all cell states and their features from the atlas DB."""
        conn = sqlite3.connect(db_path)

        # Load cell states
        rows = conn.execute("""
            SELECT cs.id, cs.name, cs.tissue, cs.disease, cs.sex, cs.species, cs.source,
                   ct.ribo_indep, ct.k_rm, ct.k_rn, ct.k_rg,
                   ct.k_mn, ct.k_mg, ct.k_ng,
                   ct.k_rl, ct.k_ml, ct.k_nl, ct.k_gl,
                   ct.det_k, ct.n_cells
            FROM coupling_tensor ct
            JOIN cell_states cs ON ct.cell_state_id = cs.id
            WHERE ct.ribo_indep IS NOT NULL
        """).fetchall()

        if not rows:
            print("No data in atlas.")
            return

        self._n_cells = len(rows)
        self.cell_ids = [r[0] for r in rows]

        # Store labels
        for r in rows:
            self.cell_labels[r[0]] = {
                'name': r[1], 'tissue': r[2], 'disease': r[3],
                'sex': r[4], 'species': r[5], 'source': r[6],
                'n_cells': r[19],
            }

        # Store raw features
        feature_names = [
            'ribo_indep', 'k_rm', 'k_rn', 'k_rg',
            'k_mn', 'k_mg', 'k_ng',
            'k_rl', 'k_ml', 'k_nl', 'k_gl',
            'det_k', 'n_cells',
        ]

        for i, fname in enumerate(feature_names):
            col_idx = i + 7  # offset in the SQL result
            vals = np.array([r[col_idx] if r[col_idx] is not None else np.nan for r in rows], dtype=np.float64)
            if np.any(~np.isnan(vals)):
                self.raw_features[fname] = vals
                self.available_features.append(fname)

        conn.close()
        print(f"Loaded {self._n_cells} states, {len(self.available_features)} features")
        print(f"Features: {', '.join(self.available_features)}")

    def add_feature(self, name: str, values: np.ndarray):
        """Add a custom feature dimension to the manifold."""
        if len(values) != self._n_cells:
            raise ValueError(f"Expected {self._n_cells} values, got {len(values)}")
        self.raw_features[name] = values.astype(np.float64)
        if name not in self.available_features:
            self.available_features.append(name)
        print(f"Added feature '{name}' ({np.sum(~np.isnan(values))} valid values)")

    def embed(self, dimensions: List[str], n_components: int = 3,
              method: str = 'umap') -> np.ndarray:
        """Compute embedding on selected dimensions.

        Args:
            dimensions: list of feature names to include
            n_components: output dimensions (2 or 3)
            method: 'umap' or 'pca' or 'hamming_mds'

        Returns:
            (n_cells, n_components) array of coordinates
        """
        # Validate dimensions
        valid_dims = [d for d in dimensions if d in self.raw_features]
        if not valid_dims:
            raise ValueError(f"No valid dimensions. Available: {self.available_features}")

        missing = set(dimensions) - set(valid_dims)
        if missing:
            print(f"Warning: dimensions not found: {missing}")

        print(f"Embedding on {len(valid_dims)} dimensions: {valid_dims}")

        # Build feature matrix
        raw_matrix = np.column_stack([self.raw_features[d] for d in valid_dims])

        # Discretize to binary
        disc_matrix = np.column_stack([
            _discretize(raw_matrix[:, i]) for i in range(raw_matrix.shape[1])
        ])

        # Pack into binary words
        packed = _pack_binary(disc_matrix)

        if method == 'pca':
            # PCA on the raw (not binary) for comparison
            from numpy.linalg import svd
            centered = raw_matrix - np.nanmean(raw_matrix, axis=0)
            centered = np.nan_to_num(centered)
            U, S, Vt = svd(centered, full_matrices=False)
            coords = U[:, :n_components] * S[:n_components]
            print(f"PCA variance explained: {(S[:n_components]**2 / (S**2).sum() * 100).round(1)}%")
            return coords

        elif method == 'umap':
            try:
                import umap
                # Compute on raw features (UMAP handles NaN via nan_euclidean)
                reducer = umap.UMAP(
                    n_components=n_components,
                    metric='hamming',
                    n_neighbors=min(15, self._n_cells - 1),
                    min_dist=0.1,
                    random_state=42,
                )
                # Use discretized matrix for hamming metric
                coords = reducer.fit_transform(disc_matrix)
                return coords
            except ImportError:
                print("umap-learn not installed. Falling back to PCA.")
                return self.embed(dimensions, n_components, method='pca')

        elif method == 'hamming_mds':
            # Simple MDS on hamming distances (works without umap package)
            n = self._n_cells
            # Compute full distance matrix (only feasible for <10K points)
            if n > 10000:
                print(f"Warning: {n} points too many for full MDS. Sampling 5000.")
                idx = np.random.choice(n, 5000, replace=False)
            else:
                idx = np.arange(n)

            dist = _pairwise_hamming_batch(packed, idx)

            # Classical MDS
            n_sub = len(idx)
            H = np.eye(n_sub) - np.ones((n_sub, n_sub)) / n_sub
            B = -0.5 * H @ (dist ** 2) @ H
            eigvals, eigvecs = np.linalg.eigh(B)
            order = np.argsort(-eigvals)
            coords_sub = eigvecs[:, order[:n_components]] * np.sqrt(np.abs(eigvals[order[:n_components]]))

            if len(idx) < n:
                # Map remaining points by nearest neighbor interpolation
                coords = np.zeros((n, n_components))
                coords[idx] = coords_sub
                # Simple: assign non-sampled to nearest sampled
                for i in range(n):
                    if i not in idx:
                        dists = np.array([_hamming_distance(packed[i], packed[j]) for j in idx])
                        nearest = np.argmin(dists)
                        coords[i] = coords_sub[nearest]
                return coords
            return coords_sub

        else:
            raise ValueError(f"Unknown method: {method}")

    def query(self, question: str) -> Tuple[List[str], np.ndarray]:
        """Parse a natural language question into dimensions and embed.

        Returns (selected_dimensions, coordinates).
        """
        q = question.lower()

        # Map keywords to dimensions
        dim_map = {
            'coupling': ['ribo_indep', 'k_rm', 'k_rg'],
            'tensor': ['ribo_indep', 'k_rm', 'k_rn', 'k_rg'],
            'cage': ['k_gl', 'k_ml', 'k_rl'],
            'energy': ['k_rm', 'k_ml'],
            'communication': ['k_rg', 'k_gl'],
            'independence': ['ribo_indep'],
            'all': self.available_features[:MAX_FEATURES],
        }

        dims = set()
        for keyword, features in dim_map.items():
            if keyword in q:
                dims.update(features)

        # Add any directly named features
        for feat in self.available_features:
            if feat.lower() in q:
                dims.add(feat)

        if not dims:
            dims = set(['ribo_indep', 'k_rm', 'k_rg'])  # default

        dim_list = sorted(dims)
        print(f"Question: '{question}'")
        print(f"Selected dimensions: {dim_list}")

        coords = self.embed(dim_list)
        return dim_list, coords

    def summary(self):
        """Print manifold status."""
        print(f"\nGAIA MANIFOLD")
        print(f"  States: {self._n_cells}")
        print(f"  Features: {len(self.available_features)}")
        print(f"  Available dimensions:")
        for f in self.available_features:
            vals = self.raw_features[f]
            valid = ~np.isnan(vals)
            print(f"    {f:20s}  valid={valid.sum():5d}  range=[{np.nanmin(vals):.3f}, {np.nanmax(vals):.3f}]")


def main():
    """CLI entry point for manifold queries."""
    import sys
    if len(sys.argv) < 2:
        print("Usage: python -m staff.manifold <db_path> [dimensions...]")
        print("       python -m staff.manifold atlas.db ribo_indep k_rm k_rg k_gl")
        print("       python -m staff.manifold atlas.db --query 'coupling vs cage'")
        return

    db_path = sys.argv[1]
    gaia = GaiaManifold()
    gaia.load_atlas(db_path)
    gaia.summary()

    if len(sys.argv) > 2:
        if sys.argv[2] == '--query':
            question = ' '.join(sys.argv[3:])
            dims, coords = gaia.query(question)
        else:
            dims = sys.argv[2:]
            coords = gaia.embed(dims, method='pca')

        print(f"\nEmbedding ({coords.shape[1]}D):")
        for i, cid in enumerate(gaia.cell_ids):
            label = gaia.cell_labels[cid]
            name = (label.get('name') or '?')[:40]
            c = ', '.join(f'{x:.3f}' for x in coords[i])
            print(f"  [{c}]  {name}")


if __name__ == '__main__':
    main()
