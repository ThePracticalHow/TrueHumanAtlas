"""
Fungal reference panel -- download, cache, and build k-mer index.

Four pathogenic fungi. ~74 Mb total uncompressed. Fits in RAM as a Python set.
Each 31-mer is stored as a 62-bit integer (2 bits per base) for speed and memory.

Usage:
    from staff.scan.panel import build_panel, load_panel
    panel = build_panel()            # downloads if needed, builds index, caches
    panel = load_panel()             # loads cached index
    species = panel.classify(seq)    # returns species name or None
"""
import os
import gzip
import json
import struct
import pickle
import hashlib
import urllib.request
from pathlib import Path

K = 31

PANEL_DIR = Path.home() / ".staff" / "fungal_panel"

SPECIES = {
    "candida_albicans": {
        "name": "Candida albicans",
        "strain": "SC5314",
        "role": "tumor colonizer",
        "assembly": "GCF_000182965.3",
        "url": "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/182/965/GCF_000182965.3_ASM18296v3/GCF_000182965.3_ASM18296v3_genomic.fna.gz",
    },
    "aspergillus_fumigatus": {
        "name": "Aspergillus fumigatus",
        "strain": "Af293",
        "role": "lung decomposer",
        "assembly": "GCF_000002655.1",
        "url": "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/002/655/GCF_000002655.1_ASM265v1/GCF_000002655.1_ASM265v1_genomic.fna.gz",
    },
    "cryptococcus_neoformans": {
        "name": "Cryptococcus neoformans",
        "strain": "JEC21",
        "role": "brain infector",
        "assembly": "GCF_000091045.1",
        "url": "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/091/045/GCF_000091045.1_ASM9104v1/GCF_000091045.1_ASM9104v1_genomic.fna.gz",
    },
    "saccharomyces_cerevisiae": {
        "name": "Saccharomyces cerevisiae",
        "strain": "S288C",
        "role": "reference yeast",
        "assembly": "GCF_000146045.2",
        "url": "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/146/045/GCF_000146045.2_R64/GCF_000146045.2_R64_genomic.fna.gz",
    },
}

BASE_ENCODE = {'A': 0, 'C': 1, 'G': 2, 'T': 3,
               'a': 0, 'c': 1, 'g': 2, 't': 3}
COMP = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C',
        'a': 't', 't': 'a', 'c': 'g', 'g': 'c', 'N': 'N', 'n': 'n'}


def _encode_kmer(seq):
    """Encode a k-mer as a 64-bit integer. Returns None if ambiguous bases."""
    val = 0
    for ch in seq:
        b = BASE_ENCODE.get(ch)
        if b is None:
            return None
        val = (val << 2) | b
    return val


def _revcomp(seq):
    return ''.join(COMP.get(c, 'N') for c in reversed(seq))


def _canonical_kmer(seq):
    """Return the lexicographically smaller of kmer and its reverse complement."""
    rc = _revcomp(seq)
    return min(seq, rc)


def _parse_fasta_gz(path):
    """Yield sequences from a gzipped FASTA."""
    current = []
    with gzip.open(path, 'rt') as f:
        for line in f:
            if line.startswith('>'):
                if current:
                    yield ''.join(current)
                    current = []
            else:
                current.append(line.strip())
    if current:
        yield ''.join(current)


def _extract_kmers(fasta_path, k=K):
    """Extract all canonical k-mers from a gzipped FASTA. Returns a set of ints."""
    kmers = set()
    for seq in _parse_fasta_gz(fasta_path):
        for i in range(len(seq) - k + 1):
            kmer = seq[i:i+k]
            canon = _canonical_kmer(kmer)
            encoded = _encode_kmer(canon)
            if encoded is not None:
                kmers.add(encoded)
    return kmers


def download_genome(species_key, force=False):
    """Download a fungal genome FASTA. Returns path to local file."""
    info = SPECIES[species_key]
    PANEL_DIR.mkdir(parents=True, exist_ok=True)
    local_path = PANEL_DIR / f"{species_key}.fna.gz"
    if local_path.exists() and not force:
        return local_path
    print(f"  Downloading {info['name']} ({info['strain']})...")
    urllib.request.urlretrieve(info['url'], str(local_path))
    print(f"    {local_path.stat().st_size:,} bytes")
    return local_path


def download_all(force=False):
    """Download all fungal reference genomes."""
    paths = {}
    for key in SPECIES:
        paths[key] = download_genome(key, force=force)
    return paths


class FungalPanel:
    """In-memory k-mer index for fungal species classification.
    
    Optimized: single merged set for fast "is fungal?" screening,
    then per-species attribution only for positive k-mers.
    """

    def __init__(self):
        self.species_kmers = {}
        self.merged = set()         # union of all species k-mers for fast screening
        self.kmer_to_species = {}   # encoded_kmer -> species_key (for unique k-mers)
        self.species_info = SPECIES
        self.k = K

    def add_species(self, species_key, kmers):
        self.species_kmers[species_key] = kmers

    def build_merged_index(self):
        """Build the merged set and per-kmer species lookup. Call after all species added."""
        self.merged = set()
        self.kmer_to_species = {}
        seen_multi = set()
        for sp, kmers in self.species_kmers.items():
            for km in kmers:
                self.merged.add(km)
                if km in seen_multi:
                    continue
                if km in self.kmer_to_species:
                    del self.kmer_to_species[km]
                    seen_multi.add(km)
                else:
                    self.kmer_to_species[km] = sp
        print(f"  Merged index: {len(self.merged):,} total, {len(self.kmer_to_species):,} species-unique")

    def classify_read(self, seq):
        """Fast classify: screen against merged set, attribute species for hits."""
        if len(seq) < self.k:
            return {}
        hits = {}
        for i in range(len(seq) - self.k + 1):
            kmer = seq[i:i+self.k]
            canon = _canonical_kmer(kmer)
            encoded = _encode_kmer(canon)
            if encoded is None:
                continue
            if encoded in self.merged:
                sp = self.kmer_to_species.get(encoded)
                if sp:
                    hits[sp] = hits.get(sp, 0) + 1
                else:
                    for sp_key, kset in self.species_kmers.items():
                        if encoded in kset:
                            hits[sp_key] = hits.get(sp_key, 0) + 1
        return hits

    def best_species(self, seq, min_hits=2):
        """Return the best-matching species for a read, or None."""
        hits = self.classify_read(seq)
        if not hits:
            return None
        best = max(hits, key=hits.get)
        if hits[best] >= min_hits:
            return best
        return None

    def save(self, path=None):
        path = path or (PANEL_DIR / "panel_index.pkl")
        PANEL_DIR.mkdir(parents=True, exist_ok=True)
        with open(path, 'wb') as f:
            pickle.dump({
                'k': self.k,
                'species_kmers': {k: v for k, v in self.species_kmers.items()},
                'merged': self.merged,
                'kmer_to_species': self.kmer_to_species,
            }, f, protocol=pickle.HIGHEST_PROTOCOL)
        print(f"  Panel index saved: {path} ({path.stat().st_size:,} bytes)")

    @classmethod
    def load(cls, path=None):
        path = path or (PANEL_DIR / "panel_index.pkl")
        if not path.exists():
            return None
        with open(path, 'rb') as f:
            data = pickle.load(f)
        panel = cls()
        panel.k = data['k']
        panel.species_kmers = data['species_kmers']
        panel.merged = data.get('merged', set())
        panel.kmer_to_species = data.get('kmer_to_species', {})
        if not panel.merged:
            panel.build_merged_index()
        return panel

    def summary(self):
        total = len(self.merged)
        unique = len(self.kmer_to_species)
        print(f"\n  Fungal Panel: {len(self.species_kmers)} species, {total:,} unique {self.k}-mers ({unique:,} species-unique)")
        for sp, kmers in sorted(self.species_kmers.items()):
            info = SPECIES.get(sp, {})
            print(f"    {info.get('name','?'):30s}  {len(kmers):>12,} k-mers  ({info.get('role','')})")


def build_panel(force=False):
    """Download genomes, build k-mer index, cache it. Returns FungalPanel."""
    cached = PANEL_DIR / "panel_index.pkl"
    if cached.exists() and not force:
        print("  Loading cached panel index...")
        panel = FungalPanel.load(cached)
        if panel:
            panel.summary()
            return panel

    print("  Building fungal reference panel...")
    paths = download_all(force=force)

    panel = FungalPanel()
    for species_key, fasta_path in paths.items():
        info = SPECIES[species_key]
        print(f"  Indexing {info['name']}...")
        kmers = _extract_kmers(fasta_path)
        panel.add_species(species_key, kmers)
        print(f"    {len(kmers):,} unique {K}-mers")

    panel.build_merged_index()
    panel.save(cached)
    panel.summary()
    return panel


def load_panel():
    """Load existing panel or build if not cached."""
    return build_panel(force=False)
