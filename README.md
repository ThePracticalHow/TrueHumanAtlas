# The True Human Atlas

> *Every cell gets a position. Zero free parameters. One measurement. Open source.*

## What This Is

A cross-kingdom coupling tensor atlas mapping the operator coupling of every human cell state — and comparing it to fungi and viruses on the same axis. The coupling tensor measures how tightly a cell's four major subsystems (RIBO, MITO, NUC, GOLGI) co-vary. A single number — RIBO independence — places any cell on a continuous axis from maximally coupled (0.039, viral infection) to maximally decoupled (0.890, glioblastoma).

**The same measurement works across kingdoms.** Yeast at active growth: 0.105. Cryptococcus: 0.224-0.567. Human cancer: 0.310-0.890. Same axis. Same physics. Zero parameters.

## Quick Start

```bash
pip install -e ".[expr]"

# Compute coupling tensor on any h5ad
staff tensor your_data.h5ad

# Query the atlas
staff ask "compare cancer vs fungus vs virus"

# Detect fungal signatures in expression data
staff detect your_data.h5ad

# Ingest results into the atlas
staff ingest results.json
```

## The Independence Axis

```
0.039   SARS-CoV-2 infected      (virus hyper-couples host)
0.105   Yeast active growth       (maximum fungal coupling)
0.154   Human proliferative       (healthy baseline)
0.220   Normal lung tissue        
0.310   Primary lung cancer       (defense activating)
0.508   Metastatic cancer         (escaping)
0.714   Yeast lag phase           (stressed fungus = cancer range)
0.743   Chronic HCV               (virus exhausts host)
0.890   Glioblastoma              (maximum decoupling)
```

## What's In the Atlas

| Kingdom | Species | Conditions | Cells/Samples |
|---------|---------|-----------|---------------|
| Human | *H. sapiens* | Cancer, normal, fetal, senescent, viral infection | 8.5M+ |
| Fungal | *S. cerevisiae*, *C. neoformans*, *C. albicans* | Active growth, lag, stress | 8,000+ |
| Viral | Influenza A, HCV, SARS-CoV-2 | Infected host cells | 250+ |

## Tools

| Command | What It Does |
|---------|-------------|
| `staff tensor` | Compute 4×4 or 5×5 coupling tensor from h5ad |
| `staff detect` | Expression-based fungal presence scoring |
| `staff scan` | Fungal transcript scanner from BAM unmapped reads |
| `staff ingest` | Auto-detect and ingest any data format into atlas |
| `staff ask` | Natural language atlas queries |
| `staff view` | Morphing visualization (any column = any axis) |
| `staff summary` | Atlas statistics |
| `staff fronts` | Three-front map (human vs fungal vs viral) |
| `staff axis` | Independence axis display |
| `staff base4` | Base-4 XOR at splice junctions |
| `staff frame` | Reading frame analysis |
| `staff splice` | Per-molecule splice decision chains |
| `staff binary` | Full binary transcriptome |

## The Coupling Tensor

Four operators define the cell's economy:

| Operator | Genes | What It Measures |
|----------|-------|-----------------|
| **RIBO** | RPS*, RPL* | Translation machinery |
| **MITO** | MT-* | Energy budget |
| **NUC** | Everything else | Transcriptional regulation |
| **GOLGI** | GOLGA/B, SEC61, COPA/B, MAN1, MGAT, GORASP | Secretory pathway |

K[i,j] = Spearman correlation between operator i and j total expression across single cells.

**RIBO independence** = 1 - mean(|K_RM|, |K_RN|, |K_RG|)

Zero parameters. Zero thresholds. Zero training. Same input = same output. Always.

## Data Sources

All data derived from public datasets:

- **SGNex** (Singapore Nanopore Expression): H9, K562, HepG2 direct RNA
- **GSE131907**: Korean NSCLC cohort (208,506 cells)
- **GSE226225**: WI-38 fetal fibroblast senescence
- **GSE250041**: Proliferative vs senescent
- **GSE131928**: Glioblastoma (28 tumors)
- **GSE174083**: Meditation retreat (PNAS 2021)
- **GSE144820**: Yeast growth conditions
- **GSE67602**: Cryptococcus atlas
- **CellxGene Census**: Pan-tissue, pan-disease human atlas
- **GSE97672, GSE84346**: Viral infection (Influenza, HCV)

## Visualization

Open `globe/index.html` in a browser for the interactive 3D battlefield map (Three.js, WGS84 oblate spheroid).

Open `globe/morphglobe.html` for the PCA-on-sphere morphing view — drag the slider to morph between geographic space and coupling tensor space.

## License

Code: MIT. Data: CC-BY 4.0.

## Citation

```
Leng, J. (2026). The True Human Atlas: Cross-Kingdom Coupling Tensor Atlas.
```

See CITATION.cff for full citation metadata.

## Related Work

- Mazan-Mamczarz K, Wind EJ, **Leng J**, Gorospe M. (2026). Intercompartmental communication in senescence. *FEBS Open Bio*. doi:10.1002/2211-5463.70236
