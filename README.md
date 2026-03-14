# 🧬 SARS-CoV-2 Mutation Tracker

A full genomic-surveillance research project that tracks how SARS-CoV-2 spike protein
mutations emerge, spread geographically, and affect viral fitness over time.

---

## Project structure

```
covid_mutation_tracker/
│
├── config.py                      # Central config: paths, domains, variants, colours
├── requirements.txt               # Python dependencies
├── app.py                         # Streamlit dashboard (8 tabs)
│
├── pipeline/
│   ├── fetch_data.py              # Download reference genome + Nextstrain / GISAID data
│   ├── alignment.py               # Mutation calling (Biopython / minimap2 / metadata)
│   └── run_pipeline.py            # End-to-end orchestrator
│
├── analysis/
│   └── frequency_analysis.py      # Frequency, convergence, co-occurrence, hotspots, sweeps
│
├── visualization/
│   └── plots.py                   # All Plotly / Folium figures
│
├── tests/
│   └── test_pipeline.py           # Full unit + integration test suite
│
└── data/
    ├── raw/                       # Raw FASTA / metadata downloads
    └── processed/                 # Parquet files (pipeline output)
```

---

## Quick start (no internet required)

```bash
# 1. Install dependencies
pip install -r requirements.txt

# 2. Generate synthetic sample data and run the full pipeline
python pipeline/run_pipeline.py --mode sample --max-sequences 5000

# 3. Launch the dashboard
streamlit run app.py
```

Open [http://localhost:8501](http://localhost:8501) in your browser.

---

## Real data — Nextstrain (recommended)

No account required. Downloads pre-processed metadata with mutation annotations.

```bash
# Download and process 50,000 sequences
python pipeline/run_pipeline.py \
    --mode nextstrain \
    --max-sequences 50000

streamlit run app.py
```

> **Note:** The full metadata file is ~500 MB compressed. The pipeline processes
> and saves a filtered Parquet (~50–100 MB) for fast dashboard loading.

---

## Real data — GISAID

Requires a free account at [gisaid.org](https://gisaid.org).

```bash
# After downloading sequences.fasta from GISAID:
python pipeline/run_pipeline.py \
    --mode gisaid \
    --gisaid-fasta /path/to/sequences.fasta \
    --max-sequences 10000 \
    --alignment-backend biopython   # or minimap2 if installed

streamlit run app.py
```

---

## Step-by-step pipeline

You can run each step individually:

```bash
# Step 1: Fetch reference genome (MN908947) + metadata
python pipeline/fetch_data.py --source nextstrain --max-sequences 10000

# Step 2: Build mutation table
python pipeline/alignment.py --backend metadata

# Step 3: Compute all analyses
python analysis/frequency_analysis.py

# Step 4: Launch dashboard
streamlit run app.py
```

---

## Dashboard tabs

| Tab | Content |
|-----|---------|
| 📊 Overview | Metrics, variant dominance chart, hotspot table |
| 📈 Mutation Timeline | Per-mutation frequency curves + sweep speed |
| 🌍 Geographic Spread | Choropleth + first-detection table |
| 🔬 Spike Protein | Linear protein map + RBD zoom + domain breakdown |
| 🧬 Convergent Evolution | Mutations in multiple lineages + co-occurrence heatmap |
| 🌳 Phylogenetic Tree | Plotly variant tree + variant definition table |
| 🗂 Data Explorer | Searchable, filterable, downloadable mutation table |
| ℹ About | Pipeline docs + data sources + citations |

---

## Key analyses

### 1. Variant dominance curves
Monthly frequency of each WHO-labelled variant (stacked area chart).
Reproduces the well-known D614G → Alpha → Delta → Omicron succession.

### 2. Mutation sweep speed
Days for each key mutation to rise from 5% → 50% global frequency.
D614G swept the world in ~90 days; Omicron mutations swept in ~60 days.

### 3. Convergent evolution
Mutations appearing independently in ≥2 phylogenetically distinct lineages.
Classic example: E484K appearing in Beta, Gamma, and Iota independently.

### 4. Co-occurrence / linkage
Chi-squared test + odds ratio for mutation pair co-occurrence.
Mutations with high odds ratios define haplotype blocks (e.g. the Omicron constellation).

### 5. Spike mutational hotspots
Positional density map across all 1,273 amino acids.
RBD (306–527) and furin cleavage site (675–695) are the highest-density regions.

### 6. Geographic spread
First-detection date per variant per country with choropleth visualisation.

---

## Alignment backends

| Backend | Speed | Requires | Best for |
|---------|-------|----------|---------|
| `metadata` | Fastest | Nothing | Nextstrain/GISAID metadata with aaSubstitutions column |
| `biopython` | Slow | Nothing | Raw FASTA without pre-called mutations |
| `minimap2` | Fast | `minimap2` binary | Large raw FASTA batches |

Install minimap2:
```bash
# macOS
brew install minimap2

# Linux
sudo apt-get install minimap2

# Conda
conda install -c bioconda minimap2
```

---

## Running tests

```bash
python -m pytest tests/ -v
```

Expected output: 20+ tests covering config, mutation parsing, frequency analysis,
visualisation, and a lightweight end-to-end smoke test.

---

## Data sources

| Source | URL | What you get |
|--------|-----|-------------|
| Nextstrain open | https://nextstrain.org/ncov | Pre-processed metadata + mutation calls |
| GISAID | https://gisaid.org | Full sequences (requires free account) |
| NCBI | https://ncbi.nlm.nih.gov | Reference genome MN908947 (Wuhan-Hu-1) |
| CoV-Spectrum | https://cov-spectrum.org | Aggregated variant statistics |

---

## Tech stack

| Layer | Libraries |
|-------|-----------|
| Data | Biopython, pandas, numpy, pyarrow |
| Alignment | Biopython PairwiseAligner, minimap2 (optional) |
| Statistics | scipy (chi-squared, correlation) |
| Visualisation | Plotly, Folium, matplotlib, seaborn |
| Dashboard | Streamlit |
| Testing | unittest, pytest |

---

## Citation

> We gratefully acknowledge all data contributors to GISAID and Nextstrain.
> This project uses open genomic surveillance data for educational and research purposes.

GISAID citation:
> Elbe S & Buckland-Merrett G (2017). Data, disease and diplomacy: GISAID's innovative
> contribution to global health. *Global Challenges*, 1(1), 33–46.

Nextstrain citation:
> Hadfield J et al. (2018). Nextstrain: real-time tracking of pathogen evolution.
> *Bioinformatics*, 34(23), 4121–4123.
