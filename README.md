# Genetic discontinuities among microbial eukaryotes suggest that the 95% ANI may be a universal reference threshold for microbial species demarcation

[![DOI](https://zenodo.org/badge/1232240139.svg)](https://zenodo.org/badge/latestdoi/1232240139)

Ethan E. Hetrick, Charlotte J. Royer, Matthew H. Seabolt, and Konstantinos T. Konstantinidis

Publication: TBD
DOI: see badge above

This repository contains the processed data, metadata, code, and figure exports needed to reproduce the manuscript analyses.

## Contents

- `ANI-microbial-eukaryotes.ipynb`: main notebook for analyses and figure generation.
- `data/`: processed pairwise comparison database and compact supporting data tables.
- `assets/genome_taxonomy.csv`: genome-to-taxonomy metadata.
- `assets/genera_classification_sources_environmental.csv`: genus annotations used for Figure 1.
- `bin/`: helper scripts called by the notebook:
  - `figure2_panel_a_rank_probabilities.R`
  - `figure2_compose_alignment_panels.py`
  - `figure3_ani_ridgeline_by_genus.R`
  - `publication_figure_style.R`
  - `publication_figure_style.py`
- `requirements.txt` and `requirements-r.txt`: Python and R dependencies.

## Reproducing the Notebook

Run the notebook from the repository root:

```bash
python -m pip install -r requirements.txt
jupyter lab ANI-microbial-eukaryotes.ipynb
```

The R helper scripts require the packages listed in `requirements-r.txt`.

## Data

The main pairwise comparison table is stored as a split Parquet dataset in `data/all_tables_processed/`. Read it directly as one dataset:

```python
import pyarrow.dataset as ds

dataset = ds.dataset("data/all_tables_processed", format="parquet")
df = dataset.to_table(columns=["Genome1", "Genome2", "ANI", "AF"]).to_pandas()
```
