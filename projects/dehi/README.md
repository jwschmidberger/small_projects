# DehI phylogenetic analysis

This project collects CMD, AhpD, and DehI-like proteins to explore the
hypothesis that DehI-like proteins arose through duplication and fusion of a
smaller ancestral unit related to the other groups.

The analysis script is [`dehi_phylogenetic_analysis.py`](dehi_phylogenetic_analysis.py).
All downloaded and derived files are written to [`results/`](results/).

## Setup

From the repository root, sync the Python environment with `uv` and install the
official MMseqs2 command-line program:

```bash
uv sync
brew install mmseqs2
```

On systems without Homebrew, install MMseqs2 using the platform's package
manager and ensure the `mmseqs` executable is on `PATH`.

## Download and cluster

The default run downloads fresh sequences and clusters each family separately
at 90% identity:

```bash
uv run python projects/dehi/dehi_phylogenetic_analysis.py
```

To cluster collections already present in `projects/dehi/results/` without
downloading them again:

```bash
uv run python projects/dehi/dehi_phylogenetic_analysis.py --cluster-only
```

The identity threshold is a fraction between `0.5` and `1.0`. For example:

```bash
# 70% identity
uv run python projects/dehi/dehi_phylogenetic_analysis.py --cluster-only --cluster-identity 0.70

# 60% identity
uv run python projects/dehi/dehi_phylogenetic_analysis.py --cluster-only --cluster-identity 0.60
```

The threshold is included in each output name. A 70% run produces files such
as `results/uniprot_PF02627_nr70.fasta` and
`results/uniprot_PF02627_nr70_clusters.tsv`. The FASTA contains one
representative per cluster; the TSV maps every input sequence to its
representative.

Clustering requires at least 80% bidirectional alignment coverage. This helps
prevent partial sequences from collapsing full-length proteins and preserves
multidomain architecture, which is important when testing the proposed DehI
duplication-and-fusion relationship.
