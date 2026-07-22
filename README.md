# small_projects
A sandbox for exploring and working on small projects

## DehI phylogenetic sequence collections

[`projects/dehi_phylogenetic_analysis1.py`](projects/dehi_phylogenetic_analysis1.py)
downloads the CMD, AhpD, and DehI protein collections and clusters each family
independently with MMseqs2. Keeping the collections separate preserves their
family labels for the later phylogenetic analysis.

### Setup

Sync the Python environment with `uv` and install the official MMseqs2 command-
line program:

```bash
uv sync
brew install mmseqs2
```

On systems without Homebrew, install MMseqs2 using the platform's package
manager and ensure the `mmseqs` executable is on `PATH`.

### Download and cluster

The default run downloads fresh sequences and clusters them at 90% identity:

```bash
uv run python projects/dehi_phylogenetic_analysis1.py
```

To cluster the existing FASTA collections without downloading them again:

```bash
uv run python projects/dehi_phylogenetic_analysis1.py --cluster-only
```

The identity threshold is a fraction between `0.5` and `1.0`. For example,
cluster at 70% or 60% identity with:

```bash
uv run python projects/dehi_phylogenetic_analysis1.py --cluster-only --cluster-identity 0.70
uv run python projects/dehi_phylogenetic_analysis1.py --cluster-only --cluster-identity 0.60
```

The threshold is included in each output name. A 70% run produces files such
as `uniprot_PF02627_nr70.fasta` and
`uniprot_PF02627_nr70_clusters.tsv`. The FASTA contains one representative per
cluster; the TSV maps every input sequence to its representative.

Clustering also requires at least 80% bidirectional alignment coverage. This
helps prevent partial sequences from collapsing full-length proteins and
preserves multidomain architecture, which is important when testing the
proposed DehI duplication-and-fusion relationship.
