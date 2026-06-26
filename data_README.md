# Data

This project uses small samples of two public biomedical sources, plus a gene identifier
map. The full source files are not committed; only the small samples in data/raw/ are.

## Sources

### DGIdb (Drug-Gene Interaction Database)
- Drug-gene interactions with interaction types (inhibitor, agonist, blocker, etc.).
- Open source. Downloaded from https://www.dgidb.org/ (Downloads section).
- Used for the drug-gene edges.

### Open Targets
- Gene-disease associations (Ensembl gene id, disease id, association score).
- Licensed CC0. Downloaded from https://platform.opentargets.org/downloads
- Used for the gene-disease edges.

### HGNC complete set (gene identifier map)
- Maps gene symbols and Ensembl ids to Entrez ids.
- Downloaded from https://www.genenames.org/download/
- Used to unify genes from both sources to a single NCBIGene identifier.

### MONDO (disease name lookup)
- Maps MONDO ids to disease labels from the official release asset `mondo_nodes.tsv`.
- Release **v2026-06-02**: https://github.com/monarch-initiative/mondo/releases/tag/v2026-06-02
- Used to replace raw MONDO ids with human-readable disease names in `output/nodes.csv`.
- See [data/mondo/README.md](mondo/README.md) for fetch instructions.

## What is committed

- data/raw/dgidb_sample.csv: a gene-anchored sample of DGIdb drug-gene rows.
- data/raw/opentargets_sample.csv: a gene-anchored sample of Open Targets gene-disease rows.

The two samples are anchored on the same set of genes (genes present in both sources), so
the integrated graph connects drug to gene to disease. They were produced by
sample_sources.py from the full downloads.

## What is not committed

The full DGIdb, Open Targets, HGNC, and MONDO `mondo_nodes.tsv` files are large and are
excluded by .gitignore.
To regenerate the samples, download those three files and run sample_sources.py (see the
top-level README). Verify each source's current license and download location at the time
you download, since these can change.

## Note

Downloaded on first build for this project. These are public research datasets used here
only to demonstrate a data-to-Biolink mapping method, not redistributed in full.
