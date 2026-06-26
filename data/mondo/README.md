# MONDO disease labels

MONDO release used for disease name lookup: **v2026-06-02**

- Release: https://github.com/monarch-initiative/mondo/releases/tag/v2026-06-02
- Asset: `mondo_nodes.tsv` (id → name, plus synonyms and xrefs)

## Fetch (not committed)

The full `mondo_nodes.tsv` (~23 MB) is gitignored. Download it once:

```
python src/fetch_mondo_release.py
```

Or pass a different release tag:

```
python src/fetch_mondo_release.py --release v2026-06-02 --out data/mondo/mondo_nodes.tsv
```

## Use in the pipeline

```
python src/map_to_biolink.py --hgnc hgnc_complete_set.txt --mondo-nodes data/mondo/mondo_nodes.tsv
```

When `--mondo-nodes` is supplied, the pipeline:

1. **Normalizes** EFO, DOID, Orphanet, and OTAR disease ids to MONDO using xrefs in the
   release (report: `output/mondo_normalization.csv`). HP, OBA, GO, MP, and unmapped ids
   are left unchanged.
2. **Enriches** all MONDO disease nodes with human-readable labels (report:
   `output/mondo_name_enrichment.csv`).
