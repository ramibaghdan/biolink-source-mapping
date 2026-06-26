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

MONDO-prefixed disease nodes in `output/nodes.csv` get human-readable names from the
release. Non-MONDO disease ids (EFO, DOID, Orphanet, OTAR) are unchanged. The enrichment
report is written to `output/mondo_name_enrichment.csv`.
