# MONDO disease labels

MONDO release used for disease normalization and naming: **v2026-06-02**

- Release: https://github.com/monarch-initiative/mondo/releases/tag/v2026-06-02
- Asset: `mondo_nodes.tsv` (id -> name, plus synonyms and xrefs)

## Fetch (not committed)

The full `mondo_nodes.tsv` (~23 MB) is gitignored. Download it once:

```
python src/fetch.py mondo
```

Or pass a different release tag:

```
python src/fetch.py mondo --release v2026-06-02 --out data/mondo/mondo_nodes.tsv
```

## Use in the pipeline

```
python src/pipeline.py --hgnc hgnc_complete_set.txt --mondo-nodes data/mondo/mondo_nodes.tsv
```

When `--mondo-nodes` is supplied, two steps use it:

1. **Step 3 normalizes** EFO, DOID, Orphanet, and OTAR disease ids to MONDO using xrefs in
   the release. Ambiguous xrefs, where one external id points at more than one MONDO term,
   are excluded and those ids stay put. HP, OBA, GO, MP, and unmapped ids are left unchanged.
   Report: `output/mondo_normalization.csv`.
2. **Step 4 names** the MONDO disease nodes from the release. Deprecated terms are skipped
   here by default, so anything the release does not name falls through to the OLS4 pass.
   Report: `output/naming.csv`, rows with `label_source=mondo_release`.

`lookups.py` parses the file. `fetch.py` downloads it. `pipeline.py` uses it.
