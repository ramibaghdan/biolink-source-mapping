# Non-MONDO ontology labels

Human-readable names for the disease nodes that the MONDO pass leaves unnamed.

## Why this exists

`mondo_normalize.py` maps EFO, DOID, Orphanet, and OTAR ids to MONDO where a
unique xref exists, and `mondo_labels.py` names every node that became MONDO.
Ids with no unique xref stay in their source ontology and keep the id as their
name:

```
EFO:0006890,biolink:Disease,EFO:0006890
HP:0040187,biolink:Disease,HP:0040187
```

A node named `HP:0040187` is unusable — unreadable to a person, and any
downstream consumer returns database codes where an answer belongs. Roughly a
quarter of disease nodes are in this state after the MONDO pass.

## Source

Labels come from the **EBI Ontology Lookup Service (OLS4)**:
<https://www.ebi.ac.uk/ols4>

One resolver covers every ontology still present after MONDO normalization
(EFO, HP, OBA, GO, MP, Orphanet), rather than a downloader per ontology.

Endpoint used: `https://www.ebi.ac.uk/ols4/api/terms?obo_id={CURIE}`

When several ontologies return the same `obo_id` (importing ontologies mirror
each other's terms), the term is chosen by ranking non-obsolete above obsolete
first, then the defining ontology above an importing one. Obsolete terms are
labelled "obsolete <name>" in OLS, so a live label from an importing ontology
is the more useful display name even though the defining ontology is more
authoritative. The chosen ontology is recorded per row.

## Fetch (committed)

```
python src/fetch_ontology_labels.py
```

Unlike `data/mondo/mondo_nodes.tsv` (~23 MB, gitignored), `labels.tsv` covers
only the ids this graph actually contains and is small enough to commit. That
keeps `ontology_labels.py` and `check_labels.py` fully offline and the pipeline
reproducible without network access.

The fetch is resumable — ids already cached are skipped, including ids
previously recorded as `not_found`. Re-query everything with `--refresh`, or
smoke-test with `--limit 10`.

## Apply

```
python src/ontology_labels.py            # writes output/nodes.csv
python src/ontology_labels.py --dry-run  # report only
python src/check_labels.py               # exits 1 if any node is still named by id
```

Report: `output/ontology_name_enrichment.csv`

## Schema

`labels.tsv`, tab-separated:

| column | meaning |
|---|---|
| `curie` | the id being resolved, e.g. `HP:0040187` |
| `label` | human-readable name; empty when unresolved |
| `ontology` | OLS ontology the label was taken from |
| `obsolete` | `true` when the chosen term is obsolete |
| `iri` | full IRI of the chosen term |
| `status` | `ok` or `not_found` |

## A category observation

Several ids reaching this stage are not diseases. Open Targets associates genes
with HP and MP phenotypes, OBA attributes, and at least one GO biological
process, all of which the ingest types as `biolink:Disease`.

`ontology_labels.py` reports this in the `suggested_category` column and changes
nothing. Correcting the categories is a mapping decision, not a labelling one —
see `mappings/mapping_decisions.md`.
