# Non-MONDO ontology labels

Human-readable names for the disease nodes the MONDO pass leaves unnamed.

## Why this exists

Step 3 of the pipeline maps EFO, DOID, Orphanet, and OTAR ids to MONDO where a unique
xref exists, and step 4 names every node that became MONDO. Ids with no unique xref stay
in their source ontology and keep the id as their name:

```
EFO:0006890,biolink:Disease,EFO:0006890
HP:0040187,biolink:Disease,HP:0040187
```

A node named `HP:0040187` is unusable. Unreadable to a person, and any downstream consumer
returns database codes where an answer belongs. Roughly a quarter of disease nodes are in
this state after the MONDO pass.

## Source

Labels come from the **EBI Ontology Lookup Service (OLS4)**: <https://www.ebi.ac.uk/ols4>

One resolver covers every ontology still present after MONDO normalization (EFO, HP, OBA,
GO, MP, Orphanet) rather than a downloader per ontology.

Endpoint: `https://www.ebi.ac.uk/ols4/api/terms?obo_id={CURIE}`

When several ontologies return the same `obo_id`, because importing ontologies mirror each
other's terms, the term is ranked non-obsolete first, then defining ontology. OLS names
obsolete terms "obsolete <name>", so a live label from an importing ontology is the better
display name even though the defining ontology is more authoritative. The chosen ontology is
recorded per row.

## Fetch (committed)

```
python src/fetch.py labels
```

Unlike `data/mondo/mondo_nodes.tsv` (~23 MB, gitignored), `labels.tsv` covers only the ids
this graph contains and is small enough to commit. That keeps step 4 offline and the
pipeline reproducible without network access.

It reads `output/nodes.csv` to decide which ids need a label, so the pipeline has to have
run at least once first.

The fetch is resumable. Ids already cached are skipped, including ids previously recorded
as `not_found`. Re-query everything with `--refresh`, or smoke-test with `--limit 10`.

## Apply

Naming happens in step 4 of the pipeline, using this cache:

```
python src/pipeline.py --hgnc hgnc_complete_set.txt --mondo-nodes data/mondo/mondo_nodes.tsv
python src/validate.py     # exits 1 if any node is still named by its id
```

Report: `output/naming.csv`, rows with `label_source=ols4`.

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

Several ids reaching this stage are not diseases. Open Targets associates genes with HP and
MP phenotypes, OBA attributes, and one GO biological process, all of which the ingest types
as `biolink:Disease`.

Step 4 reports this in the `suggested_category` column and changes nothing. Correcting the
category is a mapping decision, not a labelling one. See `mappings/mapping_decisions.md`.
