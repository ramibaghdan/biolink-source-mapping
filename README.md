# Biolink Source Ingestion and Mapping

Ingests two public biomedical sources and maps their entities and relationships onto the
[Biolink Model](https://biolink.github.io/biolink-model/), the open schema used to build
interoperable biomedical knowledge graphs. The interesting part here are
the mapping decisions: where a source does not map cleanly onto Biolink, those judgment
calls are documented in [mappings/mapping_decisions.md](mappings/mapping_decisions.md).

   ![Drug-gene-disease neighborhood for MAP2K2 (disease names from MONDO enrichment)](figures/subgraph_6416.png)

## Demo

A drug-gene interaction and a gene-disease association, from two independent sources, are
unified onto one gene node so they connect:

```
DRUG_NAME:OLAPARIB --biolink:negatively_regulates--> NCBIGene:672 --biolink:gene_associated_with_condition--> EFO:0000305
        (DGIdb)                                          (shared gene)                                      (Open Targets)
```

The raw rows that produced this, and how each field was mapped, are in
[output/before_after_examples.md](output/before_after_examples.md).


## Why this matters

Biomedical data is spread across many sources that each use their own vocabularies and
identifier schemes. Building a knowledge graph means reconciling those sources onto one
shared model so they can be queried together. Biolink is the community standard for that
model. This project shows the core of that reconciliation on a small scale: take a drug-gene
source and a gene-disease source, normalize their identifiers, map their relationship
vocabularies to Biolink predicates, and produce a connected set of Biolink triples, while
documenting the cases that required a decision.

## What it does

1. Ingests a sample of DGIdb drug-gene interactions and Open Targets gene-disease associations.
2. Unifies genes from both sources to a single NCBIGene identifier using the HGNC map, so the
   same gene is one node and the two sources can connect.
3. Maps DGIdb interaction types (inhibitor, agonist, modulator, and so on) to Biolink predicates,
   using directional predicates where the direction is known and the broad parent in cases it isn't.
4. Normalizes disease identifiers to CURIE form, mapping EFO, DOID, Orphanet, and OTAR ids to
   MONDO where a unique cross-reference exists and keeping the rest in their original ontology.
5. Names every disease node. MONDO nodes take their label from the MONDO release; ids that stay
   in EFO, HP, OBA, GO, MP, or Orphanet are resolved against OLS4, so no node displays a bare
   identifier where a name belongs.
6. Writes Biolink nodes and edges, and validates that every category and predicate is a real
   Biolink term.

## Repository layout

```
biolink_source_mapping/
├── sample_sources.py            # builds the small samples from the full downloads
├── mappings/
│   ├── category_map.csv         # source entity type -> Biolink category
│   ├── predicate_map.csv        # source relationship -> Biolink predicate
│   └── mapping_decisions.md     # the judgment calls (the core of the project)
├── data/
│   ├── raw/                     # the committed small samples
│   ├── mondo/                   # MONDO release assets for disease label lookup
│   ├── ontology_labels/         # committed OLS4 label cache for non-MONDO ids
│   └── README.md                # source provenance and licenses
├── src/
│   ├── ingest.py                # load and clean the samples
│   ├── map_to_biolink.py        # apply the mappings, produce nodes and edges
│   ├── mondo_labels.py          # MONDO id -> disease name from mondo_nodes.tsv
│   ├── mondo_normalize.py       # EFO/DOID/Orphanet -> MONDO via release xrefs
│   ├── fetch_mondo_release.py   # download mondo_nodes.tsv from a MONDO release
│   ├── fetch_ontology_labels.py # OLS4 label lookup for non-MONDO ids (network, once)
│   ├── ontology_labels.py       # apply those labels to nodes.csv
│   ├── check_labels.py          # fails if any node is still named by its id
│   ├── validate.py              # confirm Biolink terms, print a summary
│   ├── visualize_subgraph.py    # drug-gene-disease neighborhood figure
│   └── requirements.txt
├── figures/
│   └── subgraph_6416.png        # demo neighborhood figure
└── output/
    ├── nodes.csv                # Biolink entities
    ├── edges.csv                # Biolink triples
    ├── mondo_normalization.csv  # EFO/DOID/Orphanet -> MONDO xref report
    ├── mondo_name_enrichment.csv # MONDO id -> label lookup report
    ├── ontology_name_enrichment.csv # non-MONDO id -> label, with category flags
    └── before_after_examples.md # raw row -> triple, with explanation
```

## Run it

From the project root:

```
pip install -r src/requirements.txt

# 1. build the samples from the full downloads (see data/README.md for where to get them)
python sample_sources.py --dgidb dgidb_interactions.tsv --opentargets opentargets.parquet --genemap hgnc_complete_set.txt --n-genes 150

# 2. produce the Biolink nodes and edges
python src/fetch_mondo_release.py   # once: download MONDO v2026-06-02 labels
python src/map_to_biolink.py --hgnc hgnc_complete_set.txt --mondo-nodes data/mondo/mondo_nodes.tsv

# 3. name the disease nodes MONDO did not cover
python src/fetch_ontology_labels.py  # once: OLS4 lookup, writes a committed cache
python src/ontology_labels.py        # applies the labels (offline)

# 4. confirm no node is still named by its id
python src/check_labels.py

# 5. validate and summarize
python src/validate.py

# 6. optional: regenerate the demo figure (uses disease names from nodes.csv)
python src/visualize_subgraph.py --gene NCBIGene:6416
```

Install `bmt` (the Biolink Model Toolkit, included in requirements) so validation checks
against the authoritative Biolink model rather than the built-in fallback list.

## Mapping decisions

The mapping tables handle the routine cases. The decisions that needed reasoning can be found in [mappings/mapping_decisions.md](mappings/mapping_decisions.md), including reconciling
gene identifiers across two schemes, normalizing mixed disease ontology prefixes, mapping an
inconsistent interaction-type vocab to Biolink predicates, handling unknown relationships
differently by edge type, non-disease terms arriving typed as diseases, and a source
data-quality issue.

## Data and license

Uses small samples of DGIdb (open) and Open Targets (CC0), plus the HGNC gene map. Disease
labels for non-MONDO ids come from the EBI Ontology Lookup Service (OLS4). The full downloads
are not committed. See [data/README.md](data/README.md).

## Planned enhancements

- ~~Normalize disease identifiers to a single ontology (MONDO) using a cross-reference file.~~
  **Done.** EFO, DOID, Orphanet, and OTAR ids map to MONDO via release xrefs (`v2026-06-02`)
  and receive MONDO labels. Ids with no unique xref stay in their source ontology and are
  named from OLS4, so every disease node is readable regardless of which ontology it ended up
  in. Any ids that resolve to no label anywhere are documented in
  [mappings/mapping_decisions.md](mappings/mapping_decisions.md).
- Split the non-disease terms out of `biolink:Disease`. Open Targets associates genes with HP
  and MP phenotypes, OBA attributes, and at least one GO biological process, all of which the
  ingest currently types as diseases. `ontology_labels.py` reports the mismatch in a
  `suggested_category` column without changing it; correcting it affects edge semantics and
  belongs in its own pass.
- An optional LLM-assisted step that proposes Biolink mappings for source terms not yet in the
  mapping tables, with a human accepting or rejecting each suggestion. The current version is
  fully deterministic.
