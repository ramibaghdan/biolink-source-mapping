# Mapping Decisions

This document records the non-obvious choices made while mapping two public sources
(DGIdb drug-gene interactions and Open Targets gene-disease associations) onto the
Biolink Model. The straightforward mappings are in `category_map.csv` and
`predicate_map.csv`. The cases below needed a judgment call, which is described here.

## 1. Reconciling gene identifiers across two sources

DGIdb identifies genes by symbol (for example BRCA1). Open Targets identifies genes by
Ensembl ID (for example ENSG00000012048). To integrate the two sources into one graph,
the same gene has to resolve to one node, so the two identifier schemes had to be
reconciled.

I used the HGNC complete set as the bridge, mapping each DGIdb gene symbol to its
Ensembl ID, then keying genes by Ensembl ID for the join and normalizing to NCBIGene
for the final Biolink node identifier.

Two cases did not resolve cleanly:
- Some symbols mapped to more than one Ensembl ID. I kept the HGNC-approved primary
  mapping and dropped the ambiguous alternates rather than duplicating the gene.
- Some symbols had no Ensembl ID in HGNC (withdrawn or non-standard symbols). Those
  rows could not be joined to Open Targets and were excluded from the connected graph.

This reconciliation is the reason the sampling is gene anchored: rather than sampling
each source independently (which would share almost no genes), genes present in both
sources were selected first, so every gene in the graph has both a drug edge and a
disease edge.

## 2. Normalizing mixed disease identifier prefixes

Open Targets disease IDs do not all use one ontology. The sample contains EFO, MONDO,
DOID, Orphanet, and OTAR prefixes. Biolink prefers MONDO for disease nodes.

Decisions made:
- Converted the source underscore format to the CURIE colon format (EFO_0003950 becomes
  EFO:0003950), since Biolink identifiers are CURIEs.
- Originally kept each disease in its original ontology when no cross-reference was available.
  **Update:** EFO, DOID, Orphanet, and OTAR ids now normalize to MONDO when `mondo_nodes.tsv`
  provides a unique xref (see section below). Unmapped ids and phenotype ontologies (HP, OBA,
  GO, MP) remain in their source ontology.

### MONDO disease name enrichment (added)

Open Targets MONDO disease nodes were initially stored with the raw MONDO id as the display
name (for example `MONDO:0002108`). The MONDO release asset `mondo_nodes.tsv`
(v2026-06-02) provides authoritative id → label mappings. When `--mondo-nodes` is passed
to `map_to_biolink.py`, MONDO-prefixed disease nodes are enriched with the release label
(for example `breast carcinoma`). Obsolete MONDO terms are skipped for naming; unmatched
ids keep the CURIE as the name.

### OLS4 labels for non-MONDO disease nodes (added)

After MONDO naming, 273 of 1,096 disease nodes still displayed a bare CURIE (mostly EFO,
plus HP, OBA, Orphanet, GO, MP, and one MONDO id absent from the pinned release). Those
ids are resolved once against the EBI Ontology Lookup Service (OLS4) and cached in
`data/ontology_labels/labels.tsv`. `ontology_labels.py` applies the cache offline and
leaves already-named nodes untouched. All 273 resolved; `check_labels.py` confirms every
disease node has a readable name.

### Cross-ontology normalization to MONDO (added)

Open Targets uses mixed disease ontologies. Using xrefs in MONDO release `mondo_nodes.tsv`
(v2026-06-02), EFO, DOID, Orphanet, and OTAR ids are rewritten to MONDO where there is a
unique xref match. Example: `EFO:0009491` → `MONDO:0000004` (adrenocortical insufficiency).

Decisions made:
- Only unambiguous xrefs are used (one external id → exactly one non-obsolete MONDO term).
  Ambiguous xrefs are excluded from the map rather than guessing.
- Phenotype and non-disease ontologies in the sample (HP, OBA, GO, MP) are not normalized;
  they often describe traits or measurements, not disorders, and lack reliable MONDO xrefs.
- Ids with no xref match stay in their original ontology.
- After normalization, MONDO labels are applied to all MONDO disease nodes.

## 3. Mapping DGIdb interaction types to Biolink predicates

DGIdb interaction types come from many underlying sources, so the vocabulary is
inconsistent and includes near-synonyms. The full mapping is in `predicate_map.csv`;
the choices worth explaining:

- inhibitor and blocker both map to `biolink:negatively_regulates`. They are
  functionally equivalent here (both reduce target activity), so they were collapsed to
  one predicate. The source distinction is recorded in the mapping table rather than the
  graph.
- agonist, activator, and potentiator all map to `biolink:positively_regulates` for the
  same reason in the other direction.
- modulator maps to the broader parent `biolink:regulates`. A modulator has no inherent
  direction, and Biolink guidance reserves the vague parent for exactly the case where
  direction is unknown or could be both. This is the one place a non-directional
  predicate is the correct choice.
- negative modulator and positive modulator, by contrast, do carry direction, so they
  map to the directional `biolink:negatively_regulates` and `biolink:positively_regulates`
  rather than the vague parent. The contrast with plain modulator is intentional.
- binder maps to `biolink:physically_interacts_with`, not a regulatory predicate,
  because binding alone states physical contact with no functional direction.
- cleavage maps to the broad `biolink:affects`, since it names a mechanism rather than a
  clean up or down direction.

## 4. Handling "other/unknown" differently by edge type

DGIdb includes an explicit "other/unknown" interaction type. Rather than drop those rows
(dropping data is itself a silent decision), they are mapped to the most general
drug-gene predicate, `biolink:interacts_with`.

For gene-disease edges, the equivalent "no specific relationship known" case is mapped
to the gene-disease association predicate instead, because a gene and a disease are
related by association, not by interaction. Keeping these two general cases distinct by
edge type avoids implying a drug-gene-style interaction between a gene and a disease.

## 5. A source data-quality issue: "antisense oligonucleotide"

One DGIdb interaction type, "antisense oligonucleotide", is not an interaction type at
all. It describes the modality of the drug, not how the drug relates to the gene. Rather than invent a predicate
for it, the row is left unmapped and excluded from the edge output, and the issue is
recorded here. Surfacing this kind of source inconsistency is part of the value of a
mapping pass.

## 6. Non-disease terms typed as diseases

Open Targets associates genes with more than diseases. After MONDO normalization, 100 of
the 1,096 disease nodes come from ontologies that do not describe diseases: 69 HP and 1
MP phenotypes, 29 OBA biological attributes, and one GO biological process
(`GO:0009410`, response to xenobiotic stimulus). All carry `biolink:Disease` from the
ingest.

Decision: report the mismatch rather than silently rewriting it. `ontology_labels.py`
emits a `suggested_category` column (`biolink:PhenotypicFeature` for HP/MP/OBA,
`biolink:BiologicalProcessOrActivity` for GO). Splitting the category is a modelling
change that affects edge semantics — a gene-phenotype association is not a gene-disease
association — and belongs in its own pass, not in a labelling step.

## On modeling style

These mappings use simple directional predicates (negatively_regulates,
positively_regulates, physically_interacts_with, and so on). Biolink also supports a
richer qualifier-based pattern for chemical-gene effects, where the predicate is
`biolink:affects` plus separate qualifiers for direction, aspect, and causal mechanism
(for example agonism). The qualifier pattern carries more detail but is heavier to read
and produce. The simple predicates were chosen here for clarity and because they capture
the direction that matters for this graph. The qualifier-based pattern is a planned
enhancement, alongside the optional LLM-assisted mapping layer.
