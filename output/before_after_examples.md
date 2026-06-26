# Before and After: Raw Source to Biolink Triple

Real examples from the generated output, showing how a raw source row becomes a
Biolink style triple, and how two independent sources connect on a shared gene.

The shared gene in these examples is NCBIGene:6416 (MAP2K2, also called MEK2).

## Drug-gene interaction (from DGIdb)

Raw DGIdb row:

```
gene_symbol,drug_name,interaction_type
MAP2K2,TRAMETINIB,inhibitor
```

Becomes (one edge in edges.csv):

```
subject,predicate,object,subject_category,object_category,knowledge_source,source_relationship
DRUG_NAME:TRAMETINIB,biolink:negatively_regulates,NCBIGene:6416,biolink:ChemicalEntity,biolink:Gene,DGIdb,inhibitor
```

What happened in the mapping:
- The gene symbol MAP2K2 was resolved to NCBIGene:6416 via the HGNC map.
- The source interaction type "inhibitor" was mapped to biolink:negatively_regulates.
- The drug is the subject and the gene is the object, matching the drug-acts-on-gene direction.

Several MEK inhibitors in the sample map to this same gene the same way, including
Trametinib, Selumetinib, and Cobimetinib.

## Gene-disease association (from Open Targets)

Raw Open Targets row:

```
ensembl_id,disease_efo_id,association_score
ENSG00000126934,EFO_0020865,...
```

Becomes:

```
subject,predicate,object,subject_category,object_category,knowledge_source,source_relationship
NCBIGene:6416,biolink:gene_associated_with_condition,EFO:0020865,biolink:Gene,biolink:Disease,OpenTargets,associated_with
```

What happened:
- The Ensembl id ENSG00000126934 was resolved to NCBIGene:6416 via HGNC, the same id the
  gene received from the DGIdb side, so the two sources connect on one gene node.
- The disease id was converted from underscore to colon form (EFO_0020865 to EFO:0020865).
  Some disease ids in the sample are then normalized to MONDO via release xrefs (see below);
  ids without a unique xref (including EFO:0020865) stay in their original ontology.

## Cross-ontology normalization to MONDO

Open Targets also uses EFO and Orphanet ids. Using xrefs from MONDO release v2026-06-02,
matching disease ids are rewritten to MONDO before labels are applied.

Before:

```
id,category,name
EFO:0009491,biolink:Disease,EFO:0009491
```

After normalization + label enrichment:

```
id,category,name
MONDO:0000004,biolink:Disease,adrenocortical insufficiency
```

What happened:
- `mondo_nodes.tsv` lists `EFO:0009491` as an xref on `MONDO:0000004`.
- The pipeline maps the edge object and node id to MONDO, then applies the MONDO label.
- In this sample, 405 of 577 EFO/DOID/Orphanet ids normalized; 172 had no unique xref.
- HP, OBA, GO, and MP nodes are left unchanged (phenotype / non-disease ontologies).
- Full report: `output/mondo_normalization.csv`.

## MONDO disease name enrichment

Open Targets also uses MONDO disease ids directly. Before enrichment, the node name was the raw id:

```
id,category,name
MONDO:0002108,biolink:Disease,MONDO:0002108
```

After lookup against MONDO release v2026-06-02 (`mondo_nodes.tsv`):

```
id,category,name
MONDO:0002108,biolink:Disease,thyroid cancer
```

What happened:
- `mondo_nodes.tsv` maps `MONDO:0002108` → `thyroid cancer`.
- After normalization, 824 MONDO disease nodes exist; 823 received labels. One obsolete
  term (`MONDO:0700044`) keeps the raw id. Full report: `output/mondo_name_enrichment.csv`.

## Why this matters

Because both sources resolve the gene to the same NCBIGene id, the drug edge and the disease
edge share a node. That produces a connected path across two independently sourced datasets:

```
DRUG_NAME:TRAMETINIB --negatively_regulates--> NCBIGene:6416 --gene_associated_with_condition--> MONDO:0002026
       (DGIdb)                                   MAP2K2 / MEK2                                  (Open Targets)
```

(Many gene-disease edges in the graph now use MONDO ids with disease names; some EFO ids
without a unique xref, such as EFO:0020865 for this gene, remain as EFO.)

In the full output this one gene connects 8 drugs to 12 diseases, so the graph supports
questions like "which diseases are associated with a gene that this drug inhibits," which is
only possible because the two sources were unified onto the Biolink model.
