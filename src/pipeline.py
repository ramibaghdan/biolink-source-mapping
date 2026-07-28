"""
pipeline.py

Build the Biolink graph from the two sampled sources, in order:

  1  ingest        load and clean the committed samples
  2  map           apply the mapping tables, produce nodes and edges
  3  normalize     rewrite disease ids to MONDO where an xref exists
  4  name          give every disease node a readable name
  5  write         nodes.csv, edges.csv, and the reports

Each step is a function below and run() calls them in that order. Reference files
are parsed in lookups.py, downloaded by fetch.py. No network calls here.

Identifier handling matches mapping_decisions.md:
  - genes unify to a single NCBIGene id across both sources via the HGNC map
    (symbol -> entrez for DGIdb, ensembl -> entrez for Open Targets). This is what
    makes the same gene one node, so drug -> gene -> disease paths connect.
  - disease ids go to CURIE form (EFO_0003950 -> EFO:0003950), then to MONDO when
    a unique xref exists. Unmapped ids stay in their source ontology and get their
    name from OLS4 instead.

Run from the project root:
  python src/pipeline.py --hgnc hgnc_complete_set.txt --mondo-nodes data/mondo/mondo_nodes.tsv
"""

import argparse
import os

import pandas as pd

import lookups

RAW_DIR = "data/raw"
MAP_DIR = "mappings"
OUT_DIR = "output"
DEFAULT_LABEL_CACHE = "data/ontology_labels/labels.tsv"

# The single gene-disease relationship from Open Targets. Verified against the
# installed biolink-model by validate.py.
GENE_DISEASE_PREDICATE = "biolink:gene_associated_with_condition"

# Ontologies whose terms are not diseases. Reported in step 4, never rewritten.
# Splitting the category changes edge semantics, so it is its own decision.
NON_DISEASE_CATEGORY = {
    "GO": "biolink:BiologicalProcessOrActivity",
    "HP": "biolink:PhenotypicFeature",
    "MP": "biolink:PhenotypicFeature",
    "OBA": "biolink:PhenotypicFeature",
}


# ==========================================================================
# step 1  ingest
# ==========================================================================

def _clean_str_cols(df):
    for c in df.columns:
        if df[c].dtype == object:
            df[c] = df[c].astype(str).str.strip()
    return df


def ingest_dgidb(path=None):
    """Drug-gene interactions: gene_symbol, drug_name, interaction_type."""
    path = path or os.path.join(RAW_DIR, "dgidb_sample.csv")
    df = pd.read_csv(path, dtype=str)
    df.columns = [c.strip().lower() for c in df.columns]

    missing = {"gene_symbol", "drug_name"} - set(df.columns)
    if missing:
        raise SystemExit(f"DGIdb sample missing columns: {missing}. Found: {list(df.columns)}")

    df = _clean_str_cols(df)
    # Drop rows that cannot form a triple.
    df = df.dropna(subset=["gene_symbol", "drug_name"])
    df = df[(df["gene_symbol"] != "") & (df["drug_name"] != "")]
    # interaction_type may be absent. Normalize every missing form to other/unknown
    # so those rows map to a general predicate instead of being dropped.
    if "interaction_type" not in df.columns:
        df["interaction_type"] = "other/unknown"
    df["interaction_type"] = df["interaction_type"].replace(
        {"": "other/unknown", "nan": "other/unknown",
         "NaN": "other/unknown", "None": "other/unknown"}
    )
    df["interaction_type"] = df["interaction_type"].fillna("other/unknown")
    return df.drop_duplicates().reset_index(drop=True)


def ingest_opentargets(path=None):
    """Gene-disease associations: ensembl_id, disease_efo_id."""
    path = path or os.path.join(RAW_DIR, "opentargets_sample.csv")
    df = pd.read_csv(path, dtype=str)
    df.columns = [c.strip().lower() for c in df.columns]

    # The gene column arrives under several names depending on the export.
    gene_col = next((c for c in df.columns
                     if c in ("ensembl_id", "targetid", "target_id", "gene_id")), None)
    if gene_col and gene_col != "ensembl_id":
        df = df.rename(columns={gene_col: "ensembl_id"})

    missing = {"ensembl_id", "disease_efo_id"} - set(df.columns)
    if missing:
        raise SystemExit(
            f"Open Targets sample missing columns: {missing}. Found: {list(df.columns)}. "
            f"The gene column may not have been written by the sampler. "
            f"Re-run sample_sources.py, or rename the gene column to 'ensembl_id'."
        )

    df = _clean_str_cols(df)
    df = df.dropna(subset=["ensembl_id", "disease_efo_id"])
    df = df[(df["ensembl_id"] != "") & (df["disease_efo_id"] != "")]
    # One row per gene-disease pair. The sample carries duplicates by datasource.
    df = df.drop_duplicates(subset=["ensembl_id", "disease_efo_id"])
    return df.reset_index(drop=True)


# ==========================================================================
# step 2  map to Biolink
# ==========================================================================

def _disease_curie(raw):
    """EFO_0003950 -> EFO:0003950. Keeps the original ontology prefix."""
    if raw is None:
        return None
    raw = str(raw).strip()
    return raw.replace("_", ":", 1) if "_" in raw and ":" not in raw else raw


def map_to_biolink(dgidb, opentargets, sym_to_entrez, ensg_to_entrez):
    """Apply the mapping tables and return (nodes, edges, stats)."""
    pred = pd.read_csv(os.path.join(MAP_DIR, "predicate_map.csv"), dtype=str)
    pred.columns = [c.strip().lower() for c in pred.columns]
    pred_lookup = dict(zip(pred["source_relationship"].str.lower(),
                           pred["biolink_predicate"]))

    nodes = {}
    edges = []
    unmapped_predicates = set()
    genes_without_id = set()

    def gene_from_symbol(symbol):
        sym = str(symbol).upper().strip()
        entrez = sym_to_entrez.get(sym)
        if entrez:
            gid = f"NCBIGene:{entrez}"
            nodes[gid] = ("biolink:Gene", sym)
            return gid
        genes_without_id.add(sym)
        return None

    def gene_from_ensembl(ensembl_id):
        e = str(ensembl_id).strip()
        entrez = ensg_to_entrez.get(e)
        if entrez:
            gid = f"NCBIGene:{entrez}"
            # Keep the symbol if this gene was already seen, else use the ensembl id.
            if gid not in nodes:
                nodes[gid] = ("biolink:Gene", e)
            return gid
        genes_without_id.add(e)
        return None

    # Drug-gene edges from DGIdb. Drug is subject, gene is object.
    for _, r in dgidb.iterrows():
        g_id = gene_from_symbol(r["gene_symbol"])
        if g_id is None:
            continue
        drug_name = r["drug_name"]
        drug_id = f"DRUG_NAME:{drug_name}"
        itype = str(r.get("interaction_type", "other/unknown")).lower().strip()
        predicate = pred_lookup.get(itype)
        if predicate is None or str(predicate).strip() in ("", "nan"):
            unmapped_predicates.add(itype)
            continue
        nodes[drug_id] = ("biolink:ChemicalEntity", drug_name)
        edges.append({
            "subject": drug_id, "predicate": predicate, "object": g_id,
            "subject_category": "biolink:ChemicalEntity", "object_category": "biolink:Gene",
            "knowledge_source": "DGIdb", "source_relationship": itype,
        })

    # Gene-disease edges from Open Targets. Gene is subject, disease is object.
    for _, r in opentargets.iterrows():
        g_id = gene_from_ensembl(r["ensembl_id"])
        if g_id is None:
            continue
        d_curie = _disease_curie(r["disease_efo_id"])
        nodes[d_curie] = ("biolink:Disease", d_curie)
        edges.append({
            "subject": g_id, "predicate": GENE_DISEASE_PREDICATE, "object": d_curie,
            "subject_category": "biolink:Gene", "object_category": "biolink:Disease",
            "knowledge_source": "OpenTargets", "source_relationship": "associated_with",
        })

    nodes_df = pd.DataFrame([{"id": k, "category": v[0], "name": v[1]}
                             for k, v in nodes.items()])
    edges_df = pd.DataFrame(edges)
    stats = {"unmapped_predicates": sorted(unmapped_predicates),
             "genes_without_id": len(genes_without_id)}
    return nodes_df, edges_df, stats


# ==========================================================================
# step 3  normalize disease ids to MONDO
# ==========================================================================

def _dedupe_nodes(nodes_df):
    """Collapse nodes that became the same id, keeping a real name over an id."""
    rows = []
    for node_id, group in nodes_df.groupby("id", sort=False):
        names = [str(n) for n in group["name"].tolist()]
        rows.append({
            "id": node_id,
            "category": group.iloc[0]["category"],
            "name": next((n for n in names if n != node_id), names[0]),
        })
    return pd.DataFrame(rows)


def normalize_disease_ids(nodes_df, edges_df, xref_map):
    """Rewrite disease ids to MONDO where the xref lookup succeeds.

    Ids outside NORMALIZE_PREFIXES (HP, OBA, GO, MP) and ids with no xref are
    left alone. Returns (nodes, edges, report).
    """
    diseases = nodes_df[nodes_df["category"] == "biolink:Disease"]
    id_map = {}
    report_rows = []

    for source_id in diseases["id"].astype(str):
        prefix = source_id.split(":", 1)[0]
        if prefix == "MONDO":
            target, status = source_id, "already_mondo"
        elif prefix in lookups.NORMALIZE_PREFIXES:
            mondo_id = xref_map.get(source_id)
            target, status = (mondo_id, "mapped") if mondo_id else (source_id, "unmapped")
        else:
            target, status = source_id, "skipped_prefix"

        id_map[source_id] = target
        report_rows.append({
            "source_id": source_id,
            "mondo_id": target if status == "mapped" else ("" if status != "already_mondo" else target),
            "prefix": prefix,
            "normalized": status == "mapped",
            "status": status,
        })

    nodes = nodes_df.copy()
    nodes["id"] = nodes["id"].astype(str).map(lambda x: id_map.get(x, x))
    nodes = _dedupe_nodes(nodes)

    edges = edges_df.copy()
    is_disease = edges["object_category"] == "biolink:Disease"
    edges.loc[is_disease, "object"] = (
        edges.loc[is_disease, "object"].astype(str).map(lambda x: id_map.get(x, x))
    )
    return nodes, edges, pd.DataFrame(report_rows)


# ==========================================================================
# step 4  name the disease nodes
# ==========================================================================

def name_disease_nodes(nodes_df, mondo_labels=None, ontology_labels=None):
    """Give every disease node a readable name, from MONDO first then OLS4.

    MONDO nodes take their label from the release. Whatever is still named by its
    id (EFO, HP, OBA, GO, MP, Orphanet) is named from the OLS4 cache. A node with
    no label anywhere keeps its id rather than gaining a blank name.

    Returns (nodes, report). The report has one row per node either pass touched.
    """
    nodes = nodes_df.copy()
    report_rows = []

    # MONDO release labels.
    if mondo_labels:
        is_mondo = nodes["id"].astype(str).str.startswith("MONDO:")
        for idx in nodes.index[is_mondo]:
            node_id = str(nodes.at[idx, "id"])
            previous = nodes.at[idx, "name"]
            label = mondo_labels.get(node_id)
            if label:
                nodes.at[idx, "name"] = label
            report_rows.append({
                "id": node_id, "prefix": "MONDO", "previous_name": previous,
                "label": label, "label_source": "mondo_release", "label_ontology": "mondo",
                "obsolete_term": None, "matched": bool(label),
                "current_category": nodes.at[idx, "category"], "suggested_category": None,
            })

    # OLS4 labels for whatever MONDO did not cover.
    if ontology_labels:
        for idx, row in nodes.iterrows():
            node_id = str(row["id"]).strip()
            if not lookups.is_unlabeled(node_id, row["name"]):
                continue
            entry = ontology_labels.get(node_id)
            label = entry["label"] if entry else None
            if label:
                nodes.at[idx, "name"] = label
            prefix = lookups.prefix_of(node_id)
            report_rows.append({
                "id": node_id, "prefix": prefix, "previous_name": row["name"],
                "label": label, "label_source": "ols4",
                "label_ontology": entry["ontology"] if entry else None,
                "obsolete_term": entry["obsolete"] if entry else None,
                "matched": bool(label), "current_category": row["category"],
                "suggested_category": NON_DISEASE_CATEGORY.get(prefix),
            })

    return nodes, pd.DataFrame(report_rows)


# ==========================================================================
# step 5  write
# ==========================================================================

def gene_synonyms_table(nodes_df, synonyms):
    """One row per alternate symbol, for the genes actually in the graph.

    Downstream consumers need to know EGFR and ERBB1 are the same node. The
    graph itself only carries the primary symbol, so this ships alongside it.
    """
    rows = []
    genes = nodes_df[nodes_df["category"] == "biolink:Gene"]
    for node_id, primary in zip(genes["id"], genes["name"]):
        entrez = str(node_id).split(":", 1)[-1]
        for alt in synonyms.get(entrez, []):
            rows.append({"id": node_id, "primary_symbol": primary, "synonym": alt})
    return pd.DataFrame(rows, columns=["id", "primary_symbol", "synonym"])


def write_outputs(nodes_df, edges_df, norm_report=None, name_report=None, synonyms=None):
    os.makedirs(OUT_DIR, exist_ok=True)
    nodes_df.to_csv(os.path.join(OUT_DIR, "nodes.csv"), index=False)
    edges_df.to_csv(os.path.join(OUT_DIR, "edges.csv"), index=False)
    if norm_report is not None and len(norm_report):
        norm_report.to_csv(os.path.join(OUT_DIR, "mondo_normalization.csv"), index=False)
    if name_report is not None and len(name_report):
        name_report.to_csv(os.path.join(OUT_DIR, "naming.csv"), index=False)
    if synonyms is not None and len(synonyms):
        synonyms.to_csv(os.path.join(OUT_DIR, "gene_synonyms.csv"), index=False)


# ==========================================================================
# run
# ==========================================================================

def run(hgnc_path, mondo_nodes_path=None, label_cache=DEFAULT_LABEL_CACHE):
    sym_to_entrez, ensg_to_entrez = lookups.load_hgnc(hgnc_path)
    gene_synonyms = lookups.load_gene_synonyms(hgnc_path)

    # 1  ingest
    dgidb = ingest_dgidb()
    opentargets = ingest_opentargets()
    print(f"ingested: {len(dgidb)} DGIdb rows, {len(opentargets)} Open Targets rows")

    # 2  map
    nodes, edges, stats = map_to_biolink(dgidb, opentargets, sym_to_entrez, ensg_to_entrez)
    print(f"mapped:   {len(nodes)} nodes, {len(edges)} edges")
    if stats["unmapped_predicates"]:
        print(f"  unmapped interaction types (skipped): {stats['unmapped_predicates']}")
    if stats["genes_without_id"]:
        print(f"  genes with no NCBIGene id (dropped): {stats['genes_without_id']}")

    # 3  normalize
    norm_report = None
    mondo_labels = None
    if mondo_nodes_path:
        xref_map, ambiguous = lookups.load_mondo_xrefs(mondo_nodes_path)
        nodes, edges, norm_report = normalize_disease_ids(nodes, edges, xref_map)
        candidates = norm_report[norm_report["status"].isin({"mapped", "unmapped"})]
        print(f"normalized: {int(norm_report['normalized'].sum())}/{len(candidates)} "
              f"EFO/DOID/Orphanet/OTAR ids -> MONDO")
        if ambiguous:
            print(f"  ambiguous xrefs excluded from the map: {len(ambiguous)}")
        mondo_labels = lookups.load_mondo_labels(mondo_nodes_path)

    # 4  name
    ontology_labels = None
    if label_cache and os.path.exists(label_cache):
        ontology_labels = lookups.load_ontology_labels(label_cache)
    elif label_cache:
        print(f"note: no label cache at {label_cache}. run: python src/fetch.py labels")

    nodes, name_report = name_disease_nodes(nodes, mondo_labels, ontology_labels)
    if len(name_report):
        for source, group in name_report.groupby("label_source"):
            print(f"named ({source}): {int(group['matched'].sum())}/{len(group)}")
        # Check the nodes, not the report. A node the release pass missed can still
        # be named by the OLS4 pass, and it has one report row per pass.
        unresolved = lookups.find_unlabeled(nodes)
        if unresolved:
            shown = ", ".join(unresolved[:5])
            more = f" +{len(unresolved) - 5} more" if len(unresolved) > 5 else ""
            print(f"  still named by id: {shown}{more}")
        odd = name_report[name_report["suggested_category"].notna()]
        if len(odd):
            print(f"  {len(odd)} nodes typed biolink:Disease are not diseases "
                  f"({odd['prefix'].value_counts().to_dict()}), category unchanged")

    # 5  write
    syn_table = gene_synonyms_table(nodes, gene_synonyms)
    write_outputs(nodes, edges, norm_report, name_report, syn_table)
    if len(syn_table):
        print(f"gene synonyms: {len(syn_table)} across "
              f"{syn_table['id'].nunique()} genes")

    drug_genes = set(edges[edges.knowledge_source == "DGIdb"]["object"])
    disease_genes = set(edges[edges.knowledge_source == "OpenTargets"]["subject"])
    print(f"wrote:    {len(nodes)} nodes, {len(edges)} edges to {OUT_DIR}/")
    print(f"genes connecting drug -> gene -> disease: {len(drug_genes & disease_genes)}")
    return nodes, edges


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description="Build the Biolink graph")
    ap.add_argument("--hgnc", required=True, help="path to the HGNC complete set")
    ap.add_argument("--mondo-nodes", default=None,
                    help="mondo_nodes.tsv, enables MONDO normalization and naming")
    ap.add_argument("--labels", default=DEFAULT_LABEL_CACHE,
                    help="OLS4 label cache for ids MONDO does not cover")
    args = ap.parse_args()
    run(args.hgnc, args.mondo_nodes, args.labels)
