"""
map_to_biolink.py

Apply the mapping tables to the cleaned source data and produce Biolink nodes and edges.

Outputs:
  output/nodes.csv  one row per unique entity (id, category, name)
  output/edges.csv  one row per association (subject, predicate, object, categories, source)

Identifier handling matches mapping_decisions.md:
  - genes are unified to a single NCBIGene id across both sources, using the HGNC map
    (symbol -> entrez for DGIdb, ensembl -> entrez for Open Targets). This is what makes
    the same gene one node, so drug -> gene -> disease paths actually connect.
  - disease ids are converted from underscore to colon CURIE form (EFO_0003950 -> EFO:0003950)
    and, when --mondo-nodes is supplied, EFO/DOID/Orphanet/OTAR ids with unique MONDO xrefs
    are normalized to MONDO; unmapped ids stay in their original ontology.
No network calls.
"""

import os
import argparse
import pandas as pd

import ingest
from mondo_labels import apply_mondo_labels, load_mondo_labels
from mondo_normalize import apply_mondo_normalization, load_mondo_xrefs

MAP_DIR = "mappings"
OUT_DIR = "output"

# The single gene-disease relationship from Open Targets. Verify the exact slot name
# against the installed biolink-model at validate time.
GENE_DISEASE_PREDICATE = "biolink:gene_associated_with_condition"


def load_maps():
    pred = pd.read_csv(os.path.join(MAP_DIR, "predicate_map.csv"), dtype=str)
    pred.columns = [c.strip().lower() for c in pred.columns]
    return pred


def load_hgnc(hgnc_path):
    """Build symbol->entrez and ensembl->entrez from the HGNC complete set.

    Both sources are unified to NCBIGene this way: DGIdb symbols and Open Targets
    Ensembl ids both resolve to the same Entrez id, so one gene is one node.
    """
    h = pd.read_csv(hgnc_path, sep="\t", dtype=str, low_memory=False)
    h.columns = [c.strip().lower() for c in h.columns]
    sym = next((c for c in h.columns if c == "symbol"), None)
    entrez = next((c for c in h.columns if "entrez" in c), None)
    ensg = next((c for c in h.columns if "ensembl" in c), None)
    if not (sym and entrez and ensg):
        raise SystemExit(f"HGNC needs symbol+entrez+ensembl. found: {list(h.columns)}")
    h = h[[sym, entrez, ensg]].dropna(subset=[entrez])
    sym_map = {str(s).upper().strip(): str(e).strip()
               for s, e in zip(h[sym], h[entrez]) if pd.notna(s)}
    ensg_map = {str(g).strip(): str(e).strip()
                for g, e in zip(h[ensg], h[entrez]) if pd.notna(g)}
    return sym_map, ensg_map


def normalize_disease_curie(raw):
    """EFO_0003950 -> EFO:0003950. Keep the original ontology prefix."""
    if raw is None:
        return None
    raw = str(raw).strip()
    if "_" in raw and ":" not in raw:
        return raw.replace("_", ":", 1)
    return raw


def build(hgnc_path, mondo_nodes_path=None):
    os.makedirs(OUT_DIR, exist_ok=True)
    pred_map = load_maps()
    dgidb = ingest.load_dgidb()
    ot = ingest.load_opentargets()

    pred_lookup = dict(zip(pred_map["source_relationship"].str.lower(),
                           pred_map["biolink_predicate"]))
    sym_to_entrez, ensg_to_entrez = load_hgnc(hgnc_path)

    nodes = {}
    edges = []
    unmapped_predicates = set()
    genes_without_id = set()

    def gene_node_from_symbol(symbol):
        sym = str(symbol).upper().strip()
        entrez = sym_to_entrez.get(sym)
        if entrez:
            gid = f"NCBIGene:{entrez}"
            nodes[gid] = ("biolink:Gene", sym)
            return gid
        genes_without_id.add(sym)
        return None

    def gene_node_from_ensembl(ensembl_id):
        e = str(ensembl_id).strip()
        entrez = ensg_to_entrez.get(e)
        if entrez:
            gid = f"NCBIGene:{entrez}"
            # keep the symbol name if we already saw this gene; else use ensembl as name
            if gid not in nodes:
                nodes[gid] = ("biolink:Gene", e)
            return gid
        genes_without_id.add(e)
        return None

    # Drug-gene edges from DGIdb (drug is subject, gene is object).
    for _, r in dgidb.iterrows():
        g_id = gene_node_from_symbol(r["gene_symbol"])
        if g_id is None:
            continue
        drug_name = r["drug_name"]
        drug_id = f"DRUG_NAME:{drug_name}"
        itype = str(r.get("interaction_type", "other/unknown")).lower().strip()
        predicate = pred_lookup.get(itype)
        if predicate is None or str(predicate).strip() == "" or str(predicate) == "nan":
            unmapped_predicates.add(itype)
            continue
        nodes[drug_id] = ("biolink:ChemicalEntity", drug_name)
        edges.append({
            "subject": drug_id, "predicate": predicate, "object": g_id,
            "subject_category": "biolink:ChemicalEntity", "object_category": "biolink:Gene",
            "knowledge_source": "DGIdb", "source_relationship": itype,
        })

    # Gene-disease edges from Open Targets (gene is subject, disease is object).
    for _, r in ot.iterrows():
        g_id = gene_node_from_ensembl(r["ensembl_id"])
        if g_id is None:
            continue
        d_curie = normalize_disease_curie(r["disease_efo_id"])
        nodes[d_curie] = ("biolink:Disease", d_curie)
        edges.append({
            "subject": g_id, "predicate": GENE_DISEASE_PREDICATE, "object": d_curie,
            "subject_category": "biolink:Gene", "object_category": "biolink:Disease",
            "knowledge_source": "OpenTargets", "source_relationship": "associated_with",
        })

    nodes_df = pd.DataFrame([{"id": k, "category": v[0], "name": v[1]}
                             for k, v in nodes.items()])
    edges_df = pd.DataFrame(edges)

    if mondo_nodes_path:
        xref_map, ambiguous_xrefs = load_mondo_xrefs(mondo_nodes_path)
        nodes_df, edges_df, norm_report = apply_mondo_normalization(
            nodes_df, edges_df, xref_map
        )
        norm_report.to_csv(os.path.join(OUT_DIR, "mondo_normalization.csv"), index=False)
        mapped = int(norm_report["normalized"].sum())
        candidates = norm_report[norm_report["status"].isin({"mapped", "unmapped"})]
        print(f"mondo normalization: {mapped}/{len(candidates)} EFO/DOID/Orphanet/OTAR ids -> MONDO")
        if ambiguous_xrefs:
            print(f"  ambiguous xrefs excluded from map: {len(ambiguous_xrefs)}")

        mondo_labels = load_mondo_labels(mondo_nodes_path)
        nodes_df, label_report = apply_mondo_labels(nodes_df, mondo_labels)
        label_report.to_csv(os.path.join(OUT_DIR, "mondo_name_enrichment.csv"), index=False)
        matched = int(label_report["matched"].sum())
        total = len(label_report)
        print(f"mondo labels applied: {matched}/{total} MONDO disease nodes named")
        unmatched = label_report[~label_report["matched"]]["mondo_id"].tolist()
        if unmatched:
            print(f"mondo ids without label (kept as id): {unmatched[:5]}"
                  + (f" ... +{len(unmatched)-5} more" if len(unmatched) > 5 else ""))
            print("  (often obsolete terms excluded from naming, or ids absent from the release)")

    nodes_df.to_csv(os.path.join(OUT_DIR, "nodes.csv"), index=False)
    edges_df.to_csv(os.path.join(OUT_DIR, "edges.csv"), index=False)

    print(f"nodes: {len(nodes_df)}")
    print(f"edges: {len(edges_df)}")
    if unmapped_predicates:
        print(f"unmapped interaction types (skipped): {sorted(unmapped_predicates)}")
    if genes_without_id:
        print(f"genes with no NCBIGene id (dropped): {len(genes_without_id)}")

    # Connectivity check: genes that have both a drug edge and a disease edge.
    drug_genes = set(edges_df[edges_df.knowledge_source == "DGIdb"]["object"])
    dis_genes = set(edges_df[edges_df.knowledge_source == "OpenTargets"]["subject"])
    connected = drug_genes & dis_genes
    print(f"genes connecting drug -> gene -> disease: {len(connected)}")
    return nodes_df, edges_df


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--hgnc", required=True, help="path to HGNC complete set")
    ap.add_argument(
        "--mondo-nodes",
        default=None,
        help="path to mondo_nodes.tsv from a MONDO release (normalize to MONDO + enrich names)",
    )
    args = ap.parse_args()
    build(hgnc_path=args.hgnc, mondo_nodes_path=args.mondo_nodes)
