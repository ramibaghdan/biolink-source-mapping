"""
ontology_labels.py

Fill in human-readable names for disease nodes whose name is still a bare CURIE.

The MONDO pass (mondo_normalize.py + mondo_labels.py) names every node that
normalized to MONDO. Ids that had no unique MONDO xref stay in their source
ontology and keep the id as their name, which is unreadable downstream:

    EFO:0006890,biolink:Disease,EFO:0006890      <- before
    EFO:0006890,biolink:Disease,gestational age  <- after

This module applies labels collected by fetch_ontology_labels.py (OLS4) for the
ontologies MONDO does not cover here: EFO, HP, OBA, GO, MP, Orphanet.

No network calls. Run fetch_ontology_labels.py first to build the cache.

Also reports, but does not change, a category observation: several ids reaching
this stage are not diseases at all (GO biological processes, HP/MP phenotypes,
OBA attributes) yet carry biolink:Disease from the Open Targets ingest. The
suggested_category column exists so that judgment call can be made explicitly
rather than silently. See mappings/mapping_decisions.md.
"""

import os
import re
import pandas as pd

DEFAULT_CACHE = "data/ontology_labels/labels.tsv"

# A name that is just the id (or an id-shaped string) means "never labelled".
_CURIE_SHAPED = re.compile(r"^[A-Za-z][A-Za-z0-9]*[:_]\d+$")

# Ontologies whose terms are not diseases, for the report only. Nothing is
# rewritten here; this flags the mismatch so the decision can be documented.
NON_DISEASE_ONTOLOGY_CATEGORY = {
    "GO": "biolink:BiologicalProcessOrActivity",
    "HP": "biolink:PhenotypicFeature",
    "MP": "biolink:PhenotypicFeature",
    "OBA": "biolink:PhenotypicFeature",
}


def is_unlabeled(node_id, name):
    """True when a node never received a human-readable name."""
    if name is None or (isinstance(name, float) and pd.isna(name)):
        return True
    name = str(name).strip()
    if not name:
        return True
    if name == str(node_id).strip():
        return True
    return bool(_CURIE_SHAPED.match(name))


def find_unlabeled(nodes_df, category="biolink:Disease"):
    """Return the sorted list of CURIEs still showing an id where a name should be."""
    df = nodes_df
    if category:
        df = df[df["category"] == category]
    missing = [
        str(nid).strip()
        for nid, name in zip(df["id"], df["name"])
        if is_unlabeled(nid, name)
    ]
    return sorted(set(missing))


def prefix_of(curie):
    return str(curie).split(":", 1)[0] if ":" in str(curie) else ""


def load_ontology_labels(cache_path=DEFAULT_CACHE):
    """Return {'EFO:0006890': {'label': ..., 'ontology': ..., 'obsolete': bool}, ...}.

    Reads the TSV written by fetch_ontology_labels.py. Rows with an empty label
    are kept out of the map: they represent ids OLS4 could not resolve, and the
    node should keep its id rather than gain a blank name.
    """
    if not os.path.exists(cache_path):
        raise SystemExit(
            f"label cache not found: {cache_path}\n"
            "run: python src/fetch_ontology_labels.py"
        )
    df = pd.read_csv(cache_path, sep="\t", dtype=str).fillna("")
    labels = {}
    for _, r in df.iterrows():
        label = str(r.get("label", "")).strip()
        if not label:
            continue
        labels[str(r["curie"]).strip()] = {
            "label": label,
            "ontology": str(r.get("ontology", "")).strip(),
            "obsolete": str(r.get("obsolete", "")).strip().lower() == "true",
        }
    return labels


def apply_ontology_labels(nodes_df, labels):
    """Name every still-unlabeled node we have a label for.

    Only touches rows that failed is_unlabeled(); nodes already named by the
    MONDO pass are left exactly as they are. Returns (nodes_df, report).
    """
    nodes = nodes_df.copy()
    report_rows = []

    for idx, row in nodes.iterrows():
        nid = str(row["id"]).strip()
        if not is_unlabeled(nid, row["name"]):
            continue

        entry = labels.get(nid)
        prefix = prefix_of(nid)
        new_name = entry["label"] if entry else None
        if new_name:
            nodes.at[idx, "name"] = new_name

        report_rows.append({
            "id": nid,
            "prefix": prefix,
            "previous_name": row["name"],
            "label": new_name,
            "label_ontology": entry["ontology"] if entry else None,
            "obsolete_term": entry["obsolete"] if entry else None,
            "matched": bool(new_name),
            "current_category": row["category"],
            "suggested_category": NON_DISEASE_ONTOLOGY_CATEGORY.get(prefix),
        })

    report = pd.DataFrame(report_rows)
    return nodes, report


def main():
    import argparse

    ap = argparse.ArgumentParser(
        description="Apply OLS4-derived labels to disease nodes still named by id"
    )
    ap.add_argument("--nodes", default="output/nodes.csv")
    ap.add_argument("--cache", default=DEFAULT_CACHE)
    ap.add_argument("--report", default="output/ontology_name_enrichment.csv")
    ap.add_argument(
        "--dry-run",
        action="store_true",
        help="report what would change without writing nodes.csv",
    )
    args = ap.parse_args()

    nodes_df = pd.read_csv(args.nodes, dtype=str)
    before = find_unlabeled(nodes_df)
    if not before:
        print("all disease nodes already have names; nothing to do")
        return

    labels = load_ontology_labels(args.cache)
    nodes_df, report = apply_ontology_labels(nodes_df, labels)

    matched = int(report["matched"].sum())
    total = len(report)
    print(f"ontology labels applied: {matched}/{total} previously unnamed nodes")

    by_prefix = report.groupby("prefix")["matched"].agg(["sum", "count"])
    for prefix, r in by_prefix.iterrows():
        print(f"  {prefix:10} {int(r['sum'])}/{int(r['count'])}")

    unresolved = report[~report["matched"]]["id"].tolist()
    if unresolved:
        shown = unresolved[:5]
        more = f" ... +{len(unresolved) - 5} more" if len(unresolved) > 5 else ""
        print(f"still unresolved (kept as id): {shown}{more}")

    obsolete = report[report["obsolete_term"] == True]  # noqa: E712
    if len(obsolete):
        print(f"labels from obsolete terms: {len(obsolete)} (see report)")

    miscategorized = report[report["suggested_category"].notna()]
    if len(miscategorized):
        counts = miscategorized["prefix"].value_counts().to_dict()
        print(
            f"note: {len(miscategorized)} nodes typed biolink:Disease are not diseases "
            f"({counts}); category unchanged, see suggested_category in the report"
        )

    os.makedirs(os.path.dirname(args.report) or ".", exist_ok=True)
    report.to_csv(args.report, index=False)
    print(f"report -> {args.report}")

    if args.dry_run:
        print("dry run: nodes.csv not written")
        return

    nodes_df.to_csv(args.nodes, index=False)
    print(f"nodes  -> {args.nodes}")


if __name__ == "__main__":
    main()
