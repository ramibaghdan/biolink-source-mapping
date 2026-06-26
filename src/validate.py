"""
validate.py

Confirm the produced nodes and edges use real Biolink categories and predicates, and
print a summary. If the biolink-model package is installed, validate against it. If not,
fall back to a small vendored list of the terms this project uses and warn that the
authoritative check was skipped.
"""

import os
import pandas as pd

OUT_DIR = "output"

# Vendored fallback list: the terms this project maps to. The authoritative source is
# the biolink-model package; this list only exists so validate runs without a network.
FALLBACK_CATEGORIES = {
    "biolink:Gene", "biolink:ChemicalEntity", "biolink:Disease",
}
FALLBACK_PREDICATES = {
    "biolink:negatively_regulates", "biolink:positively_regulates",
    "biolink:regulates", "biolink:affects", "biolink:physically_interacts_with",
    "biolink:interacts_with", "biolink:gene_associated_with_condition",
    "biolink:associated_with",
}


def get_biolink_terms():
    """Return (categories, predicates, source_label). Try the real package first."""
    try:
        from biolink_model.datamodel.model import ClassDefinition  # noqa
        # The package layout varies by version; rather than depend on internals,
        # use the toolkit if present.
        from bmt import Toolkit
        tk = Toolkit()
        cats = set(tk.get_all_classes(formatted=True))
        preds = set(tk.get_all_predicates(formatted=True))
        return cats, preds, "biolink-model toolkit (authoritative)"
    except Exception:
        return FALLBACK_CATEGORIES, FALLBACK_PREDICATES, "vendored fallback (install bmt for the authoritative check)"


def main():
    nodes = pd.read_csv(os.path.join(OUT_DIR, "nodes.csv"))
    edges = pd.read_csv(os.path.join(OUT_DIR, "edges.csv"))
    cats, preds, src = get_biolink_terms()

    print(f"validating against: {src}\n")

    used_cats = set(nodes["category"].dropna())
    used_preds = set(edges["predicate"].dropna())

    bad_cats = sorted(c for c in used_cats if c not in cats)
    bad_preds = sorted(p for p in used_preds if p not in preds)

    print("=== summary ===")
    print(f"nodes: {len(nodes)}")
    print(f"edges: {len(edges)}")
    print(f"distinct categories: {len(used_cats)} -> {sorted(used_cats)}")
    print(f"distinct predicates: {len(used_preds)} -> {sorted(used_preds)}")
    print()

    if bad_cats:
        print(f"CATEGORIES NOT FOUND IN BIOLINK: {bad_cats}")
    else:
        print("all categories valid")
    if bad_preds:
        print(f"PREDICATES NOT FOUND IN BIOLINK: {bad_preds}")
    else:
        print("all predicates valid")

    # id prefix sanity: flag fallback gene ids that did not resolve.
    unresolved = nodes[nodes["id"].astype(str).str.startswith(("GENE_SYMBOL:", "ENSEMBL:"))]
    if len(unresolved):
        print(f"\nnote: {len(unresolved)} gene nodes did not resolve to NCBIGene")

    ok = not bad_cats and not bad_preds
    print(f"\nvalidation {'passed' if ok else 'found issues (see above)'}")


if __name__ == "__main__":
    main()
