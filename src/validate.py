"""
validate.py

Check the generated graph and print a summary. Two things are checked:

  biolink terms     every category and predicate is a real Biolink term
  naming            no node is left showing its identifier where a name belongs

Exits 1 on failure so it can gate a build. If the biolink-model toolkit is
installed the term check is authoritative, otherwise it falls back to a vendored
list of the terms this project uses and says so.

  python src/validate.py
  python src/validate.py --allow 6      tolerate a documented naming remainder
"""

import argparse
import collections
import os
import sys

import pandas as pd

import lookups

OUT_DIR = "output"

# The terms this project maps to. Only here so validate runs without a network.
# The authoritative source is the biolink-model toolkit.
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
    """Return (categories, predicates, source_label). Try the real toolkit first."""
    try:
        from bmt import Toolkit
        tk = Toolkit()
        return (set(tk.get_all_classes(formatted=True)),
                set(tk.get_all_predicates(formatted=True)),
                "biolink-model toolkit (authoritative)")
    except Exception:
        return (FALLBACK_CATEGORIES, FALLBACK_PREDICATES,
                "vendored fallback (install bmt for the authoritative check)")


def check_biolink_terms(nodes, edges):
    """Flag categories and predicates that are not real Biolink terms."""
    cats, preds, source = get_biolink_terms()
    used_cats = set(nodes["category"].dropna())
    used_preds = set(edges["predicate"].dropna())
    bad_cats = sorted(c for c in used_cats if c not in cats)
    bad_preds = sorted(p for p in used_preds if p not in preds)

    print(f"validating against: {source}\n")
    print("=== summary ===")
    print(f"nodes: {len(nodes)}")
    print(f"edges: {len(edges)}")
    print(f"distinct categories: {len(used_cats)} -> {sorted(used_cats)}")
    print(f"distinct predicates: {len(used_preds)} -> {sorted(used_preds)}")
    print()

    print(f"CATEGORIES NOT IN BIOLINK: {bad_cats}" if bad_cats else "all categories valid")
    print(f"PREDICATES NOT IN BIOLINK: {bad_preds}" if bad_preds else "all predicates valid")

    # Gene ids that never resolved to NCBIGene.
    unresolved = nodes[nodes["id"].astype(str).str.startswith(("GENE_SYMBOL:", "ENSEMBL:"))]
    if len(unresolved):
        print(f"note: {len(unresolved)} gene nodes did not resolve to NCBIGene")

    return not bad_cats and not bad_preds


def check_naming(nodes, allow=0, category="biolink:Disease"):
    """Flag nodes still showing their identifier instead of a name."""
    subset = nodes if category is None else nodes[nodes["category"] == category]
    unlabeled = lookups.find_unlabeled(nodes, category=category)
    named = len(subset) - len(unlabeled)
    pct = (named / len(subset) * 100) if len(subset) else 100.0

    print(f"\n=== naming ===")
    print(f"{category or 'all'}: {named}/{len(subset)} named ({pct:.1f}%)")

    if not unlabeled:
        print("no node is named by its id")
        return True

    by_prefix = collections.Counter(lookups.prefix_of(c) for c in unlabeled)
    print(f"{len(unlabeled)} still named by id:")
    for prefix, n in by_prefix.most_common():
        print(f"  {prefix:10} {n}")
    print("  " + ", ".join(unlabeled[:10]) + (" ..." if len(unlabeled) > 10 else ""))

    if len(unlabeled) <= allow:
        print(f"within the allowed remainder of {allow}")
        return True

    print("run: python src/fetch.py labels, then rebuild the pipeline")
    print("if some ids genuinely have no label, document them in")
    print("mappings/mapping_decisions.md and re-run with --allow N")
    return False


def main():
    ap = argparse.ArgumentParser(description="Validate the generated graph")
    ap.add_argument("--allow", type=int, default=0,
                    help="pass with this many nodes still named by id")
    ap.add_argument("--category", default="biolink:Disease",
                    help="category to name-check, or 'all'")
    args = ap.parse_args()

    nodes = pd.read_csv(os.path.join(OUT_DIR, "nodes.csv"))
    edges = pd.read_csv(os.path.join(OUT_DIR, "edges.csv"))

    terms_ok = check_biolink_terms(nodes, edges)
    naming_ok = check_naming(nodes, args.allow,
                             None if args.category == "all" else args.category)

    ok = terms_ok and naming_ok
    print(f"\nvalidation {'passed' if ok else 'found issues (see above)'}")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
