"""
check_labels.py

Verify that no node is still displaying its identifier where a name belongs.

Makes "every disease node is readable" a check rather than something eyeballed
in the demo figure. Exits 1 when nodes remain unnamed, so it can gate a build.

    python src/check_labels.py
    python src/check_labels.py --allow 6      # tolerate a documented remainder
"""

import argparse
import collections
import sys

import pandas as pd

from ontology_labels import find_unlabeled, prefix_of


def main():
    ap = argparse.ArgumentParser(description="Check for nodes still named by id")
    ap.add_argument("--nodes", default="output/nodes.csv")
    ap.add_argument(
        "--category",
        default="biolink:Disease",
        help="category to check, or 'all' for every node",
    )
    ap.add_argument(
        "--allow",
        type=int,
        default=0,
        help="pass even with this many unnamed nodes (document them first)",
    )
    args = ap.parse_args()

    nodes = pd.read_csv(args.nodes, dtype=str)
    category = None if args.category == "all" else args.category
    total = len(nodes if category is None else nodes[nodes["category"] == category])

    unlabeled = find_unlabeled(nodes, category=category)
    named = total - len(unlabeled)
    pct = (named / total * 100) if total else 100.0

    label = category or "all"
    print(f"{label}: {named}/{total} named ({pct:.1f}%)")

    if not unlabeled:
        print("PASS - no node is named by its id")
        return 0

    by_prefix = collections.Counter(prefix_of(c) for c in unlabeled)
    print(f"\n{len(unlabeled)} still named by id:")
    for prefix, n in by_prefix.most_common():
        print(f"  {prefix:10} {n}")
    print("\n  " + ", ".join(unlabeled[:10]) + (" ..." if len(unlabeled) > 10 else ""))

    if len(unlabeled) <= args.allow:
        print(f"\nPASS - within the allowed remainder of {args.allow}")
        return 0

    print("\nFAIL - run src/fetch_ontology_labels.py then src/ontology_labels.py")
    print("       if some ids genuinely have no label, document them in")
    print("       mappings/mapping_decisions.md and re-run with --allow N")
    return 1


if __name__ == "__main__":
    sys.exit(main())
