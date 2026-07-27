"""
fetch.py

Download the reference files the pipeline needs. The only module that touches the
network. Run once, then everything else works offline.

  python src/fetch.py mondo     mondo_nodes.tsv from a pinned MONDO release (~23 MB, gitignored)
  python src/fetch.py labels    OLS4 labels for the ids MONDO does not cover (~30 KB, committed)

The MONDO file is large and pinned to a release, so it stays out of git. The label
cache only covers ids this graph actually contains, so it is small enough to commit
and keeps the naming step reproducible without a network.
"""

import argparse
import csv
import json
import os
import sys
import time
import urllib.error
import urllib.parse
import urllib.request

import pandas as pd

import lookups

OLS4_TERMS_URL = "https://www.ebi.ac.uk/ols4/api/terms"
DEFAULT_MONDO_OUT = "data/mondo/mondo_nodes.tsv"
DEFAULT_LABELS_OUT = "data/ontology_labels/labels.tsv"
DEFAULT_NODES = "output/nodes.csv"
LABEL_FIELDS = ["curie", "label", "ontology", "obsolete", "iri", "status"]


# --------------------------------------------------------------------------
# mondo
# --------------------------------------------------------------------------

def fetch_mondo(release, dest):
    url = ("https://github.com/monarch-initiative/mondo/releases/download/"
           f"{release}/mondo_nodes.tsv")
    os.makedirs(os.path.dirname(dest) or ".", exist_ok=True)
    print(f"downloading {url}")
    print(f"         -> {dest}")
    urllib.request.urlretrieve(url, dest)
    print(f"done ({os.path.getsize(dest) / (1024 * 1024):.1f} MB)")


# --------------------------------------------------------------------------
# ols4 labels
# --------------------------------------------------------------------------

def _read_cache(path):
    if not os.path.exists(path):
        return {}
    with open(path, newline="", encoding="utf-8") as fh:
        return {r["curie"]: r for r in csv.DictReader(fh, delimiter="\t")}


def _write_cache(path, records):
    os.makedirs(os.path.dirname(path) or ".", exist_ok=True)
    with open(path, "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=LABEL_FIELDS, delimiter="\t")
        w.writeheader()
        for curie in sorted(records):
            w.writerow({k: records[curie].get(k, "") for k in LABEL_FIELDS})


def _get_json(url, timeout, retries, backoff):
    """GET with retries. Returns parsed JSON, or None on 404 or persistent failure."""
    last = None
    for attempt in range(retries + 1):
        try:
            req = urllib.request.Request(url, headers={
                "Accept": "application/json",
                "User-Agent": "biolink-source-mapping/fetch",
            })
            with urllib.request.urlopen(req, timeout=timeout) as resp:
                return json.loads(resp.read().decode("utf-8"))
        except urllib.error.HTTPError as e:
            if e.code == 404:
                return None
            last = e
        except Exception as e:
            last = e
        if attempt < retries:
            time.sleep(backoff * (2 ** attempt))
    print(f"    giving up: {last}", file=sys.stderr)
    return None


def _pick_term(terms):
    """Choose one term when several ontologies return the same obo_id.

    Rank non-obsolete first, then the defining ontology. OLS names obsolete terms
    "obsolete <name>", so a live label from an importing ontology is the better
    display name even though the defining ontology is more authoritative.
    """
    labelled = [t for t in terms if (t.get("label") or "").strip()]
    if not labelled:
        return None
    return sorted(labelled, key=lambda t: (
        1 if t.get("is_obsolete") else 0,
        0 if t.get("is_defining_ontology") else 1,
    ))[0]


def resolve(curie, timeout, retries, backoff):
    """Look up one CURIE. status is ok or not_found."""
    url = f"{OLS4_TERMS_URL}?{urllib.parse.urlencode({'obo_id': curie, 'size': 20})}"
    data = _get_json(url, timeout, retries, backoff)
    term = _pick_term(((data or {}).get("_embedded") or {}).get("terms") or [])
    if not term:
        return {"curie": curie, "label": "", "ontology": "", "obsolete": "",
                "iri": "", "status": "not_found"}
    return {
        "curie": curie,
        "label": (term.get("label") or "").strip(),
        "ontology": (term.get("ontology_name") or "").strip(),
        "obsolete": "true" if term.get("is_obsolete") else "false",
        "iri": (term.get("iri") or "").strip(),
        "status": "ok",
    }


def fetch_labels(nodes_path, out, refresh, sleep, timeout, retries, backoff, limit):
    nodes = pd.read_csv(nodes_path, dtype=str)
    curies = lookups.find_unlabeled(nodes)
    if not curies:
        print("no unnamed disease nodes. nothing to fetch")
        return

    cache = {} if refresh else _read_cache(out)
    todo = [c for c in curies if c not in cache]
    if limit:
        todo = todo[:limit]

    print(f"unnamed disease nodes: {len(curies)}")
    print(f"already cached:        {len(curies) - len([c for c in curies if c not in cache])}")
    print(f"to fetch:              {len(todo)}")
    if not todo:
        print("cache is complete")
        return

    ok = 0
    for i, curie in enumerate(todo, 1):
        rec = resolve(curie, timeout, retries, backoff)
        cache[curie] = rec
        ok += rec["status"] == "ok"
        print(f"  [{i}/{len(todo)}] {curie:20} {rec['label'] or '(not found)'}")
        # Save as we go so an interrupted run loses nothing.
        if i % 25 == 0:
            _write_cache(out, cache)
        time.sleep(sleep)

    _write_cache(out, cache)
    print(f"\nresolved {ok}/{len(todo)}")
    print(f"cache -> {out}")


# --------------------------------------------------------------------------

def main():
    ap = argparse.ArgumentParser(description="Download reference files")
    sub = ap.add_subparsers(dest="what", required=True)

    m = sub.add_parser("mondo", help="mondo_nodes.tsv from a MONDO release")
    m.add_argument("--release", default=lookups.MONDO_RELEASE)
    m.add_argument("--out", default=DEFAULT_MONDO_OUT)

    l = sub.add_parser("labels", help="OLS4 labels for ids MONDO does not cover")
    l.add_argument("--nodes", default=DEFAULT_NODES)
    l.add_argument("--out", default=DEFAULT_LABELS_OUT)
    l.add_argument("--refresh", action="store_true", help="re-query ids already cached")
    l.add_argument("--limit", type=int, default=None, help="stop after N lookups")
    l.add_argument("--sleep", type=float, default=0.15)
    l.add_argument("--timeout", type=float, default=30.0)
    l.add_argument("--retries", type=int, default=3)
    l.add_argument("--backoff", type=float, default=1.0)

    args = ap.parse_args()
    if args.what == "mondo":
        fetch_mondo(args.release, args.out)
    else:
        fetch_labels(args.nodes, args.out, args.refresh, args.sleep,
                     args.timeout, args.retries, args.backoff, args.limit)


if __name__ == "__main__":
    main()
