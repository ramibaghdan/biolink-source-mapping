"""
fetch_ontology_labels.py

Resolve the disease-node CURIEs that the MONDO pass left unnamed to
human-readable labels, using the EBI Ontology Lookup Service (OLS4).

One resolver for every ontology still present after MONDO normalization
(EFO, HP, OBA, GO, MP, Orphanet) instead of a downloader per ontology.

Writes a small TSV cache so the labelling step itself needs no network and the
pipeline stays reproducible offline, mirroring how fetch_mondo_release.py pins a
MONDO release. The cache is small enough to commit.

Resumable: ids already in the cache are skipped, including ids previously
recorded as unresolved. Use --refresh to re-query everything.

No dependencies beyond the standard library.

    python src/fetch_ontology_labels.py
    python src/fetch_ontology_labels.py --nodes output/nodes.csv --out data/ontology_labels/labels.tsv
"""

import argparse
import csv
import json
import os
import re
import sys
import time
import urllib.error
import urllib.parse
import urllib.request

OLS4_TERMS_URL = "https://www.ebi.ac.uk/ols4/api/terms"
DEFAULT_OUT = "data/ontology_labels/labels.tsv"
DEFAULT_NODES = "output/nodes.csv"
FIELDS = ["curie", "label", "ontology", "obsolete", "iri", "status"]

_CURIE_SHAPED = re.compile(r"^[A-Za-z][A-Za-z0-9]*[:_]\d+$")


def is_unlabeled(node_id, name):
    """Mirror of ontology_labels.is_unlabeled, kept stdlib-only here on purpose."""
    if name is None:
        return True
    name = str(name).strip()
    if not name or name == str(node_id).strip():
        return True
    return bool(_CURIE_SHAPED.match(name))


def read_unlabeled_curies(nodes_path, category="biolink:Disease"):
    with open(nodes_path, newline="", encoding="utf-8") as fh:
        rows = list(csv.DictReader(fh))
    out = set()
    for r in rows:
        if category and r.get("category") != category:
            continue
        nid = (r.get("id") or "").strip()
        if nid and is_unlabeled(nid, r.get("name")):
            out.add(nid)
    return sorted(out)


def read_cache(path):
    if not os.path.exists(path):
        return {}
    with open(path, newline="", encoding="utf-8") as fh:
        return {r["curie"]: r for r in csv.DictReader(fh, delimiter="\t")}


def write_cache(path, records):
    os.makedirs(os.path.dirname(path) or ".", exist_ok=True)
    with open(path, "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=FIELDS, delimiter="\t")
        w.writeheader()
        for curie in sorted(records):
            row = records[curie]
            w.writerow({k: row.get(k, "") for k in FIELDS})


def _get_json(url, timeout, retries, backoff):
    """GET with retries. Returns parsed JSON, or None on 404 / persistent failure."""
    last = None
    for attempt in range(retries + 1):
        try:
            req = urllib.request.Request(
                url,
                headers={
                    "Accept": "application/json",
                    "User-Agent": "biolink-source-mapping/label-fetch",
                },
            )
            with urllib.request.urlopen(req, timeout=timeout) as resp:
                return json.loads(resp.read().decode("utf-8"))
        except urllib.error.HTTPError as e:
            if e.code == 404:
                return None
            last = e
        except Exception as e:  # network hiccup, timeout, bad JSON
            last = e
        if attempt < retries:
            time.sleep(backoff * (2 ** attempt))
    print(f"    giving up: {last}", file=sys.stderr)
    return None


def _pick_term(terms):
    """Choose the best term for a CURIE.

    OLS4 can return the same obo_id from several ontologies that import it.
    Rank non-obsolete above obsolete first, then the defining ontology above an
    importing one. Obsolete wins nothing here: OLS labels obsolete terms
    "obsolete <name>", so a live label from an importing ontology is the more
    useful display name even though the defining ontology is more authoritative.
    """
    labelled = [t for t in terms if (t.get("label") or "").strip()]
    if not labelled:
        return None

    def rank(t):
        return (
            1 if t.get("is_obsolete") else 0,
            0 if t.get("is_defining_ontology") else 1,
        )

    return sorted(labelled, key=rank)[0]


def resolve(curie, timeout, retries, backoff):
    """Return a cache record for one CURIE. status is one of ok / not_found."""
    url = f"{OLS4_TERMS_URL}?{urllib.parse.urlencode({'obo_id': curie, 'size': 20})}"
    data = _get_json(url, timeout, retries, backoff)
    terms = ((data or {}).get("_embedded") or {}).get("terms") or []
    term = _pick_term(terms)
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


def main():
    ap = argparse.ArgumentParser(
        description="Fetch labels for unnamed disease CURIEs from OLS4"
    )
    ap.add_argument("--nodes", default=DEFAULT_NODES)
    ap.add_argument("--out", default=DEFAULT_OUT)
    ap.add_argument("--refresh", action="store_true",
                    help="re-query CURIEs already present in the cache")
    ap.add_argument("--sleep", type=float, default=0.15,
                    help="seconds between requests (default: 0.15)")
    ap.add_argument("--timeout", type=float, default=30.0)
    ap.add_argument("--retries", type=int, default=3)
    ap.add_argument("--backoff", type=float, default=1.0)
    ap.add_argument("--limit", type=int, default=None,
                    help="stop after N lookups (useful for a first smoke test)")
    args = ap.parse_args()

    curies = read_unlabeled_curies(args.nodes)
    if not curies:
        print("no unnamed disease nodes found; nothing to fetch")
        return

    cache = {} if args.refresh else read_cache(args.out)
    todo = [c for c in curies if c not in cache]
    if args.limit:
        todo = todo[: args.limit]

    print(f"unnamed disease nodes: {len(curies)}")
    print(f"already cached:        {len(curies) - len([c for c in curies if c not in cache])}")
    print(f"to fetch:              {len(todo)}")
    if not todo:
        print("cache is complete")
        return

    ok = 0
    for i, curie in enumerate(todo, 1):
        rec = resolve(curie, args.timeout, args.retries, args.backoff)
        cache[curie] = rec
        if rec["status"] == "ok":
            ok += 1
        marker = rec["label"] or "(not found)"
        print(f"  [{i}/{len(todo)}] {curie:20} {marker}")
        # Save as we go so an interrupted run loses nothing.
        if i % 25 == 0:
            write_cache(args.out, cache)
        time.sleep(args.sleep)

    write_cache(args.out, cache)
    print(f"\nresolved {ok}/{len(todo)}")
    print(f"cache -> {args.out}")
    print("next: python src/ontology_labels.py")


if __name__ == "__main__":
    main()
