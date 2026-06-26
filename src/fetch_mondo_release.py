"""
fetch_mondo_release.py

Download mondo_nodes.tsv from a pinned MONDO GitHub release.
No dependencies beyond the standard library.
"""

import argparse
import os
import urllib.request

from mondo_labels import DEFAULT_MONDO_NODES_URL, DEFAULT_MONDO_RELEASE


def fetch(url, dest_path):
    os.makedirs(os.path.dirname(dest_path) or ".", exist_ok=True)
    print(f"downloading {url}")
    print(f"         -> {dest_path}")
    urllib.request.urlretrieve(url, dest_path)
    size_mb = os.path.getsize(dest_path) / (1024 * 1024)
    print(f"done ({size_mb:.1f} MB)")


def main():
    ap = argparse.ArgumentParser(description="Fetch MONDO release assets for label lookup")
    ap.add_argument(
        "--release",
        default=DEFAULT_MONDO_RELEASE,
        help=f"MONDO release tag (default: {DEFAULT_MONDO_RELEASE})",
    )
    ap.add_argument(
        "--out",
        default="data/mondo/mondo_nodes.tsv",
        help="output path for mondo_nodes.tsv",
    )
    args = ap.parse_args()

    url = (
        "https://github.com/monarch-initiative/mondo/releases/download/"
        f"{args.release}/mondo_nodes.tsv"
    )
    fetch(url, args.out)


if __name__ == "__main__":
    main()
