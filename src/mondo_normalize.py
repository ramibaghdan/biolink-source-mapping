"""
mondo_normalize.py

Map external disease CURIEs (EFO, DOID, Orphanet, OTAR) to MONDO using xrefs from
mondo_nodes.tsv (MONDO release v2026-06-02).

Unmapped ids and non-disease ontologies (HP, OBA, GO, MP) are left unchanged.
"""

import re
from collections import defaultdict

import pandas as pd

# Open Targets disease ontologies we attempt to normalize via MONDO xrefs.
DEFAULT_NORMALIZE_PREFIXES = frozenset({"EFO", "DOID", "Orphanet", "OTAR"})

_OBO_CURIE_RE = re.compile(r"/obo/([A-Za-z]+)_(\S+)$")
_EFO_URL_RE = re.compile(r"/efo/(EFO_\d+)", re.I)


def _token_to_curie(token):
    """Parse one xref or same_as token into a CURIE like EFO:0009491."""
    token = str(token).strip()
    if not token:
        return None

    if token.startswith("http"):
        obo = _OBO_CURIE_RE.search(token)
        if obo:
            return f"{obo.group(1)}:{obo.group(2)}"
        efo = _EFO_URL_RE.search(token)
        if efo:
            return efo.group(1).replace("_", ":", 1)
        if "orphanet" in token.lower() or "/ORDO/" in token:
            orphanet_id = token.rstrip("/").rsplit("/", 1)[-1]
            if orphanet_id.isdigit():
                return f"Orphanet:{orphanet_id}"
        return None

    if ":" in token:
        return token
    if "_" in token:
        return token.replace("_", ":", 1)
    return None


def _xref_tokens(row):
    for col in ("xref", "same_as"):
        val = row.get(col)
        if pd.isna(val):
            continue
        for part in str(val).split("|"):
            curie = _token_to_curie(part)
            if curie:
                yield curie


def load_mondo_xrefs(mondo_nodes_path, include_obsolete=False):
    """Build external CURIE -> MONDO id from mondo_nodes.tsv xrefs.

    Ambiguous xrefs (one external id -> multiple MONDO terms) are excluded.
    """
    df = pd.read_csv(
        mondo_nodes_path,
        sep="\t",
        dtype=str,
        usecols=["id", "xref", "same_as", "deprecated"],
        low_memory=False,
    )
    df = df[df["id"].str.startswith("MONDO:", na=False)].copy()
    if not include_obsolete:
        obsolete = df["deprecated"].fillna("").str.strip().str.lower() == "true"
        df = df[~obsolete]

    xref_to_mondos = defaultdict(set)
    for _, row in df.iterrows():
        mondo_id = str(row["id"]).strip()
        for curie in _xref_tokens(row):
            xref_to_mondos[curie].add(mondo_id)

    xref_map = {}
    ambiguous = {}
    for curie, mondo_ids in xref_to_mondos.items():
        if len(mondo_ids) == 1:
            xref_map[curie] = next(iter(mondo_ids))
        else:
            ambiguous[curie] = sorted(mondo_ids)

    return xref_map, ambiguous


def _dedupe_nodes(nodes_df):
    rows = []
    for node_id, group in nodes_df.groupby("id", sort=False):
        category = group.iloc[0]["category"]
        names = [str(n) for n in group["name"].tolist()]
        best = next((n for n in names if n != node_id), names[0])
        rows.append({"id": node_id, "category": category, "name": best})
    return pd.DataFrame(rows)


def apply_mondo_normalization(
    nodes_df,
    edges_df,
    xref_map,
    normalize_prefixes=None,
):
    """Rewrite disease nodes/edge objects to MONDO where xref lookup succeeds."""
    normalize_prefixes = set(normalize_prefixes or DEFAULT_NORMALIZE_PREFIXES)

    disease_nodes = nodes_df[nodes_df["category"] == "biolink:Disease"]
    id_map = {}
    report_rows = []

    for source_id in disease_nodes["id"].astype(str):
        prefix = source_id.split(":", 1)[0]
        if prefix == "MONDO":
            id_map[source_id] = source_id
            report_rows.append({
                "source_id": source_id,
                "mondo_id": source_id,
                "prefix": prefix,
                "normalized": False,
                "status": "already_mondo",
            })
        elif prefix in normalize_prefixes:
            mondo_id = xref_map.get(source_id)
            if mondo_id:
                id_map[source_id] = mondo_id
                report_rows.append({
                    "source_id": source_id,
                    "mondo_id": mondo_id,
                    "prefix": prefix,
                    "normalized": True,
                    "status": "mapped",
                })
            else:
                id_map[source_id] = source_id
                report_rows.append({
                    "source_id": source_id,
                    "mondo_id": "",
                    "prefix": prefix,
                    "normalized": False,
                    "status": "unmapped",
                })
        else:
            id_map[source_id] = source_id
            report_rows.append({
                "source_id": source_id,
                "mondo_id": "",
                "prefix": prefix,
                "normalized": False,
                "status": "skipped_prefix",
            })

    nodes = nodes_df.copy()
    nodes["id"] = nodes["id"].astype(str).map(lambda x: id_map.get(x, x))
    nodes = _dedupe_nodes(nodes)

    edges = edges_df.copy()
    disease_edge = edges["object_category"] == "biolink:Disease"
    edges.loc[disease_edge, "object"] = (
        edges.loc[disease_edge, "object"].astype(str).map(lambda x: id_map.get(x, x))
    )

    report = pd.DataFrame(report_rows)
    return nodes, edges, report
