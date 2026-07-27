"""
lookups.py

Parse the reference files the pipeline needs into plain dicts.

  HGNC complete set    symbol -> entrez, ensembl -> entrez
  mondo_nodes.tsv      external CURIE -> MONDO id, and MONDO id -> name
  labels.tsv           CURIE -> label, for ids MONDO does not cover

Everything here reads a file and returns a dict. No transformation of the graph,
no network calls. fetch.py downloads these files, pipeline.py uses them.
"""

import re
from collections import defaultdict

import pandas as pd

MONDO_RELEASE = "v2026-06-02"
MONDO_NODES_URL = (
    "https://github.com/monarch-initiative/mondo/releases/download/"
    f"{MONDO_RELEASE}/mondo_nodes.tsv"
)

# Open Targets disease ontologies we try to normalize via MONDO xrefs.
NORMALIZE_PREFIXES = frozenset({"EFO", "DOID", "Orphanet", "OTAR"})

# A name that is just the id, or id-shaped, means the node was never labelled.
_CURIE_SHAPED = re.compile(r"^[A-Za-z][A-Za-z0-9]*[:_]\d+$")

_OBO_CURIE_RE = re.compile(r"/obo/([A-Za-z]+)_(\S+)$")
_EFO_URL_RE = re.compile(r"/efo/(EFO_\d+)", re.I)


# --------------------------------------------------------------------------
# shared helpers
# --------------------------------------------------------------------------

def prefix_of(curie):
    """EFO:0006890 -> EFO"""
    curie = str(curie)
    return curie.split(":", 1)[0] if ":" in curie else ""


def is_unlabeled(node_id, name):
    """True when a node still shows an identifier where a name belongs."""
    if name is None or (isinstance(name, float) and pd.isna(name)):
        return True
    name = str(name).strip()
    if not name or name == str(node_id).strip():
        return True
    return bool(_CURIE_SHAPED.match(name))


def find_unlabeled(nodes_df, category="biolink:Disease"):
    """Sorted CURIEs still named by their id."""
    df = nodes_df if category is None else nodes_df[nodes_df["category"] == category]
    return sorted({
        str(nid).strip()
        for nid, name in zip(df["id"], df["name"])
        if is_unlabeled(nid, name)
    })


# --------------------------------------------------------------------------
# HGNC
# --------------------------------------------------------------------------

def load_hgnc(hgnc_path):
    """Build symbol->entrez and ensembl->entrez from the HGNC complete set.

    Both sources unify to NCBIGene this way: DGIdb symbols and Open Targets
    Ensembl ids resolve to the same Entrez id, so one gene is one node.
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


# --------------------------------------------------------------------------
# MONDO release
# --------------------------------------------------------------------------

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


def _read_mondo(mondo_nodes_path, usecols, include_obsolete):
    df = pd.read_csv(mondo_nodes_path, sep="\t", dtype=str,
                     usecols=usecols, low_memory=False)
    df = df[df["id"].str.startswith("MONDO:", na=False)].copy()
    if not include_obsolete:
        obsolete = df["deprecated"].fillna("").str.strip().str.lower() == "true"
        df = df[~obsolete]
    return df


def load_mondo_xrefs(mondo_nodes_path, include_obsolete=False):
    """External CURIE -> MONDO id. Ambiguous xrefs are excluded.

    Returns (xref_map, ambiguous). An external id pointing at more than one MONDO
    term is not a safe rewrite, so it goes in ambiguous and the id stays put.
    """
    df = _read_mondo(mondo_nodes_path, ["id", "xref", "same_as", "deprecated"],
                     include_obsolete)

    xref_to_mondos = defaultdict(set)
    for _, row in df.iterrows():
        mondo_id = str(row["id"]).strip()
        for curie in _xref_tokens(row):
            xref_to_mondos[curie].add(mondo_id)

    xref_map, ambiguous = {}, {}
    for curie, mondo_ids in xref_to_mondos.items():
        if len(mondo_ids) == 1:
            xref_map[curie] = next(iter(mondo_ids))
        else:
            ambiguous[curie] = sorted(mondo_ids)
    return xref_map, ambiguous


def load_mondo_labels(mondo_nodes_path, include_obsolete=False):
    """MONDO id -> disease name."""
    df = _read_mondo(mondo_nodes_path, ["id", "name", "deprecated"], include_obsolete)
    return {
        str(mid).strip(): str(name).strip()
        for mid, name in zip(df["id"], df["name"])
        if pd.notna(name) and str(name).strip()
    }


# --------------------------------------------------------------------------
# OLS4 label cache
# --------------------------------------------------------------------------

def load_ontology_labels(cache_path):
    """CURIE -> {label, ontology, obsolete} from the cache fetch.py writes.

    Rows with an empty label are dropped. They are ids OLS4 could not resolve,
    and the node should keep its id rather than gain a blank name.
    """
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
