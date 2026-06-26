"""
mondo_labels.py

Load MONDO id -> disease name from a MONDO release asset (mondo_nodes.tsv).

Source: https://github.com/monarch-initiative/mondo/releases/tag/v2026-06-02
Asset: mondo_nodes.tsv (columns id, name, deprecated, ...)
"""

import pandas as pd

DEFAULT_MONDO_RELEASE = "v2026-06-02"
DEFAULT_MONDO_NODES_URL = (
    "https://github.com/monarch-initiative/mondo/releases/download/"
    f"{DEFAULT_MONDO_RELEASE}/mondo_nodes.tsv"
)


def load_mondo_labels(mondo_nodes_path, include_obsolete=False):
    """Return {MONDO:0000004: 'adrenocortical insufficiency', ...} from mondo_nodes.tsv."""
    df = pd.read_csv(
        mondo_nodes_path,
        sep="\t",
        dtype=str,
        usecols=["id", "name", "deprecated"],
        low_memory=False,
    )
    df = df[df["id"].str.startswith("MONDO:", na=False)].copy()
    if not include_obsolete:
        obsolete = df["deprecated"].fillna("").str.strip().str.lower() == "true"
        df = df[~obsolete]
    labels = {}
    for mid, name in zip(df["id"], df["name"]):
        if pd.notna(name) and str(name).strip():
            labels[str(mid).strip()] = str(name).strip()
    return labels


def apply_mondo_labels(nodes_df, labels):
    """Replace MONDO disease node names (currently the raw id) with MONDO labels."""
    nodes = nodes_df.copy()
    is_mondo = nodes["id"].astype(str).str.startswith("MONDO:")
    old_names = nodes.loc[is_mondo, "name"].copy()
    mapped = nodes.loc[is_mondo, "id"].map(labels)
    nodes.loc[is_mondo, "name"] = mapped.fillna(nodes.loc[is_mondo, "name"])

    report = pd.DataFrame({
        "mondo_id": nodes.loc[is_mondo, "id"].values,
        "previous_name": old_names.values,
        "mondo_label": mapped.values,
        "matched": mapped.notna().values,
    })
    return nodes, report
