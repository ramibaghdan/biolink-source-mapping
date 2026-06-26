"""
ingest.py

Load the two sampled source files, validate that the columns we need are present,
clean them, and return tidy dataframes ready for mapping. No network calls.

DGIdb sample: drug-gene interactions (gene_symbol, drug_name, interaction_type, ensembl_id)
Open Targets sample: gene-disease associations (ensembl_id, disease_efo_id, association_score)
"""

import os
import pandas as pd

RAW_DIR = "data/raw"


def _clean_str_cols(df):
    for c in df.columns:
        if df[c].dtype == object:
            df[c] = df[c].astype(str).str.strip()
    return df


def load_dgidb(path=None):
    path = path or os.path.join(RAW_DIR, "dgidb_sample.csv")
    df = pd.read_csv(path, dtype=str)
    df.columns = [c.strip().lower() for c in df.columns]

    needed = {"gene_symbol", "drug_name"}
    missing = needed - set(df.columns)
    if missing:
        raise SystemExit(f"DGIdb sample missing columns: {missing}. Found: {list(df.columns)}")

    df = _clean_str_cols(df)
    # Drop rows that cannot form a triple.
    df = df.dropna(subset=["gene_symbol", "drug_name"])
    df = df[(df["gene_symbol"] != "") & (df["drug_name"] != "")]
    # An interaction_type may be absent; normalize all the missing forms to
    # other/unknown so those rows map to a general predicate instead of being dropped.
    if "interaction_type" not in df.columns:
        df["interaction_type"] = "other/unknown"
    df["interaction_type"] = (df["interaction_type"]
                              .replace({"": "other/unknown", "nan": "other/unknown",
                                        "NaN": "other/unknown", "None": "other/unknown"}))
    df["interaction_type"] = df["interaction_type"].fillna("other/unknown")
    df = df.drop_duplicates()
    return df.reset_index(drop=True)


def load_opentargets(path=None):
    path = path or os.path.join(RAW_DIR, "opentargets_sample.csv")
    df = pd.read_csv(path, dtype=str)
    df.columns = [c.strip().lower() for c in df.columns]

    # The gene column may arrive under several names depending on the sampler/export.
    gene_col = next((c for c in df.columns
                     if c in ("ensembl_id", "targetid", "target_id", "gene_id")), None)
    if gene_col and gene_col != "ensembl_id":
        df = df.rename(columns={gene_col: "ensembl_id"})

    needed = {"ensembl_id", "disease_efo_id"}
    missing = needed - set(df.columns)
    if missing:
        raise SystemExit(
            f"Open Targets sample missing columns: {missing}. Found: {list(df.columns)}. "
            f"The gene/target column may not have been written by the sampler. "
            f"Re-run sample_sources.py, or rename the gene column to 'ensembl_id'."
        )

    df = _clean_str_cols(df)
    df = df.dropna(subset=["ensembl_id", "disease_efo_id"])
    df = df[(df["ensembl_id"] != "") & (df["disease_efo_id"] != "")]
    # One row per gene-disease pair (the sample may carry duplicates by datasource).
    df = df.drop_duplicates(subset=["ensembl_id", "disease_efo_id"])
    return df.reset_index(drop=True)


if __name__ == "__main__":
    d = load_dgidb()
    o = load_opentargets()
    print(f"DGIdb rows: {len(d)}")
    print(f"Open Targets rows: {len(o)}")
    print("\nDGIdb columns:", list(d.columns))
    print("Open Targets columns:", list(o.columns))
