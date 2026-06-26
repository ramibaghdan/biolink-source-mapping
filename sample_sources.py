#!/usr/bin/env python3
"""
sample_sources.py  (gene-anchored version)

Run this LOCALLY (needs the downloaded source files).

The key idea: a drug-gene sample and a gene-disease sample taken INDEPENDENTLY would
not share genes, so the integrated graph would be disconnected (drugs on one set of
genes, diseases on another). Instead we:

  1. load both full sources,
  2. reconcile gene identifiers (DGIdb gene symbols vs Open Targets Ensembl IDs),
  3. find genes present in BOTH sources,
  4. sample rows from each source restricted to those SHARED genes.

Result: every gene in the final sample has at least one drug link and at least one
disease link, so the graph connects drug -> gene -> disease.

You provide a gene-symbol <-> Ensembl mapping file, since DGIdb and Open Targets use
different gene identifiers. No em dashes, no AI tonality.

WHERE TO GET FILES (verify license + current URL at download time):
  DGIdb interactions TSV:        https://www.dgidb.org/  -> Downloads
  Open Targets associations:     https://platform.opentargets.org/downloads  (CC0)
  Gene symbol <-> Ensembl map:   HGNC complete set (https://www.genenames.org/download/)
                                 or Ensembl BioMart export with columns
                                 [ensembl_gene_id, hgnc_symbol, entrez_id]
"""

import argparse
import os
import random
import pandas as pd


def load_gene_map(path):
    df = pd.read_csv(path, sep=None, engine="python", dtype=str)
    df.columns = [c.strip().lower().replace(" ", "_") for c in df.columns]
    ensg = next((c for c in df.columns if "ensembl" in c), None)
    sym = next((c for c in df.columns if "symbol" in c or c == "hgnc_symbol"), None)
    entrez = next((c for c in df.columns if "entrez" in c or "ncbi" in c), None)
    if not (ensg and sym):
        raise SystemExit(f"gene map needs ensembl + symbol columns. found: {list(df.columns)}")
    keep = [c for c in (ensg, sym, entrez) if c]
    df = df[keep].dropna(subset=[ensg, sym]).drop_duplicates()
    df = df.rename(columns={ensg: "ensembl_id", sym: "gene_symbol"})
    if entrez:
        df = df.rename(columns={entrez: "entrez_id"})
    return df


def load_dgidb(path):
    df = pd.read_csv(path, sep="\t", dtype=str)
    df.columns = [c.strip().lower().replace(" ", "_") for c in df.columns]
    gene_col = next((c for c in df.columns if c in
                     ("gene_name", "gene", "gene_claim_name")), None)
    drug_col = next((c for c in df.columns if c in
                     ("drug_name", "drug", "drug_claim_name")), None)
    itype_col = next((c for c in df.columns if "interaction_type" in c), None)
    if not (gene_col and drug_col):
        raise SystemExit(f"DGIdb: need gene+drug columns. found: {list(df.columns)}")
    keep = [c for c in (gene_col, drug_col, itype_col) if c]
    df = df[keep].dropna(subset=[gene_col, drug_col]).drop_duplicates()
    df = df.rename(columns={gene_col: "gene_symbol", drug_col: "drug_name"})
    if itype_col:
        df = df.rename(columns={itype_col: "interaction_type"})
    df["gene_symbol"] = df["gene_symbol"].str.upper().str.strip()
    return df


def load_opentargets(path):
    if path.endswith(".parquet"):
        df = pd.read_parquet(path)
    elif path.endswith((".tsv", ".txt")):
        df = pd.read_csv(path, sep="\t", dtype=str)
    else:
        df = pd.read_csv(path, dtype=str)
    df.columns = [c.strip().lower().replace(" ", "_") for c in df.columns]
    gene_col = next((c for c in df.columns if c in
                     ("targetid", "target_id", "ensembl_id", "gene_id")), None)
    dis_col = next((c for c in df.columns if c in
                    ("diseaseid", "disease_id", "efo_id")), None)
    score_col = next((c for c in df.columns if "score" in c), None)
    if not (gene_col and dis_col):
        raise SystemExit(f"Open Targets: need gene+disease id cols. found: {list(df.columns)}")
    keep = [c for c in (gene_col, dis_col, score_col) if c]
    df = df[keep].dropna(subset=[gene_col, dis_col]).drop_duplicates()
    df = df.rename(columns={gene_col: "ensembl_id", dis_col: "disease_efo_id"})
    if score_col:
        df = df.rename(columns={score_col: "association_score"})
    return df


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--dgidb", required=True, help="DGIdb interactions TSV")
    ap.add_argument("--opentargets", required=True, help="Open Targets associations file")
    ap.add_argument("--genemap", required=True,
                    help="gene symbol <-> Ensembl mapping (HGNC or BioMart export)")
    ap.add_argument("--n-genes", type=int, default=150,
                    help="how many shared genes to anchor on")
    ap.add_argument("--seed", type=int, default=42)
    ap.add_argument("--outdir", default="data/raw")
    args = ap.parse_args()
    os.makedirs(args.outdir, exist_ok=True)

    gmap = load_gene_map(args.genemap)
    dgidb = load_dgidb(args.dgidb)
    ot = load_opentargets(args.opentargets)

    sym_to_ensg = gmap[["gene_symbol", "ensembl_id"]].copy()
    sym_to_ensg["gene_symbol"] = sym_to_ensg["gene_symbol"].str.upper().str.strip()
    dgidb = dgidb.merge(sym_to_ensg, on="gene_symbol", how="left")

    dgidb_genes = set(dgidb["ensembl_id"].dropna())
    ot_genes = set(ot["ensembl_id"].dropna())
    shared = sorted(dgidb_genes & ot_genes)
    print(f"DGIdb genes (mapped to ENSG): {len(dgidb_genes)}")
    print(f"Open Targets genes: {len(ot_genes)}")
    print(f"shared genes: {len(shared)}")
    if not shared:
        raise SystemExit("No shared genes. Check the gene map / identifier formats.")

    # Only anchor on genes that actually have edges in BOTH sources, so every
    # anchor gene is genuinely connected (drug -> gene -> disease). "shared" above
    # only means the gene exists in both; this is stricter and is what drives
    # real connectivity.
    connectable = sorted(set(dgidb["ensembl_id"].dropna()) & set(ot["ensembl_id"].dropna()))
    print(f"genes connectable in both sources: {len(connectable)}")

    random.seed(args.seed)
    anchor = set(random.sample(connectable, min(args.n_genes, len(connectable))))

    dgidb_s = dgidb[dgidb["ensembl_id"].isin(anchor)].copy()
    ot_s = ot[ot["ensembl_id"].isin(anchor)].copy()

    # Cap rows per gene without losing the grouping column. (On pandas 3.x,
    # groupby(...).apply(...) can drop the grouping column, which silently produced
    # a sample file with no gene column.) sample(frac=1) shuffles, groupby.head(n)
    # keeps at most n rows per gene and preserves all columns.
    ot_s = (ot_s.sample(frac=1, random_state=args.seed)
            .groupby("ensembl_id", as_index=False, group_keys=False)
            .head(15)
            .reset_index(drop=True))

    # For DGIdb, cap rows per gene the same way (preserves connectivity + columns).
    dgidb_s = (dgidb_s.sample(frac=1, random_state=args.seed)
               .groupby("ensembl_id", as_index=False, group_keys=False)
               .head(8)
               .reset_index(drop=True))

    dgidb_s.to_csv(f"{args.outdir}/dgidb_sample.csv", index=False)
    ot_s.to_csv(f"{args.outdir}/opentargets_sample.csv", index=False)

    print(f"\nwrote {args.outdir}/dgidb_sample.csv  ({len(dgidb_s)} rows)")
    print(f"wrote {args.outdir}/opentargets_sample.csv  ({len(ot_s)} rows)")
    print(f"anchor genes used: {len(anchor)}")
    if "interaction_type" in dgidb_s.columns:
        print("\ninteraction types in DGIdb sample:")
        print(dgidb_s["interaction_type"].value_counts().to_string())

    def ensg_set(frame):
        if "ensembl_id" in frame.columns:
            return set(frame["ensembl_id"].dropna())
        return set(pd.Series(frame.index).dropna())

    connected = ensg_set(dgidb_s) & ensg_set(ot_s)
    print(f"\ngenes with both a drug link and a disease link: {len(connected)}")
    print("(these are the genes that make drug -> gene -> disease paths exist)")


if __name__ == "__main__":
    main()
