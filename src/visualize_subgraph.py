"""
visualize_subgraph.py

Draw the drug -> gene -> disease neighborhood for one gene from the generated edges,
and save it as a PNG in figures/. This makes the connected structure visible: drugs on
the left, the anchor gene in the middle, diseases on the right.

Disease labels come from output/nodes.csv (MONDO names after normalization/enrichment).
Unmapped ids (where name equals the CURIE) are shown with their ontology prefix.

Run from the project root, for example:
  python src/visualize_subgraph.py --gene NCBIGene:6416
  python src/visualize_subgraph.py --gene NCBIGene:6416 --max-diseases 10 --max-drugs 8

NCBIGene:6416 is MAP2K2 (MEK2), a well-connected gene in the sample.
"""

import argparse
import os
import textwrap

import matplotlib
matplotlib.use("Agg")  # no display needed; write straight to file
import matplotlib.pyplot as plt
import networkx as nx
import pandas as pd

EDGES = "output/edges.csv"
NODES = "output/nodes.csv"
FIG_DIR = "figures"


def short(label, keep_prefix=False):
    """Trim the source prefix so labels read cleanly."""
    if keep_prefix:
        return label
    return label.split(":", 1)[1] if ":" in label else label


def load_node_names():
    nodes = pd.read_csv(NODES, dtype=str)
    return dict(zip(nodes["id"], nodes["name"].fillna(nodes["id"])))


def display_label(node_id, node_names, kind, max_len=32):
    """Prefer human-readable names from nodes.csv; fall back to CURIE form."""
    name = node_names.get(node_id, node_id)
    if kind == "disease":
        if name and name != node_id:
            label = name
        else:
            label = short(node_id, keep_prefix=True)
    else:
        label = name if name else short(node_id)

    if len(label) > max_len:
        return textwrap.shorten(label, width=max_len, placeholder="…")
    return label


def build(gene, max_drugs, max_diseases, gene_name):
    edges = pd.read_csv(EDGES, dtype=str)
    node_names = load_node_names()

    drug_edges = edges[(edges["object"] == gene) &
                       (edges["knowledge_source"] == "DGIdb")].head(max_drugs)
    dis_edges = edges[(edges["subject"] == gene) &
                      (edges["knowledge_source"] == "OpenTargets")].head(max_diseases)

    if drug_edges.empty and dis_edges.empty:
        raise SystemExit(f"No edges found for {gene}. Check the id against output/edges.csv.")

    g = nx.DiGraph()
    center = gene_name or display_label(gene, node_names, "gene")
    g.add_node(gene, kind="gene", label=center)

    pos = {}
    pos[gene] = (0, 0)

    # drugs on the left
    drugs = list(drug_edges["subject"])
    for i, (_, row) in enumerate(drug_edges.iterrows()):
        drug_id = row["subject"]
        g.add_node(drug_id, kind="drug",
                   label=display_label(drug_id, node_names, "drug"))
        g.add_edge(drug_id, gene, label=row["source_relationship"])
        y = (i - (len(drugs) - 1) / 2) * 1.0
        pos[drug_id] = (-3, y)

    # diseases on the right
    diseases = list(dis_edges["object"])
    for i, (_, row) in enumerate(dis_edges.iterrows()):
        disease_id = row["object"]
        g.add_node(disease_id, kind="disease",
                   label=display_label(disease_id, node_names, "disease"))
        g.add_edge(gene, disease_id, label="associated_with")
        y = (i - (len(diseases) - 1) / 2) * 1.0
        pos[disease_id] = (3, y)

    labels = {n: g.nodes[n]["label"] for n in g.nodes}
    color = {"drug": "#5B8FF9", "gene": "#F6BD16", "disease": "#5AD8A6"}
    node_colors = [color[g.nodes[n]["kind"]] for n in g.nodes]

    plt.figure(figsize=(14, max(6, 0.55 * max(len(drugs), len(diseases)) + 2)))
    nx.draw_networkx_nodes(g, pos, node_color=node_colors, node_size=1800, alpha=0.95)
    nx.draw_networkx_edges(g, pos, arrows=True, arrowsize=14, edge_color="#999999",
                           node_size=1800)
    nx.draw_networkx_labels(g, pos, labels=labels, font_size=8)

    plt.title(f"Drug -> Gene -> Disease neighborhood for {center}", fontsize=12)
    plt.axis("off")
    plt.tight_layout()

    os.makedirs(FIG_DIR, exist_ok=True)
    out = os.path.join(FIG_DIR, f"subgraph_{short(gene)}.png")
    plt.savefig(out, dpi=150, bbox_inches="tight")
    plt.close()
    print(f"wrote {out}")
    print(f"drugs shown: {len(drugs)}, diseases shown: {len(diseases)}")


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--gene", default="NCBIGene:6416", help="gene node id from edges.csv")
    ap.add_argument("--gene-name", default="MAP2K2", help="readable name for the center node")
    ap.add_argument("--max-drugs", type=int, default=8)
    ap.add_argument("--max-diseases", type=int, default=10)
    args = ap.parse_args()
    build(args.gene, args.max_drugs, args.max_diseases, args.gene_name)
