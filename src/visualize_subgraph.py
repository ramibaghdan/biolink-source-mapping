"""
visualize_subgraph.py

Draw the drug -> gene -> disease neighborhood for one gene from the generated edges,
and save it as a PNG in figures/. This makes the connected structure visible: drugs on
the left, the anchor gene in the middle, diseases on the right.

Run from the project root, for example:
  python src/visualize_subgraph.py --gene NCBIGene:6416
  python src/visualize_subgraph.py --gene NCBIGene:6416 --max-diseases 10 --max-drugs 8

NCBIGene:6416 is MAP2K2 (MEK2), a well-connected gene in the sample.
"""

import argparse
import os
import pandas as pd
import networkx as nx
import matplotlib
matplotlib.use("Agg")  # no display needed; write straight to file
import matplotlib.pyplot as plt

EDGES = "output/edges.csv"
FIG_DIR = "figures"


def short(label, keep_prefix=False):
    """Trim the source prefix so labels read cleanly.

    DRUG_NAME:TRAMETINIB -> TRAMETINIB. For diseases the ontology prefix is meaningful
    (EFO vs MONDO), so keep_prefix=True leaves EFO:0020865 intact.
    """
    if keep_prefix:
        return label
    return label.split(":", 1)[1] if ":" in label else label


def build(gene, max_drugs, max_diseases, gene_name):
    edges = pd.read_csv(EDGES, dtype=str)

    drug_edges = edges[(edges["object"] == gene) &
                       (edges["knowledge_source"] == "DGIdb")].head(max_drugs)
    dis_edges = edges[(edges["subject"] == gene) &
                      (edges["knowledge_source"] == "OpenTargets")].head(max_diseases)

    if drug_edges.empty and dis_edges.empty:
        raise SystemExit(f"No edges found for {gene}. Check the id against output/edges.csv.")

    g = nx.DiGraph()
    center = gene_name or short(gene)
    g.add_node(center, kind="gene")

    pos = {}
    pos[center] = (0, 0)

    # drugs on the left
    drugs = list(drug_edges["subject"])
    for i, (_, row) in enumerate(drug_edges.iterrows()):
        d = short(row["subject"])
        g.add_node(d, kind="drug")
        g.add_edge(d, center, label=row["source_relationship"])
        y = (i - (len(drugs) - 1) / 2) * 1.0
        pos[d] = (-3, y)

    # diseases on the right
    diseases = list(dis_edges["object"])
    for i, (_, row) in enumerate(dis_edges.iterrows()):
        ds = short(row["object"], keep_prefix=True)
        g.add_node(ds, kind="disease")
        g.add_edge(center, ds, label="associated_with")
        y = (i - (len(diseases) - 1) / 2) * 1.0
        pos[ds] = (3, y)

    # colors by kind
    color = {"drug": "#5B8FF9", "gene": "#F6BD16", "disease": "#5AD8A6"}
    node_colors = [color[g.nodes[n]["kind"]] for n in g.nodes]

    plt.figure(figsize=(12, max(6, 0.5 * max(len(drugs), len(diseases)) + 2)))
    nx.draw_networkx_nodes(g, pos, node_color=node_colors, node_size=1600, alpha=0.95)
    nx.draw_networkx_edges(g, pos, arrows=True, arrowsize=14, edge_color="#999999",
                           node_size=1600)
    nx.draw_networkx_labels(g, pos, font_size=8)

    plt.title(f"Drug -> Gene -> Disease neighborhood for {center}", fontsize=12)
    plt.axis("off")
    plt.tight_layout()

    os.makedirs(FIG_DIR, exist_ok=True)
    out = os.path.join(FIG_DIR, f"subgraph_{short(gene)}.png")
    plt.savefig(out, dpi=150, bbox_inches="tight")
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
