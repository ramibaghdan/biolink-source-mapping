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

WRAP_WIDTH = {"drug": 11, "gene": 10, "disease": 16}
FONT_SIZE = 7


def short(label, keep_prefix=False):
    """Trim the source prefix so labels read cleanly."""
    if keep_prefix:
        return label
    return label.split(":", 1)[1] if ":" in label else label


def load_node_names():
    nodes = pd.read_csv(NODES, dtype=str)
    return dict(zip(nodes["id"], nodes["name"].fillna(nodes["id"])))


def raw_label(node_id, node_names, kind):
    """Full display text before wrapping."""
    name = node_names.get(node_id, node_id)
    if kind == "disease":
        if name and name != node_id:
            return name
        return short(node_id, keep_prefix=True)
    return name if name else short(node_id)


def wrap_label(text, kind):
    """Wrap label to multiple lines so it fits inside the node circle."""
    return textwrap.fill(
        str(text),
        width=WRAP_WIDTH[kind],
        break_long_words=False,
        break_on_hyphens=True,
    )


def node_size(label, kind):
    """Scale circle area from wrapped line count and longest line."""
    lines = label.split("\n")
    max_len = max(len(line) for line in lines)
    n_lines = len(lines)
    base = {"drug": 1500, "gene": 1700, "disease": 1700}[kind]
    return base + (n_lines - 1) * 450 + max(0, max_len - 8) * 90


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
    center_text = gene_name or raw_label(gene, node_names, "gene")
    center_label = wrap_label(center_text, "gene")
    g.add_node(gene, kind="gene", label=center_label)

    pos = {}
    pos[gene] = (0, 0)

    row_count = max(len(drug_edges), len(dis_edges), 1)
    y_step = max(1.15, 5.5 / row_count)

    drugs = list(drug_edges["subject"])
    for i, (_, row) in enumerate(drug_edges.iterrows()):
        drug_id = row["subject"]
        label = wrap_label(raw_label(drug_id, node_names, "drug"), "drug")
        g.add_node(drug_id, kind="drug", label=label)
        g.add_edge(drug_id, gene, label=row["source_relationship"])
        y = (i - (len(drugs) - 1) / 2) * y_step
        pos[drug_id] = (-4.2, y)

    diseases = list(dis_edges["object"])
    for i, (_, row) in enumerate(dis_edges.iterrows()):
        disease_id = row["object"]
        label = wrap_label(raw_label(disease_id, node_names, "disease"), "disease")
        g.add_node(disease_id, kind="disease", label=label)
        g.add_edge(gene, disease_id, label="associated_with")
        y = (i - (len(diseases) - 1) / 2) * y_step
        pos[disease_id] = (4.2, y)

    node_list = list(g.nodes())
    labels = {n: g.nodes[n]["label"] for n in node_list}
    sizes = [node_size(labels[n], g.nodes[n]["kind"]) for n in node_list]
    colors = [{"drug": "#5B8FF9", "gene": "#F6BD16", "disease": "#5AD8A6"}[g.nodes[n]["kind"]]
              for n in node_list]

    fig_h = max(7, row_count * 0.85 + 2)
    plt.figure(figsize=(15, fig_h))
    nx.draw_networkx_nodes(
        g, pos, nodelist=node_list, node_color=colors, node_size=sizes, alpha=0.95,
    )
    nx.draw_networkx_edges(
        g, pos, arrows=True, arrowsize=12, edge_color="#999999",
        node_size=sizes, min_source_margin=15, min_target_margin=20,
    )
    nx.draw_networkx_labels(
        g, pos, labels=labels, font_size=FONT_SIZE, font_weight="regular",
    )

    plt.title(f"Drug -> Gene -> Disease neighborhood for {center_text}", fontsize=12)
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
