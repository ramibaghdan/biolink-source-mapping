"""
visualize_subgraph.py

Draw the drug -> gene -> disease neighborhood for one gene from the generated edges,
and save it as a PNG in figures/. Drugs on the left, anchor gene in the center, diseases
on the right. Labels come from output/nodes.csv (MONDO names after normalization).

Run from the project root:
  python src/visualize_subgraph.py --gene NCBIGene:6416
"""

import argparse
import os
import textwrap

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch
import pandas as pd

EDGES = "output/edges.csv"
NODES = "output/nodes.csv"
FIG_DIR = "figures"

COLORS = {"drug": "#5B8FF9", "gene": "#F6BD16", "disease": "#5AD8A6"}
FONT_SIZE = 8
WRAP_WIDTH = {"drug": 14, "gene": 12, "disease": 22}
PAD_X = 0.35
PAD_Y = 0.22
V_GAP = 0.35
COL_GAP = 1.2


def short(label, keep_prefix=False):
    if keep_prefix:
        return label
    return label.split(":", 1)[1] if ":" in label else label


def load_node_names():
    nodes = pd.read_csv(NODES, dtype=str)
    return dict(zip(nodes["id"], nodes["name"].fillna(nodes["id"])))


def raw_label(node_id, node_names, kind):
    name = node_names.get(node_id, node_id)
    if kind == "disease":
        if name and name != node_id:
            return name
        return short(node_id, keep_prefix=True)
    return name if name else short(node_id)


def wrap_label(text, kind):
    return textwrap.fill(
        str(text),
        width=WRAP_WIDTH[kind],
        break_long_words=False,
        break_on_hyphens=True,
    )


class LabeledNode:
    """Rounded box drawn to fit its label text."""

    def __init__(self, node_id, kind, label):
        self.node_id = node_id
        self.kind = kind
        self.label = label
        self.x = 0.0
        self.y = 0.0
        self.w = 0.0
        self.h = 0.0
        self.patch = None
        self.text = None

    def measure(self, ax, renderer):
        self.text = ax.text(
            0, 0, self.label, ha="center", va="center",
            fontsize=FONT_SIZE, linespacing=1.2, visible=False,
        )
        bbox = self.text.get_window_extent(renderer).transformed(ax.transData.inverted())
        self.text.remove()
        self.text = None
        self.w = bbox.width + 2 * PAD_X
        self.h = bbox.height + 2 * PAD_Y

    def draw(self, ax):
        rounding = min(self.w, self.h) * 0.18
        self.patch = FancyBboxPatch(
            (self.x - self.w / 2, self.y - self.h / 2),
            self.w, self.h,
            boxstyle=f"round,pad=0,rounding_size={rounding}",
            facecolor=COLORS[self.kind],
            edgecolor="white",
            linewidth=0.8,
            alpha=0.95,
            zorder=2,
        )
        ax.add_patch(self.patch)
        self.text = ax.text(
            self.x, self.y, self.label, ha="center", va="center",
            fontsize=FONT_SIZE, linespacing=1.2, zorder=3,
        )

    @property
    def left(self):
        return self.x - self.w / 2

    @property
    def right(self):
        return self.x + self.w / 2


def stack_column(nodes, x):
    """Vertically stack nodes at column x without overlap."""
    if not nodes:
        return
    total = sum(n.h for n in nodes) + V_GAP * (len(nodes) - 1)
    y = total / 2
    for node in nodes:
        node.x = x
        node.y = y - node.h / 2
        y -= node.h + V_GAP


def column_x(nodes, side, ref_x, ref_w):
    """Place a column beside a reference node."""
    max_w = max(n.w for n in nodes)
    if side == "left":
        return ref_x - ref_w / 2 - COL_GAP - max_w / 2
    return ref_x + ref_w / 2 + COL_GAP + max_w / 2


def draw_arrow(ax, src, dst):
    ax.add_patch(FancyArrowPatch(
        (src.right, src.y), (dst.left, dst.y),
        arrowstyle="-|>", mutation_scale=12,
        color="#888888", linewidth=1.0,
        shrinkA=0, shrinkB=0, zorder=1,
    ))


def build(gene, max_drugs, max_diseases, gene_name):
    edges = pd.read_csv(EDGES, dtype=str)
    node_names = load_node_names()

    drug_edges = edges[(edges["object"] == gene) &
                       (edges["knowledge_source"] == "DGIdb")].head(max_drugs)
    dis_edges = edges[(edges["subject"] == gene) &
                      (edges["knowledge_source"] == "OpenTargets")].head(max_diseases)

    if drug_edges.empty and dis_edges.empty:
        raise SystemExit(f"No edges found for {gene}. Check the id against output/edges.csv.")

    center_text = gene_name or raw_label(gene, node_names, "gene")
    gene_node = LabeledNode(gene, "gene", wrap_label(center_text, "gene"))

    drug_nodes = [
        LabeledNode(row["subject"], "drug", wrap_label(raw_label(row["subject"], node_names, "drug"), "drug"))
        for _, row in drug_edges.iterrows()
    ]
    disease_nodes = [
        LabeledNode(row["object"], "disease", wrap_label(raw_label(row["object"], node_names, "disease"), "disease"))
        for _, row in dis_edges.iterrows()
    ]

    # Measure all label boxes on a scratch axis (same font/size as final figure).
    fig, ax = plt.subplots(figsize=(8, 8))
    ax.set_xlim(-10, 10)
    ax.set_ylim(-10, 10)
    renderer = fig.canvas.get_renderer()
    for node in [gene_node] + drug_nodes + disease_nodes:
        node.measure(ax, renderer)
    plt.close(fig)

    gene_node.x = 0.0
    gene_node.y = 0.0
    stack_column(drug_nodes, column_x(drug_nodes, "left", gene_node.x, gene_node.w))
    stack_column(disease_nodes, column_x(disease_nodes, "right", gene_node.x, gene_node.w))

    # Center all columns vertically on the gene node.
    for column in (drug_nodes, disease_nodes):
        if not column:
            continue
        col_center = (column[0].y + column[-1].y) / 2
        shift = gene_node.y - col_center
        for node in column:
            node.y += shift

    all_nodes = drug_nodes + [gene_node] + disease_nodes
    xs = [n.left for n in all_nodes] + [n.right for n in all_nodes]
    ys = [n.y - n.h / 2 for n in all_nodes] + [n.y + n.h / 2 for n in all_nodes]

    fig_h = max(8, (max(ys) - min(ys)) + 2)
    fig_w = max(14, (max(xs) - min(xs)) + 3)
    fig, ax = plt.subplots(figsize=(fig_w, fig_h))
    ax.set_xlim(min(xs) - 0.8, max(xs) + 0.8)
    ax.set_ylim(min(ys) - 0.6, max(ys) + 0.6)

    for drug in drug_nodes:
        draw_arrow(ax, drug, gene_node)
    for disease in disease_nodes:
        draw_arrow(ax, gene_node, disease)

    for node in all_nodes:
        node.draw(ax)

    ax.set_title(f"Drug -> Gene -> Disease neighborhood for {center_text}", fontsize=12, pad=12)
    ax.axis("off")

    os.makedirs(FIG_DIR, exist_ok=True)
    out = os.path.join(FIG_DIR, f"subgraph_{short(gene)}.png")
    fig.savefig(out, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"wrote {out}")
    print(f"drugs shown: {len(drug_nodes)}, diseases shown: {len(disease_nodes)}")


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--gene", default="NCBIGene:6416")
    ap.add_argument("--gene-name", default="MAP2K2")
    ap.add_argument("--max-drugs", type=int, default=8)
    ap.add_argument("--max-diseases", type=int, default=10)
    args = ap.parse_args()
    build(args.gene, args.max_drugs, args.max_diseases, args.gene_name)
