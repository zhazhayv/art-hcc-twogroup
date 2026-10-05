# -*- coding: utf-8 -*-
"""Figure 1: locked two-level intersections."""
from __future__ import annotations

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Circle

from constants import FIGURES, FULL67, NETWORK6, STRICT14


def venn(ax, left, right, title):
    a, b = set(left), set(right)
    ax.set_xlim(-1.6, 1.6)
    ax.set_ylim(-1.3, 1.3)
    ax.set_aspect("equal")
    ax.add_patch(Circle((-0.45, 0), 0.9, fill=True, alpha=0.35, color="#4C72B0"))
    ax.add_patch(Circle((0.45, 0), 0.9, fill=True, alpha=0.35, color="#C44E52"))
    ax.text(-0.95, 0.95, "Artesunate\nfiltered", ha="center", fontsize=8)
    ax.text(0.95, 0.95, "GeneCards\nHCC >2", ha="center", fontsize=8)
    ax.text(-0.7, 0, str(len(a - b)), ha="center", va="center", fontsize=11)
    ax.text(0.7, 0, str(len(b - a)), ha="center", va="center", fontsize=11)
    ax.text(0, 0, str(len(a & b)), ha="center", va="center", fontsize=12, fontweight="bold")
    ax.axis("off")
    ax.set_title(title, fontsize=10)


def main():
    fig, axes = plt.subplots(1, 2, figsize=(10, 4.2), dpi=300)
    # Left: Strict_14 as artesunate-filtered vs implied disease (intersection is 14)
    venn(axes[0], STRICT14, STRICT14, "A  Strict_14 (primary)")
    axes[0].texts  # placeholder
    # Custom labels for A: we only have intersection list; draw counts from manuscript: 14 shared
    axes[0].cla()
    venn_strict(axes[0])
    venn_full(axes[1])
    fig.tight_layout()
    fig.savefig(FIGURES / "fig1_two_level_intersection.png", bbox_inches="tight")
    print("fig1", FIGURES / "fig1_two_level_intersection.png")


def venn_strict(ax):
    ax.set_xlim(-1.7, 1.7)
    ax.set_ylim(-1.35, 1.35)
    ax.set_aspect("equal")
    ax.add_patch(Circle((-0.45, 0), 0.95, fill=True, alpha=0.35, color="#4C72B0"))
    ax.add_patch(Circle((0.45, 0), 0.95, fill=True, alpha=0.35, color="#C44E52"))
    ax.text(-1.05, 1.05, "CTD + Swiss P>=0.1", ha="center", fontsize=8)
    ax.text(1.05, 1.05, "GeneCards >2", ha="center", fontsize=8)
    ax.text(-0.75, 0, "1", ha="center", fontsize=11)  # PTGES2 only on compound side among high-p
    ax.text(0.85, 0, "3741", ha="center", fontsize=9)
    ax.text(0.0, 0.0, "14", ha="center", fontsize=14, fontweight="bold")
    ax.text(0, -1.22, "Network cassette inside Strict_14: " + ",".join(NETWORK6), fontsize=7, ha="center")
    ax.axis("off")
    ax.set_title("A  Strict_14", fontsize=11)


def venn_full(ax):
    ax.set_xlim(-1.7, 1.7)
    ax.set_ylim(-1.35, 1.35)
    ax.set_aspect("equal")
    ax.add_patch(Circle((-0.45, 0), 0.95, fill=True, alpha=0.35, color="#4C72B0"))
    ax.add_patch(Circle((0.45, 0), 0.95, fill=True, alpha=0.35, color="#C44E52"))
    ax.text(-1.05, 1.05, "Swiss all + CTD", ha="center", fontsize=8)
    ax.text(1.05, 1.05, "GeneCards >2", ha="center", fontsize=8)
    ax.text(-0.75, 0, "41", ha="center", fontsize=11)
    ax.text(0.85, 0, "3688", ha="center", fontsize=10)
    ax.text(0.0, 0.0, "67", ha="center", fontsize=14, fontweight="bold")
    ax.text(0, -1.22, f"Full_67 n={len(FULL67)}; STRING 660 edges only", fontsize=7, ha="center")
    ax.axis("off")
    ax.set_title("B  Full_67 (exploratory)", fontsize=11)


if __name__ == "__main__":
    main()
