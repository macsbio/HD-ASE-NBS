# -*- coding: utf-8 -*-

import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


def plot_cluster_stability(cluster_stab_df, out_fn=None, title=None):
    """
    Barplot of within-cluster vs between-cluster consensus.

    Input:
        cluster_stab_df = output of cluster_stability_table()

    Columns required:
        cluster, within_mean, between_mean
    """
    df = cluster_stab_df.copy()
    df = df.sort_values("cluster").reset_index(drop=True)

    x = np.arange(len(df))
    width = 0.35

    fig, ax = plt.subplots(figsize=(7, 5))

    ax.bar(x - width / 2.0, df["within_mean"].values, width, label="Within-cluster")
    ax.bar(x + width / 2.0, df["between_mean"].values, width, label="Between-cluster")

    ax.set_xticks(x)
    ax.set_xticklabels([str(c) for c in df["cluster"].values])
    ax.set_xlabel("Cluster")
    ax.set_ylabel("Mean consensus")
    ax.set_ylim(0, 1)

    if title is None:
        title = "Cluster stability summary"
    ax.set_title(title)

    ax.legend(frameon=False)
    plt.tight_layout()

    if out_fn is not None:
        plt.savefig(out_fn, dpi=300, bbox_inches="tight")
        plt.close()
    else:
        plt.show()


def plot_cluster_separation(cluster_stab_df, out_fn=None, title=None):
    """
    Simple barplot of separation = within_mean - between_mean.
    """
    df = cluster_stab_df.copy()
    df = df.sort_values("cluster").reset_index(drop=True)

    x = np.arange(len(df))

    fig, ax = plt.subplots(figsize=(6, 4.5))
    ax.bar(x, df["separation"].values)

    ax.set_xticks(x)
    ax.set_xticklabels([str(c) for c in df["cluster"].values])
    ax.set_xlabel("Cluster")
    ax.set_ylabel("Separation")
    ax.set_ylim(bottom=0)

    if title is None:
        title = "Cluster separation"
    ax.set_title(title)

    plt.tight_layout()

    if out_fn is not None:
        plt.savefig(out_fn, dpi=300, bbox_inches="tight")
        plt.close()
    else:
        plt.show()


def plot_sample_confidence(sample_conf_df, out_fn=None, title=None, sort_within_cluster=True):
    """
    Dotplot of sample assignment margin by cluster.

    Input:
        sample_conf_df = output of sample_assignment_confidence()

    Columns required:
        cluster, margin
    """
    df = sample_conf_df.copy()
    df["sample_id"] = df.index.astype(str)

    if sort_within_cluster:
        df = df.sort_values(["cluster", "margin"], ascending=[True, False])
    else:
        df = df.sort_values(["cluster", "sample_id"])

    fig, ax = plt.subplots(figsize=(10, 4.5))

    xpos = np.arange(len(df))
    ax.scatter(xpos, df["margin"].values)

    # vertical separators between clusters
    cluster_changes = []
    prev = None
    for i, cl in enumerate(df["cluster"].values):
        if prev is None:
            prev = cl
        elif cl != prev:
            cluster_changes.append(i - 0.5)
            prev = cl

    for xc in cluster_changes:
        ax.axvline(x=xc, linestyle="--", linewidth=0.8)

    ax.axhline(y=0, linestyle="--", linewidth=0.8)

    ax.set_xlabel("Samples")
    ax.set_ylabel("Assignment margin")
    if title is None:
        title = "Sample-level assignment confidence"
    ax.set_title(title)

    # label cluster centers on x-axis
    centers = []
    labels = []
    for cl in sorted(df["cluster"].dropna().unique()):
        idx = np.where(df["cluster"].values == cl)[0]
        centers.append(idx.mean())
        labels.append(str(cl))

    ax.set_xticks(centers)
    ax.set_xticklabels(labels)

    plt.tight_layout()

    if out_fn is not None:
        plt.savefig(out_fn, dpi=300, bbox_inches="tight")
        plt.close()
    else:
        plt.show()


def plot_sample_confidence_boxplot(sample_conf_df, out_fn=None, title=None):
    """
    Boxplot of assignment margin per cluster.
    """
    df = sample_conf_df.copy()
    clusters = sorted(df["cluster"].dropna().unique())
    data = [df.loc[df["cluster"] == cl, "margin"].dropna().values for cl in clusters]

    fig, ax = plt.subplots(figsize=(6.5, 4.5))
    ax.boxplot(data, labels=[str(c) for c in clusters])

    ax.axhline(y=0, linestyle="--", linewidth=0.8)
    ax.set_xlabel("Cluster")
    ax.set_ylabel("Assignment margin")

    if title is None:
        title = "Assignment margin by cluster"
    ax.set_title(title)

    plt.tight_layout()

    if out_fn is not None:
        plt.savefig(out_fn, dpi=300, bbox_inches="tight")
        plt.close()
    else:
        plt.show()


def plot_silhouette_from_samples(sil_sample_df, out_fn=None, title=None):
    """
    Standard silhouette plot from silhouette_sample output.

    Input:
        sil_sample_df = sample-level output from consensus_silhouette()

    Columns required:
        cluster, silhouette
    """
    df = sil_sample_df.copy()
    df = df.sort_values(["cluster", "silhouette"], ascending=[True, False])

    clusters = sorted(df["cluster"].dropna().unique())

    fig, ax = plt.subplots(figsize=(7, 5))

    y_lower = 10
    yticks = []

    for cl in clusters:
        vals = df.loc[df["cluster"] == cl, "silhouette"].dropna().values
        vals.sort()
        vals = vals[::-1]  # descending
        size_cl = len(vals)
        y_upper = y_lower + size_cl

        ax.fill_betweenx(
            np.arange(y_lower, y_upper),
            0,
            vals
        )

        yticks.append(y_lower + 0.5 * size_cl)
        y_lower = y_upper + 10

    ax.axvline(x=np.mean(df["silhouette"].dropna().values), linestyle="--", linewidth=1.0)

    ax.set_xlabel("Silhouette value")
    ax.set_ylabel("Cluster")
    ax.set_yticks(yticks)
    ax.set_yticklabels([str(c) for c in clusters])
    ax.set_xlim([-1, 1])

    if title is None:
        title = "Silhouette plot"
    ax.set_title(title)

    plt.tight_layout()

    if out_fn is not None:
        plt.savefig(out_fn, dpi=300, bbox_inches="tight")
        plt.close()
    else:
        plt.show()


def plot_pairwise_cluster_consensus(pairwise_df, out_fn=None, title=None):
    """
    Barplot of mean consensus between cluster pairs.
    Lower means better separation.
    """
    df = pairwise_df.copy()
    df["pair"] = df["cluster_1"].astype(str) + " vs " + df["cluster_2"].astype(str)
    df = df.sort_values("pair").reset_index(drop=True)

    x = np.arange(len(df))

    fig, ax = plt.subplots(figsize=(6.5, 4.5))
    ax.bar(x, df["mean_consensus"].values)

    ax.set_xticks(x)
    ax.set_xticklabels(df["pair"].values, rotation=45, ha="right")
    ax.set_ylabel("Mean between-cluster consensus")
    ax.set_xlabel("Cluster pair")
    ax.set_ylim(0, 1)

    if title is None:
        title = "Pairwise cluster similarity"
    ax.set_title(title)

    plt.tight_layout()

    if out_fn is not None:
        plt.savefig(out_fn, dpi=300, bbox_inches="tight")
        plt.close()
    else:
        plt.show()


def save_all_stability_plots(report, outdir, prefix="pyNBS"):
    """
    Convenience wrapper for all major plots.

    report = output of full_stability_report()
    """
    if not os.path.exists(outdir):
        os.makedirs(outdir)

    plot_cluster_stability(
        report["cluster_stability"],
        out_fn=os.path.join(outdir, "%s_cluster_stability_barplot.png" % prefix)
    )

    plot_cluster_separation(
        report["cluster_stability"],
        out_fn=os.path.join(outdir, "%s_cluster_separation_barplot.png" % prefix)
    )

    plot_sample_confidence(
        report["sample_confidence"],
        out_fn=os.path.join(outdir, "%s_sample_confidence_dotplot.png" % prefix)
    )

    plot_sample_confidence_boxplot(
        report["sample_confidence"],
        out_fn=os.path.join(outdir, "%s_sample_confidence_boxplot.png" % prefix)
    )

    plot_pairwise_cluster_consensus(
        report["pairwise_cluster_consensus"],
        out_fn=os.path.join(outdir, "%s_pairwise_cluster_consensus.png" % prefix)
    )

    if "silhouette_sample" in report:
        plot_silhouette_from_samples(
            report["silhouette_sample"],
            out_fn=os.path.join(outdir, "%s_silhouette_plot.png" % prefix)
        )