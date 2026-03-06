# pyNBS_stability.py
# Compatibility-oriented version for older pyNBS environments
# No f-strings, no modern-only pandas options

import os
import numpy as np
import pandas as pd
from sklearn.metrics import silhouette_samples, silhouette_score


def _check_inputs(cc_table, cluster_assign):
    if not isinstance(cc_table, pd.DataFrame):
        raise TypeError("cc_table must be a pandas DataFrame")
    if not isinstance(cluster_assign, pd.Series):
        raise TypeError("cluster_assign must be a pandas Series")

    # Keep only overlapping samples and align order
    common = cc_table.index.intersection(cc_table.columns).intersection(cluster_assign.index)
    if len(common) == 0:
        raise ValueError("No overlapping sample IDs between cc_table and cluster_assign")

    cc_table = cc_table.loc[common, common].copy()
    cluster_assign = cluster_assign.loc[common].copy()

    # Basic checks
    if not np.allclose(cc_table.values, cc_table.values.T, atol=1e-8):
        raise ValueError("cc_table must be symmetric")
    if (cc_table.values < 0).any() or (cc_table.values > 1).any():
        raise ValueError("cc_table values should be in [0, 1]")

    return cc_table, cluster_assign


def _distance_from_consensus(cc_table):
    dist = 1.0 - cc_table.astype(float)
    np.fill_diagonal(dist.values, 0.0)
    return dist


def _upper_triangle_values(mat):
    arr = np.asarray(mat)
    iu = np.triu_indices_from(arr, k=1)
    return arr[iu]


def summarize_consensus_matrix(cc_table):
    """
    Global summary of the consensus matrix.
    """
    vals = _upper_triangle_values(cc_table.values)
    return pd.Series({
        "n_samples": cc_table.shape[0],
        "mean_consensus": np.mean(vals),
        "median_consensus": np.median(vals),
        "sd_consensus": np.std(vals, ddof=1),
        "min_consensus": np.min(vals),
        "max_consensus": np.max(vals)
    }, name="global_consensus")


def cluster_stability_table(cc_table, cluster_assign):
    """
    Per-cluster within/between consensus and separation.
    Higher within-cluster consensus and lower between-cluster consensus are better.
    """
    cc_table, cluster_assign = _check_inputs(cc_table, cluster_assign)
    clusters = sorted(cluster_assign.dropna().unique())

    rows = []
    for cl in clusters:
        in_idx = cluster_assign[cluster_assign == cl].index
        out_idx = cluster_assign[cluster_assign != cl].index

        n_in = len(in_idx)
        if n_in < 2:
            rows.append({
                "cluster": cl,
                "n": n_in,
                "within_mean": np.nan,
                "within_median": np.nan,
                "between_mean": np.nan,
                "between_median": np.nan,
                "separation": np.nan
            })
            continue

        within = cc_table.loc[in_idx, in_idx].values
        within_vals = _upper_triangle_values(within)

        if len(out_idx) > 0:
            between_vals = cc_table.loc[in_idx, out_idx].values.flatten()
            between_mean = float(np.mean(between_vals))
            between_median = float(np.median(between_vals))
        else:
            between_mean = np.nan
            between_median = np.nan

        within_mean = float(np.mean(within_vals))
        within_median = float(np.median(within_vals))

        if pd.notnull(between_mean):
            separation = within_mean - between_mean
        else:
            separation = np.nan

        rows.append({
            "cluster": cl,
            "n": n_in,
            "within_mean": within_mean,
            "within_median": within_median,
            "between_mean": between_mean,
            "between_median": between_median,
            "separation": separation
        })

    return pd.DataFrame(rows).sort_values("cluster")


def sample_assignment_confidence(cc_table, cluster_assign):
    """
    Sample-level confidence:
    - mean consensus to own cluster
    - best competing cluster consensus
    - margin = own - best_other
    """
    cc_table, cluster_assign = _check_inputs(cc_table, cluster_assign)
    clusters = sorted(cluster_assign.dropna().unique())

    rows = []
    for sample in cc_table.index:
        own_cluster = cluster_assign.loc[sample]
        own_members = cluster_assign[cluster_assign == own_cluster].index.difference([sample])

        if len(own_members) > 0:
            own_mean = float(cc_table.loc[sample, own_members].mean())
        else:
            own_mean = np.nan

        other_means = {}
        for cl in clusters:
            if cl == own_cluster:
                continue
            members = cluster_assign[cluster_assign == cl].index
            if len(members) == 0:
                continue
            other_means[cl] = float(cc_table.loc[sample, members].mean())

        if len(other_means) > 0:
            best_other_cluster = max(other_means, key=other_means.get)
            best_other_mean = other_means[best_other_cluster]
        else:
            best_other_cluster = np.nan
            best_other_mean = np.nan

        if pd.notnull(own_mean) and pd.notnull(best_other_mean):
            margin = own_mean - best_other_mean
        else:
            margin = np.nan

        rows.append({
            "sample": sample,
            "cluster": own_cluster,
            "own_cluster_mean_consensus": own_mean,
            "best_other_cluster": best_other_cluster,
            "best_other_mean_consensus": best_other_mean,
            "margin": margin
        })

    out = pd.DataFrame(rows).set_index("sample")
    return out.sort_values(["cluster", "margin"], ascending=[True, False])


def consensus_silhouette(cc_table, cluster_assign):
    """
    Silhouette score using precomputed distance = 1 - consensus.
    """
    cc_table, cluster_assign = _check_inputs(cc_table, cluster_assign)
    dist = _distance_from_consensus(cc_table)

    labels = cluster_assign.loc[dist.index].values
    if len(np.unique(labels)) < 2:
        raise ValueError("Need at least 2 clusters for silhouette analysis")

    sil_samples = silhouette_samples(dist.values, labels, metric="precomputed")
    sil_mean = silhouette_score(dist.values, labels, metric="precomputed")

    sample_df = pd.DataFrame({
        "cluster": labels,
        "silhouette": sil_samples
    }, index=dist.index)

    cluster_df = (
        sample_df.groupby("cluster")["silhouette"]
        .agg(["count", "mean", "median", "std", "min", "max"])
        .reset_index()
        .rename(columns={
            "count": "n",
            "mean": "silhouette_mean",
            "median": "silhouette_median"
        })
    )

    global_df = pd.Series({"global_silhouette": sil_mean}, name="silhouette_summary")

    return global_df, cluster_df, sample_df


def pairwise_cluster_consensus(cc_table, cluster_assign):
    """
    Mean/median consensus for each cluster pair.
    Useful to see which clusters are poorly separated.
    """
    cc_table, cluster_assign = _check_inputs(cc_table, cluster_assign)
    clusters = sorted(cluster_assign.dropna().unique())

    rows = []
    for i, c1 in enumerate(clusters):
        idx1 = cluster_assign[cluster_assign == c1].index
        for c2 in clusters[i + 1:]:
            idx2 = cluster_assign[cluster_assign == c2].index
            vals = cc_table.loc[idx1, idx2].values.flatten()
            rows.append({
                "cluster_1": c1,
                "cluster_2": c2,
                "mean_consensus": float(np.mean(vals)),
                "median_consensus": float(np.median(vals)),
                "min_consensus": float(np.min(vals)),
                "max_consensus": float(np.max(vals))
            })

    return pd.DataFrame(rows)


def seed_reproducibility(cluster_runs, method="ari"):
    """
    Compare cluster labels across repeated pyNBS runs.
    cluster_runs: dict of {run_name: pd.Series(sample -> cluster)}
    method: 'ari' or 'nmi'

    This requires you to save labels from repeated runs with different seeds/subsampling.
    """
    from sklearn.metrics import adjusted_rand_score, normalized_mutual_info_score

    run_names = list(cluster_runs.keys())
    rows = []

    for i, r1 in enumerate(run_names):
        s1 = cluster_runs[r1]
        for r2 in run_names[i + 1:]:
            s2 = cluster_runs[r2]
            common = s1.index.intersection(s2.index)
            if len(common) == 0:
                continue

            y1 = s1.loc[common].values
            y2 = s2.loc[common].values

            if method.lower() == "ari":
                score = adjusted_rand_score(y1, y2)
            elif method.lower() == "nmi":
                score = normalized_mutual_info_score(y1, y2)
            else:
                raise ValueError("method must be 'ari' or 'nmi'")

            rows.append({
                "run_1": r1,
                "run_2": r2,
                "n_samples": len(common),
                method.lower(): score
            })

    return pd.DataFrame(rows)


def full_stability_report(cc_table, cluster_assign, outdir=None, prefix="pyNBS"):
    """
    Run the main stability summaries and optionally save them.
    """
    cc_table, cluster_assign = _check_inputs(cc_table, cluster_assign)

    global_cons = summarize_consensus_matrix(cc_table)
    cluster_stab = cluster_stability_table(cc_table, cluster_assign)
    sample_conf = sample_assignment_confidence(cc_table, cluster_assign)
    sil_global, sil_cluster, sil_sample = consensus_silhouette(cc_table, cluster_assign)
    pairwise = pairwise_cluster_consensus(cc_table, cluster_assign)

    report = {
        "global_consensus": global_cons,
        "cluster_stability": cluster_stab,
        "sample_confidence": sample_conf,
        "silhouette_global": sil_global,
        "silhouette_cluster": sil_cluster,
        "silhouette_sample": sil_sample,
        "pairwise_cluster_consensus": pairwise
    }

    if outdir is not None:
        if not os.path.exists(outdir):
            os.makedirs(outdir)

        global_cons.to_csv(
            os.path.join(outdir, "{}_global_consensus.csv".format(prefix)),
            header=True
        )
        cluster_stab.to_csv(
            os.path.join(outdir, "{}_cluster_stability.csv".format(prefix)),
            index=False
        )
        sample_conf.to_csv(
            os.path.join(outdir, "{}_sample_confidence.csv".format(prefix))
        )
        sil_global.to_csv(
            os.path.join(outdir, "{}_silhouette_global.csv".format(prefix)),
            header=True
        )
        sil_cluster.to_csv(
            os.path.join(outdir, "{}_silhouette_cluster.csv".format(prefix)),
            index=False
        )
        sil_sample.to_csv(
            os.path.join(outdir, "{}_silhouette_sample.csv".format(prefix))
        )
        pairwise.to_csv(
            os.path.join(outdir, "{}_pairwise_cluster_consensus.csv".format(prefix)),
            index=False
        )

    return report
