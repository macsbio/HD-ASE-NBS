# -*- coding: utf-8 -*-

# pyNBS_stability.py
# Python 2.7-compatible stability metrics for pyNBS consensus clustering

import os
import numpy as np
import pandas as pd

try:
    from sklearn.metrics import silhouette_samples, silhouette_score
    SKLEARN_SILHOUETTE_AVAILABLE = True
except ImportError:
    SKLEARN_SILHOUETTE_AVAILABLE = False


def _check_inputs(cc_table, cluster_assign):
    if not isinstance(cc_table, pd.DataFrame):
        raise TypeError("cc_table must be a pandas DataFrame")
    if not isinstance(cluster_assign, pd.Series):
        raise TypeError("cluster_assign must be a pandas Series")

    common = cc_table.index.intersection(cc_table.columns).intersection(cluster_assign.index)
    if len(common) == 0:
        raise ValueError("No overlapping sample IDs between cc_table and cluster_assign")

    cc_table = cc_table.loc[common, common].copy()
    cluster_assign = cluster_assign.loc[common].copy()

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
    vals = _upper_triangle_values(cc_table.values)

    return pd.Series({
        "n_samples": int(cc_table.shape[0]),
        "mean_consensus": float(np.mean(vals)),
        "median_consensus": float(np.median(vals)),
        "sd_consensus": float(np.std(vals, ddof=1)),
        "min_consensus": float(np.min(vals)),
        "max_consensus": float(np.max(vals))
    }, name="global_consensus")


def cluster_stability_table(cc_table, cluster_assign):
    """
    Per-cluster within-cluster consensus, between-cluster consensus,
    and separation = within_mean - between_mean
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
                "n": int(n_in),
                "within_mean": np.nan,
                "within_median": np.nan,
                "between_mean": np.nan,
                "between_median": np.nan,
                "separation": np.nan
            })
            continue

        within = cc_table.loc[in_idx, in_idx].values
        within_vals = _upper_triangle_values(within)

        within_mean = float(np.mean(within_vals))
        within_median = float(np.median(within_vals))

        if len(out_idx) > 0:
            between_vals = cc_table.loc[in_idx, out_idx].values.flatten()
            between_mean = float(np.mean(between_vals))
            between_median = float(np.median(between_vals))
            separation = within_mean - between_mean
        else:
            between_mean = np.nan
            between_median = np.nan
            separation = np.nan

        rows.append({
            "cluster": cl,
            "n": int(n_in),
            "within_mean": within_mean,
            "within_median": within_median,
            "between_mean": between_mean,
            "between_median": between_median,
            "separation": separation
        })

    return pd.DataFrame(rows).sort_values("cluster")


def sample_assignment_confidence(cc_table, cluster_assign):
    """
    For each sample:
    - mean consensus to its own cluster
    - highest mean consensus to another cluster
    - margin = own - best_other
    """
    cc_table, cluster_assign = _check_inputs(cc_table, cluster_assign)
    clusters = sorted(cluster_assign.dropna().unique())

    rows = []

    for sample in cc_table.index:
        own_cluster = cluster_assign.loc[sample]
        own_members = cluster_assign[cluster_assign == own_cluster].index
        own_members = own_members.difference([sample])

        if len(own_members) > 0:
            own_mean = float(cc_table.loc[sample, own_members].mean())
        else:
            own_mean = np.nan

        best_other_cluster = np.nan
        best_other_mean = np.nan

        for cl in clusters:
            if cl == own_cluster:
                continue

            members = cluster_assign[cluster_assign == cl].index
            if len(members) == 0:
                continue

            this_mean = float(cc_table.loc[sample, members].mean())

            if pd.isnull(best_other_mean) or (this_mean > best_other_mean):
                best_other_mean = this_mean
                best_other_cluster = cl

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

    out = pd.DataFrame(rows)
    out = out.set_index("sample")
    out = out.sort_values(["cluster", "margin"], ascending=[True, False])

    return out


def consensus_silhouette(cc_table, cluster_assign):
    """
    Silhouette using precomputed distance = 1 - consensus.

    Returns:
    - global_df: Series with global silhouette
    - cluster_df: summary per cluster
    - sample_df: silhouette per sample
    """
    if not SKLEARN_SILHOUETTE_AVAILABLE:
        raise ImportError("scikit-learn silhouette functions are not available in this environment")

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

    cluster_df = sample_df.groupby("cluster")["silhouette"].agg([
        "count", "mean", "median", "std", "min", "max"
    ]).reset_index()

    cluster_df = cluster_df.rename(columns={
        "count": "n",
        "mean": "silhouette_mean",
        "median": "silhouette_median"
    })

    global_df = pd.Series({
        "global_silhouette": float(sil_mean)
    }, name="silhouette_summary")

    return global_df, cluster_df, sample_df


def pairwise_cluster_consensus(cc_table, cluster_assign):
    """
    Mean/median consensus between each pair of clusters.
    Lower values indicate better separation.
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
    Compare repeated clustering runs.

    cluster_runs: dict
        {run_name: pandas Series indexed by sample, values = cluster labels}

    method:
        'ari' or 'nmi'
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
                "n_samples": int(len(common)),
                method.lower(): float(score)
            })

    return pd.DataFrame(rows)


def full_stability_report(cc_table, cluster_assign, outdir=None, prefix="pyNBS", do_silhouette=True):
    """
    Main wrapper.

    Returns a dict with:
    - global_consensus
    - cluster_stability
    - sample_confidence
    - pairwise_cluster_consensus
    - silhouette_global      (if available)
    - silhouette_cluster     (if available)
    - silhouette_sample      (if available)
    """
    cc_table, cluster_assign = _check_inputs(cc_table, cluster_assign)

    global_cons = summarize_consensus_matrix(cc_table)
    cluster_stab = cluster_stability_table(cc_table, cluster_assign)
    sample_conf = sample_assignment_confidence(cc_table, cluster_assign)
    pairwise = pairwise_cluster_consensus(cc_table, cluster_assign)

    report = {
        "global_consensus": global_cons,
        "cluster_stability": cluster_stab,
        "sample_confidence": sample_conf,
        "pairwise_cluster_consensus": pairwise
    }

    if do_silhouette:
        try:
            sil_global, sil_cluster, sil_sample = consensus_silhouette(cc_table, cluster_assign)
            report["silhouette_global"] = sil_global
            report["silhouette_cluster"] = sil_cluster
            report["silhouette_sample"] = sil_sample
        except Exception as e:
            report["silhouette_error"] = str(e)

    if outdir is not None:
        if not os.path.exists(outdir):
            os.makedirs(outdir)

        global_cons.to_csv(
            os.path.join(outdir, "%s_global_consensus.csv" % prefix),
            header=True
        )

        cluster_stab.to_csv(
            os.path.join(outdir, "%s_cluster_stability.csv" % prefix),
            index=False
        )

        sample_conf.to_csv(
            os.path.join(outdir, "%s_sample_confidence.csv" % prefix)
        )

        pairwise.to_csv(
            os.path.join(outdir, "%s_pairwise_cluster_consensus.csv" % prefix),
            index=False
        )

        if "silhouette_global" in report:
            report["silhouette_global"].to_csv(
                os.path.join(outdir, "%s_silhouette_global.csv" % prefix),
                header=True
            )
            report["silhouette_cluster"].to_csv(
                os.path.join(outdir, "%s_silhouette_cluster.csv" % prefix),
                index=False
            )
            report["silhouette_sample"].to_csv(
                os.path.join(outdir, "%s_silhouette_sample.csv" % prefix)
            )

        if "silhouette_error" in report:
            fh = open(os.path.join(outdir, "%s_silhouette_error.txt" % prefix), "w")
            fh.write(report["silhouette_error"])
            fh.close()

    return report
