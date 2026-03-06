import pandas as pd
from pyNBS.pyNBS_stability import full_stability_report

cc_table = pd.read_csv("pynbs_results/myjob_cc_matrix.csv", index_col=0)
cluster_assign = pd.read_csv("pynbs_results/myjob_cluster_assignments.csv", index_col=0).iloc[:, 0]
cluster_assign.name = "cluster"

report = full_stability_report(
    cc_table=cc_table,
    cluster_assign=cluster_assign,
    outdir="pynbs_results/stability/",
    prefix="myjob"
)

print(report["global_consensus"])
print(report["cluster_stability"])
print(report["silhouette_global"])