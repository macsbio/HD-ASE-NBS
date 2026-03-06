# -*- coding: utf-8 -*-

import pandas as pd
from pyNBS.pyNBS_stability import full_stability_report
from pyNBS.pyNBS_stability import save_all_stability_plots

cc_file = "pynbs_results/YOUR_CC_MATRIX.csv"
cluster_file = "pynbs_results/YOUR_CLUSTER_ASSIGNMENTS.csv"
prefix = "YOUR_PREFIX"

cc_table = pd.read_csv(cc_file, index_col=0)

cluster_assign = pd.read_csv(
    cluster_file,
    index_col=0,
    squeeze=True
)
cluster_assign.name = "cluster"

report = full_stability_report(
    cc_table=cc_table,
    cluster_assign=cluster_assign,
    outdir="pynbs_results/stability",
    prefix=prefix,
    do_silhouette=True
)

save_all_stability_plots(
    report=report,
    outdir="pynbs_results/stability",
    prefix=prefix
)

print report["global_consensus"]
print report["cluster_stability"]

if "silhouette_global" in report:
    print report["silhouette_global"]
else:
    print "Silhouette not available."
    if "silhouette_error" in report:
        print report["silhouette_error"]
