from pyNBS import pyNBS_plotting as plot
import os
import pandas as pd
import numpy as np
from scipy.spatial.distance import squareform
from scipy.cluster.hierarchy import linkage

outdir = "pynbs_results/"
job_name = "HTT_900_wasp_replot"
save_args = {"outdir": outdir, "job_name": job_name}
if not os.path.exists(outdir):
    os.makedirs(outdir)

# 1) consensus matrix (csv)
cc = pd.read_csv("HTT_900_wasp_cc_matrix.csv", index_col=0)
cc = cc.loc[cc.index, cc.index].astype(float)  # enforce square + order

# 2) linkage from consensus (similarity -> distance)
dist = 1.0 - cc.values
np.fill_diagonal(dist, 0.0)
Z = linkage(squareform(dist, checks=False), method="average")

# 3) cluster assignments (csv)
a = pd.read_csv("HTT_900_wasp_cluster_assignments.csv")

# expected format: first col = sample, second col = cluster (or columns named sample/cluster)
if "sample" in a.columns and "cluster" in a.columns:
    assign = pd.Series(a["cluster"].values, index=a["sample"].astype(str))
else:
    assign = pd.Series(a.iloc[:, 1].values, index=a.iloc[:, 0].astype(str))

assign = assign.reindex(cc.index)
cmap = plot.cluster_color_assign(assign, name="pyNBS Cluster Assignments")

# 4) plot
plot.plot_cc_map(cc, Z, col_color_map=cmap, **save_args)
print("Wrote:", os.path.join(outdir, job_name + "_cc_map.png"))
