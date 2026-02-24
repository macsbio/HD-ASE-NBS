# pynbs_parallel.py  (Python 2.7)

import os
# Prevent BLAS/OpenMP oversubscription when using multiprocessing
os.environ["OMP_NUM_THREADS"] = "1"
os.environ["MKL_NUM_THREADS"] = "1"
os.environ["OPENBLAS_NUM_THREADS"] = "1"
os.environ["NUMEXPR_NUM_THREADS"] = "1"

from pyNBS import data_import_tools as dit
from pyNBS import network_propagation as prop
from pyNBS import pyNBS_core as core
from pyNBS import pyNBS_single
from pyNBS import consensus_clustering as cc
from pyNBS import pyNBS_plotting as plot

import time
import pandas as pd
import numpy as np
import multiprocessing as mp
from IPython.display import Image

# -------- worker plumbing --------
def _init_worker(_sm_mat, _knnGlap, _network, _kernel, _clusters):
    global sm_mat, knnGlap, network, kernel, clusters
    sm_mat = _sm_mat
    knnGlap = _knnGlap
    network = _network
    kernel = _kernel
    clusters = _clusters

def _run_one(i):
    # Optional: make randomness reproducible per-iteration
    np.random.seed(12345 + i)

    t0 = time.time()
    H = pyNBS_single.NBS_single(sm_mat, knnGlap, propNet=network, propNet_kernel=kernel, k=clusters)
    dt = time.time() - t0
    print("NBS iteration: %d complete: %.2f seconds" % (i + 1, dt))
    return H

def main():
    sm_data_filepath = 'mutationmatrix_900_wasp.txt'
    network_filepath = 'ppi_network_900_wasp.txt'

    sm_mat_local = dit.load_binary_mutation_data(sm_data_filepath, filetype='list', delimiter='\t')
    network_local = dit.load_network_file(network_filepath)

    outdir = 'pynbs_results/'
    if not os.path.exists(outdir):
        os.makedirs(outdir)
    job_name = 'HTT_900_wasp'
    save_args = {'outdir': outdir, 'job_name': job_name}

    # Graph Laplacian regulariser
    knnGlap_local = core.network_inf_KNN_glap(network_local, **save_args)

    # Network node order
    network_nodes = list(network_local.nodes())

    # IMPORTANT: make sure mutation matrix has a column for every network node (missing -> 0)
    sm_mat_local = sm_mat_local.reindex(columns=network_nodes).fillna(0).astype(int)

    # Identity matrix
    network_I = pd.DataFrame(np.identity(len(network_nodes)), index=network_nodes, columns=network_nodes)

    alpha = 0.7
    save_args['iteration_label'] = 'kernel'

    # Build propagation kernel (your script overwrote kernel; here we keep only the identity-kernel)
    kernel_local = prop.network_propagation(network_local, network_I, alpha=alpha, symmetric_norm=True)

    clusters_local = 2
    niter = 1000

    nproc = 12  # set cores here

    # Parallel map
    pool = mp.Pool(
        processes=nproc,
        initializer=_init_worker,
        initargs=(sm_mat_local, knnGlap_local, network_local, kernel_local, clusters_local),
        maxtasksperchild=50  # helps if there are small memory leaks
    )

    Hlist = []
    try:
        # chunksize reduces IPC overhead for large niter
        for H in pool.imap_unordered(_run_one, range(niter), chunksize=5):
            Hlist.append(H)
    finally:
        pool.close()
        pool.join()

    # Consensus clustering + plot
    NBS_cc_table, NBS_cc_linkage, NBS_cluster_assign = cc.consensus_hclust_hard(Hlist, k=clusters_local, **save_args)
    clust_cmap = plot.cluster_color_assign(NBS_cluster_assign, name='pyNBS HTT Cluster Assignments')
    plot.plot_cc_map(NBS_cc_table, NBS_cc_linkage, col_color_map=clust_cmap, **save_args)
    Image(filename=save_args['outdir'] + save_args['job_name'] + '_cc_map.png', width=600, height=600)

if __name__ == '__main__':
    main()
