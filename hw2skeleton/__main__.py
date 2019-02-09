import sys
from .io import read_active_sites, write_clustering, write_mult_clusterings
from .cluster import cluster_by_partitioning, cluster_hierarchically, format_data, make_distance_matrix,cluster_randomly
import numpy as np
from sklearn import decomposition
import pandas as pd
import seaborn as sns
import scipy.cluster.hierarchy as shc
from .similarity import run_similarity_evals

# Some quick stuff to make sure the program is called correctly
if len(sys.argv) < 4:
    print("Usage: python -m hw2skeleton [-P| -H | -R] <pdb directory> <output file>")
    sys.exit(0)

# read active sites
active_sites = read_active_sites(sys.argv[2])

# set number of clusters
k = 3

# transform the active sites
df_sites = format_data(active_sites)
matrix_sites = make_distance_matrix(df_sites)

# Choose clustering algorithm
if sys.argv[1][0:2] == '-P':
    print("Clustering using Partitioning method")
    clustering = cluster_by_partitioning(matrix_sites,k)
    # write_clustering(sys.argv[3], clustering)

if sys.argv[1][0:2] == '-H':
    print("Clustering using hierarchical method")
    clusterings = cluster_hierarchically(matrix_sites,k)
    # write_mult_clusterings(sys.argv[3], clusterings)
    # shc.dendrogram(shc.linkage(matrix_sites, method='ward'))

if sys.argv[1][0:2] == '-R':
    print("Clustering using randomly")
    clusterings = cluster_randomly(matrix_sites,k)
    # write_mult_clusterings(sys.argv[3], clusterings)

run_similarity_evals(matrix_sites,k)
