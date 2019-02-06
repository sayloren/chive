import sys
from .io import read_active_sites, write_clustering, write_mult_clusterings
from .cluster import cluster_by_partitioning, cluster_hierarchically
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import numpy as np

# Some quick stuff to make sure the program is called correctly
if len(sys.argv) < 4:
    print("Usage: python -m hw2skeleton [-P| -H] <pdb directory> <output file>")
    sys.exit(0)

active_sites = read_active_sites(sys.argv[2])

# Choose clustering algorithm
if sys.argv[1][0:2] == '-P':
    print("Clustering using Partitioning method")
    clustering = cluster_by_partitioning(active_sites)
    write_clustering(sys.argv[3], clustering)

if sys.argv[1][0:2] == '-H':
    print("Clustering using hierarchical method")
    clusterings = cluster_hierarchically(active_sites)
    write_mult_clusterings(sys.argv[3], clusterings)

# make a random grouping to compare against
# eval clustering alg with davis boulding in
# plot in pc space for similarity
