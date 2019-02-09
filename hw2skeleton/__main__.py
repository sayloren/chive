import sys
from .io import read_active_sites, write_clustering, write_mult_clusterings
from .cluster import cluster_by_partitioning, cluster_hierarchically, format_data, make_distance_matrix
import matplotlib.pyplot as plt
import numpy as np
from sklearn import decomposition
import pandas as pd
import seaborn as sns
import scipy.cluster.hierarchy as shc
# Some quick stuff to make sure the program is called correctly
if len(sys.argv) < 4:
    print("Usage: python -m hw2skeleton [-P| -H] <pdb directory> <output file>")
    sys.exit(0)

# read active sites
active_sites = read_active_sites(sys.argv[2])

# set number of clusters
k = 5

# transform the active sites
df_sites = format_data(active_sites)
matrix_sites = make_distance_matrix(df_sites)

# Choose clustering algorithm
if sys.argv[1][0:2] == '-P':
    print("Clustering using Partitioning method")
    clustering = cluster_by_partitioning(matrix_sites,k)
    write_clustering(sys.argv[3], clustering)

if sys.argv[1][0:2] == '-H':
    print("Clustering using hierarchical method")
    clusterings = cluster_hierarchically(matrix_sites,k)
    write_mult_clusterings(sys.argv[3], clusterings)

# plotting
pca = decomposition.PCA(n_components=4)
pc = pca.fit_transform(matrix_sites)
pc_df = pd.DataFrame(data = pc,columns=['PC1', 'PC2','PC3','PC4'])
pc_df['Cluster'] = list(matrix_sites.columns.values)
df = pd.DataFrame({'var':pca.explained_variance_ratio_,'PC':['PC1','PC2','PC3','PC4']})
sns.barplot(x='PC',y="var",data=df, color="c");
sns.lmplot( x="PC1", y="PC2",data=pc_df,fit_reg=False,hue='Cluster',legend=True,scatter_kws={"s": 80}) # specify the point size
shc.dendrogram(shc.linkage(matrix_sites, method='ward'))
plt.show()

# make a random grouping to compare against
# eval clustering alg with davis boulding in
# plot in pc space for similarity
