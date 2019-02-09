from hw2skeleton import cluster
from hw2skeleton import io
import os
import pandas as pd
import numpy as np

def test_similarity():
    active_sites = io.read_active_sites("data")

    # make the similarity matrix
    df_sites = cluster.format_data(active_sites)
    matrix_sites = cluster.make_distance_matrix(df_sites)

    # get the sum of diagonal
    sum_diagonal = int(pd.Series(np.diag(matrix_sites), index=[matrix_sites.index, matrix_sites.columns]).sum())

    # get the length of number of sites in the matrix
    num_sites = len(matrix_sites.index)

    assert sum_diagonal == num_sites

# def test_partition_clustering():
#     # tractable subset
#     pdb_ids = [276, 4629, 10701]
#
#     # active_sites = io.read_active_sites("data")
#
#     active_sites = []
#     for id in pdb_ids:
#         filepath = os.path.join("data", "%i.pdb"%id)
#         active_sites.append(io.read_active_site(filepath))
#
#     # make the similarity matrix
#     df_sites = cluster.format_data(active_sites)
#     matrix_sites = cluster.make_distance_matrix(df_sites)
#     k = 1
#
#     # update this assertion
#     assert sorted(cluster.cluster_by_partitioning(df_sites,k)) == sorted([276, 4629, 10701])
#
# def test_hierarchical_clustering():
#     # tractable subset
#     pdb_ids = [276, 4629, 10701]
#
#     active_sites = []
#     for id in pdb_ids:
#         filepath = os.path.join("data", "%i.pdb"%id)
#         active_sites.append(io.read_active_site(filepath))
#
#     df_sites = format_data(active_sites)
#     matrix_sites = make_distance_matrix(df_sites)
#     k = 1
#
#     update this assertion
#     assert sorted(cluster.cluster_hierarchically(df_sites,k)) == sorted([276, 4629, 10701])

# make a random grouping to compare against
# eval clustering alg with davis boulding in

# plotting
# pca = decomposition.PCA(n_components=4)
# pc = pca.fit_transform(matrix_sites)
# pc_df = pd.DataFrame(data = pc,columns=['PC1', 'PC2','PC3','PC4'])
# pc_df['Cluster'] = list(matrix_sites.columns.values)
# df = pd.DataFrame({'var':pca.explained_variance_ratio_,'PC':['PC1','PC2','PC3','PC4']})
# sns.barplot(x='PC',y="var",data=df, color="c");
# sns.lmplot( x="PC1", y="PC2",data=pc_df,fit_reg=False,hue='Cluster',legend=True,scatter_kws={"s": 80}) # specify the point size
# plt.show()

# rand_score = adjusted_rand_score(clusterings,)
