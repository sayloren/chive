from hw2skeleton import cluster
from hw2skeleton import io
import os

def test_similarity():
    filename_a = os.path.join("data", "276.pdb")
    filename_b = os.path.join("data", "4629.pdb")

    activesite_a = io.read_active_site(filename_a)
    activesite_b = io.read_active_site(filename_b)

    # update this assertion
    assert cluster.compute_similarity(activesite_a, activesite_b) == 0.0

def test_partition_clustering():
    # tractable subset
    pdb_ids = [276, 4629, 10701]

    active_sites = []
    for id in pdb_ids:
        filepath = os.path.join("data", "%i.pdb"%id)
        active_sites.append(io.read_active_site(filepath))

    # update this assertion
    assert cluster.cluster_by_partitioning(active_sites) == []

def test_hierarchical_clustering():
    # tractable subset
    pdb_ids = [276, 4629, 10701]

    active_sites = []
    for id in pdb_ids:
        filepath = os.path.join("data", "%i.pdb"%id)
        active_sites.append(io.read_active_site(filepath))

    # update this assertion
    assert cluster.cluster_hierarchically(active_sites) == []






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
