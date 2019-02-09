import matplotlib.pyplot as plt
from .cluster import cluster_by_partitioning, cluster_hierarchically, cluster_randomly
from hw2skeleton import io
from sklearn.metrics.cluster import adjusted_rand_score

def run_similarity_evals(matrix_sites,k):
    part_clustering = cluster_by_partitioning(matrix_sites,k)
    hier_clusterings = cluster_hierarchically(matrix_sites,k)
    rand_clusterings = cluster_randomly(matrix_sites,k)

    # compare my partitioning and hierarchical clustering methods
    # rand_score = adjusted_rand_score()
    # print(rand_score)



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
