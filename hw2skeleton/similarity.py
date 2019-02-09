import matplotlib.pyplot as plt
from .cluster import cluster_by_partitioning, cluster_hierarchically, cluster_randomly
from hw2skeleton import io
from sklearn.metrics.cluster import adjusted_rand_score
from sklearn.metrics import silhouette_score
import pandas as pd
import seaborn as sns
import pathlib
import numpy as np
from matplotlib_venn import venn2

def run_similarity_evals(matrix_sites,k):

    # run clustering algs
    part_clusterings = cluster_by_partitioning(matrix_sites,k)
    hier_clusterings = cluster_hierarchically(matrix_sites,k)
    rand_clusterings = cluster_randomly(matrix_sites,k)

    # get silhouette score for each cluster alg
    print(silhouette_score(matrix_sites,part_clusterings))
    print(silhouette_score(matrix_sites,hier_clusterings))
    print(silhouette_score(matrix_sites,rand_clusterings))

    # compare my partitioning and hierarchical clustering methods
    print(adjusted_rand_score(part_clusterings,hier_clusterings))
    print(adjusted_rand_score(part_clusterings,rand_clusterings))
    print(adjusted_rand_score(hier_clusterings,rand_clusterings))

    # get the histograms for the number of active sites in each clusters per method
    fig = plt.figure(figsize=(10,10))
    sns.distplot(part_clusterings,color='#f21e5e',label='partition')
    sns.distplot(hier_clusterings,color='#d409cc',label='hierarchy')
    sns.distplot(rand_clusterings,color='#e7df1b',label='random')
    plt.title('Clustering with k = {0}'.format(k))
    sns.despine()
    plt.legend()

    outdir = pathlib.Path('images')
    outfile = outdir / "Cluster_hist.png"
    outdir.mkdir(parents=True, exist_ok=True)
    plt.savefig(str(outfile),format='png')
    plt.close()

    site_labels = np.array(list(matrix_sites.columns.values))

    # get the largest cluster for each method
    largest_part = np.argmax(np.bincount(part_clusterings))
    largest_hier = np.argmax(np.bincount(hier_clusterings))
    largest_rand = np.argmax(np.bincount(rand_clusterings))

    # the list indeces where the largest cluster values are
    loc_part = [i for i, e in enumerate(part_clusterings) if e == largest_part]
    loc_hier = [i for i, e in enumerate(hier_clusterings) if e == largest_hier]
    loc_rand = [i for i, e in enumerate(rand_clusterings) if e == largest_rand]

    # the active sites that are in the largest group
    active_part = list(site_labels[loc_part])
    active_hier = list(site_labels[loc_hier])
    active_rand = list(site_labels[loc_rand])

    # get the number of elements shared by the largest clusters between hierarchical and partition clustering
    fig2 = plt.figure(figsize=(10,10))
    venn2([set(active_part), set(active_hier)],set_labels=('Partitioning','Hierarchical'))
    plt.title('Overlap between largest clusters')
    outdir = pathlib.Path('images')
    outfile = outdir / "Cluster_venn.png"
    outdir.mkdir(parents=True, exist_ok=True)
    plt.savefig(str(outfile),format='png')
    plt.close()

    print(active_part)

    # plotting
    # pca = decomposition.PCA(n_components=4)
    # pc = pca.fit_transform(matrix_sites)
    # pc_df = pd.DataFrame(data = pc,columns=['PC1', 'PC2','PC3','PC4'])
    # pc_df['Cluster'] = list(matrix_sites.columns.values)
    # df = pd.DataFrame({'var':pca.explained_variance_ratio_,'PC':['PC1','PC2','PC3','PC4']})
    # sns.barplot(x='PC',y="var",data=df, color="c");
    # sns.lmplot( x="PC1", y="PC2",data=pc_df,fit_reg=False,hue='Cluster',legend=True,scatter_kws={"s": 80}) # specify the point size
    # plt.show()
