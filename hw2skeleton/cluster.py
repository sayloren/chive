from .utils import Atom, Residue, ActiveSite
import pandas as pd
import numpy as np
from scipy import spatial
import collections
from sklearn.cluster import AgglomerativeClustering
from scipy.cluster.hierarchy import dendrogram, linkage
from sklearn.metrics.cluster import adjusted_rand_score

# format the active site data into a df with the active sites as
# index and a list of unit vectors for the residues
def format_data(active_sites):
    # active sites contain residues, contain atoms, contain coords

    # collect the residue, max/min coords for each active site
    collect_vector_coords = []
    for s in active_sites:
        for r in s.residues:
            collect_bb = [str(s),r.type]
            collect_a = []
            for a in r.atoms[:3]: # r.atoms[:3] are N, CA, C
                collect_a.append(a.coords)
            x_vals = [item[0] for item in collect_a]
            y_vals = [item[1] for item in collect_a]
            z_vals = [item[2] for item in collect_a]
            collect_vals = [max(x_vals),max(y_vals),max(z_vals),min(x_vals),min(y_vals),min(z_vals)]
            collect_bb.append(collect_vals)
            collect_vector_coords.append(collect_bb)

    # collect the active site, residue, and max/min coords into df
    pd_backbone = pd.DataFrame(collect_vector_coords,columns=['activesite','aminoacid','mmcoords'])
    pd_backbone[['x_max','y_max','z_max','x_min','y_min','z_min']] = pd.DataFrame(pd_backbone.mmcoords.values.tolist(), index= pd_backbone.index)
    pd_backbone.drop('mmcoords',inplace=True,axis=1)

    # get the distance between the max and min for each residue
    pd_backbone['x_dist'] = pd_backbone['x_max']-pd_backbone['x_min']
    pd_backbone['y_dist'] = pd_backbone['y_max']-pd_backbone['y_min']
    pd_backbone['z_dist'] = pd_backbone['z_max']-pd_backbone['z_min']
    pd_backbone.drop(['x_max','x_min','y_max','y_min','z_max','z_min'],inplace=True,axis=1)

    # make a vector for the min/max coords
    pd_backbone['vector'] = pd_backbone[['x_dist','y_dist','z_dist']].values.tolist()
    pd_backbone.drop(['x_dist','y_dist','z_dist'],inplace=True,axis=1)

    # convert to unit vector
    pd_backbone['unit_v'] = pd_backbone['vector'].apply(lambda x: x/np.linalg.norm(x))
    pd_backbone.drop('vector',inplace=True,axis=1)

    # list the residues and max/min coords for each active site
    group_bb = pd_backbone.groupby(['activesite']).agg(lambda x: list(x))
    group_bb.drop(['aminoacid'],inplace=True,axis=1)

    # get average of vectors
    group_bb['ave'] = group_bb['unit_v'].apply(lambda x: [float(sum(col))/len(col) for col in zip(*x)])
    group_ave = group_bb[['ave']]

    return group_ave

# make matrix comparing every active site to every other active site
def make_distance_matrix(df_sites):
    collect_comp = []
    collect_names = []

    # compare every active site distance matrix value to every other
    for indexone, rowone in df_sites.iterrows():
        collect_names.append(indexone)
        collect = []
        site_a = rowone['ave']
        for indextwo, rowtwo in df_sites.iterrows():
            site_b = rowtwo['ave']
            collect.append(compute_similarity(site_a, site_b))
        collect_comp.append(collect)

    df_matrix = pd.DataFrame(collect_comp,columns=collect_names,index=collect_names)
    return df_matrix

# compute the cosine similarity between each pair of active sites
def compute_similarity(site_a, site_b):
    """
    Compute the similarity between two given ActiveSite instances.

    Input: two ActiveSite instances
    Output: the similarity between them (a floating point number)
    """
    # get cosine similarity between each pair of unit vectors from each list
    # cosine simliarity computes distance, subtract from 1 for similarity
    # 0.0 # 0 if no similarity, 1 if the same
    similarity = 1 - spatial.distance.cosine(site_a, site_b)
    return similarity

# cluster the data to closest value in centers
def create_partition_clusters(df,centers,k):
    clusters = [[] for i in range(k)]
    clusters_vals = [[] for i in range(k)]

    # iterate through each active site, get the value from the matrix
    # corresponding to the center and append to cluster
    for index, row in df.iterrows():
        get_values = [df.loc[c,index] for c in centers]
        max_index = get_values.index(max(get_values))
        max_val = max(get_values)
        clusters[max_index].append(index)
        clusters_vals[max_index].append(max_val)

    # for each cluster, get the median value and return as the new centers
    new_centers = []
    for c,d in zip(clusters,centers):
        new_c = sorted(c)
        index = int((len(new_c)-1)/2)
        new_centers.append(new_c[index])

    return clusters,new_centers

# k means
def cluster_by_partitioning(matrix_sites,k):
    """
    Cluster a given set of ActiveSite instances using a partitioning method.

    Input: a list of ActiveSite instances
    Output: a clustering of ActiveSite instances
            (this is really a list of clusters, each of which is list of
            ActiveSite instances)
    """

    # initialize the clusters randomly
    np.random.seed(3)
    names = matrix_sites.index.tolist()
    centers = [np.random.choice(names) for i in range(k)]
    test = np.empty(k)

    # while the centers are not the same as the previous iteration, cluster
    count = 0
    while (collections.Counter(centers) != collections.Counter(test) or count < 20):
        test = centers
        clusters,centers = create_partition_clusters(matrix_sites,centers,k)
        count += 1

    return clusters

# agglomerative clustering
def cluster_hierarchically(matrix_sites,k):
    """
    Cluster the given set of ActiveSite instances using a hierarchical algorithm.                                                                  #

    Input: a list of ActiveSite instances
    Output: a list of clusterings
            (each clustering is a list of lists of Sequence objects)
    """

    # make the list of clusters to return
    clusters = [[] for i in range(k)]

    # cluster using the agglomerative method;
    # each active site starts as its own cluster,
    # merge each pair of active sites into the same cluster if
    # they have the smallest distance between the closest points (single)
    # in euclidean space
    cluster = AgglomerativeClustering(n_clusters=k, affinity='euclidean', linkage='single')
    cluster.fit_predict(matrix_sites)
    cluster_labels = cluster.labels_
    site_labels = list(matrix_sites.columns.values)

    # using the index of cluster labels created, get the active site names to
    # populate the return clustered list
    for s,c in zip(site_labels,cluster_labels):
        clusters[c].append(s)

    return clusters

# random clustering
def cluster_randomly(matrix_sites,k):
    """
    Cluster the given set of ActiveSite instances randomly.                                                                  #

    Input: a list of ActiveSite instances
    Output: a list of clusterings
            (each clustering is a list of lists of Sequence objects)
    """
    # make empty list with the number of clusters desired
    clusters = [[] for i in range(k)]

    # collect the active site labels
    site_labels = list(matrix_sites.columns.values)

    # randomly make a list of the same size as the number of active sites
    # labels, with a random index for one of the k number of clusters
    cluster_labels = np.random.choice(k, len(site_labels))

    # populate the clusters
    for s,c in zip(site_labels,cluster_labels):
        clusters[c].append(s)

    return clusters

# this is code from when I tried to make my own hierarchical clustering alg,
# but couldn't get it put together in time

    # # get the average distance from the matrix
    # # average_distance = matrix_sites.sum().sum()/len(matrix_sites.index)
    #
    # # 1 - each point is its own cluster
    # # 2 - merge points that are closest (use the centroid of the cluster)
    # # 3 - stop when get k clusters | when the max distance is reached
    #
    # # clusters = [[] for i in range(k)] # make k empty clusters
    # clusters = [[]]
    # # while [not i for i in clusters]: # while any cluster is empty
    #
    # # while len(clusters) < k:
    # while any(matrix_sites.values.flatten() < 1.1):
    #     print(clusters)
    #     i = matrix_sites.min().idxmin()
    #     j = matrix_sites.idxmin()[i]
    #     print(i,j)
    #     # smallest_value = matrix_sites.loc[i,j]
    #
    #
    #     list_index_i = next((c for c,v in enumerate(clusters) if i in v), -1)
    #     list_index_j = next((c for c,v in enumerate(clusters) if j in v), -1)
    #     if list_index_i != -1:
    #         clusters[list_index_i].append(j)
    #     elif list_index_j != -1:
    #         clusters[list_index_j].append(i)
    #     else:
    #         clusters.append([i,j])
    #
    #     matrix_sites.loc[i,j] = 1.1
    #     matrix_sites.loc[j,i] = 1.1
    #     print(matrix_sites.loc[i,j])
    #
    # print(clusters)
    # # pd_sort = matrix_sites.sort_values(i).sort_values(j,axis=1)
