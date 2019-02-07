from .utils import Atom, Residue, ActiveSite
import pandas as pd
import numpy as np
from scipy import spatial
# import itertools
import collections
# import math

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
            for a in r.atoms:
                collect_a.append(a.coords) # r.atoms[:3] are N, CA, C
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
def create_clusters(df,centers,k):
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
def cluster_by_partitioning(active_sites):
    """
    Cluster a given set of ActiveSite instances using a partitioning method.

    Input: a list of ActiveSite instances
    Output: a clustering of ActiveSite instances
            (this is really a list of clusters, each of which is list of
            ActiveSite instances)
    """
    df_sites = format_data(active_sites)
    matrix_sites = make_distance_matrix(df_sites)

    # initialize the clusters randomly
    k = 5
    np.random.seed(3)
    names = matrix_sites.index.tolist()
    centers = [np.random.choice(names) for i in range(k)]
    test = np.empty(k)

    # while the centers are not the same as the previous iteration
    # creat the clusters
    while collections.Counter(centers) != collections.Counter(test):
        test = centers
        clusters,centers = create_clusters(matrix_sites,centers,k)

    return clusters

# ward hierarchical
def cluster_hierarchically(active_sites):
    """
    Cluster the given set of ActiveSite instances using a hierarchical algorithm.                                                                  #

    Input: a list of ActiveSite instances
    Output: a list of clusterings
            (each clustering is a list of lists of Sequence objects)
    """

    # Fill in your code here!

    return []
