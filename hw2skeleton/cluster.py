from .utils import Atom, Residue, ActiveSite
import pandas as pd
import numpy as np

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
    # get angle between each residue

    # return np.arccos(np.clip(np.dot(v1_u, v2_u), -1.0, 1.0))
    print(pd_backbone.head())
    return pd_backbone


def compute_similarity(site_a, site_b):
    """
    Compute the similarity between two given ActiveSite instances.

    Input: two ActiveSite instances
    Output: the similarity between them (a floating point number)
    """
    similarity = 0.0 # 0 if no similarity, 1 if the same
    # cosine similarity

    # jaccard similarity
    # kabash algorithm

    return similarity


def cluster_by_partitioning(active_sites):
    """
    Cluster a given set of ActiveSite instances using a partitioning method.

    Input: a list of ActiveSite instances
    Output: a clustering of ActiveSite instances
            (this is really a list of clusters, each of which is list of
            ActiveSite instances)
    """

    # len = active_sites.len()
    # k = 5
    format_data(active_sites)



    # compare the similarites

    # pick random starting points for k clusters
    # assign to nearest cluster centroid
    # new cluster center by taking the average of the assigned points
    # recursive 2/3 until non change (bistable? - or cap)





    return []


def cluster_hierarchically(active_sites):
    """
    Cluster the given set of ActiveSite instances using a hierarchical algorithm.                                                                  #

    Input: a list of ActiveSite instances
    Output: a list of clusterings
            (each clustering is a list of lists of Sequence objects)
    """

    # Fill in your code here!

    return []
