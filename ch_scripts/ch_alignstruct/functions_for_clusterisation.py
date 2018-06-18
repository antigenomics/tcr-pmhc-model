import os

import pandas as pd
import numpy as np

from matplotlib import pyplot as plt
from scipy.cluster.hierarchy import *
import scipy.spatial.distance as ssd

class filename_maker():
    """
    just class for creating fixed names for files
    """
    def __init__(self, chain, segment, specie):
        self.chain = chain
        self.segment = segment
        self.specie = specie
        self.workdir = ''

    def add_dir(self, workdir):
        self.workdir = workdir

    def get_name(self, name):
        return os.path.join(self.workdir, 'TR{}{}_{}.{}.txt'.format(self.chain, self.segment, self.specie, name))


def find_distances(dists, segment, specie, chain):
    w_dists = dists[(dists['segment'] == segment) & (dists['species']==specie)]
    w_dists = pd.DataFrame(w_dists[w_dists['gene'] == 'TR{}'.format(chain.upper())], columns=['id.1', 'id.2', 'segm.score'])
    w_dists['segm.score'] = w_dists['segm.score']*-1
    w_dists = w_dists.pivot(index='id.1', columns='id.2', values='segm.score')
    return w_dists


def compute_clusters(dist, cutoff=350, method='ward', use_plot=True, figsize_x=8, figsize_y=14):
    work_index = dist.index
    work_Z = linkage(ssd.squareform(dist), method)
    work_clusters = fcluster(work_Z, cutoff, criterion='distance')
    df_clusters = pd.DataFrame(sorted(list(zip(work_clusters, work_index)),
                                      key=lambda x: x[0]), columns=['cluster', 'id'])
    if use_plot is True:
        fig, ax = plt.subplots(figsize=(figsize_x, figsize_y))
        ax = dendrogram(work_Z,
                        color_threshold=cutoff,
                        leaf_font_size=10,
                        labels = work_index,
                        orientation="left")
        plt.show()
    return df_clusters