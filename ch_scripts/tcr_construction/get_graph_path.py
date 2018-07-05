import pandas as pd
import numpy as np
import os


def levensteindist(seq1, seq2):
    return sum(1 if seq1[i] != seq2[i] else 0 for i in range(len(seq1)))


def maxdist_finder(cdr, cdrs):
    maxdist = 0
    for compare_cdr in cdrs:
        dist = levensteindist(cdr, compare_cdr)
        if maxdist < dist:
            maxdist = dist
    return maxdist


def find_centroids(cluster_path="path_finder/clusters.txt", outpath="path_finder/centroids.txt"):
    centroids = pd.DataFrame()
    clusters = pd.read_table(cluster_path)
    for cluster, cdrs in clusters.groupby('cid'):
        cdrs = cdrs.reset_index(drop=True)
        cdrs['maxdist'] = cdrs['cdr3aa'].apply(maxdist_finder, cdrs=list(cdrs['cdr3aa']))
        centroids = centroids.append(cdrs.iloc[[cdrs['maxdist'].idxmin('maxdist')]])

    centroids.reset_index(drop=True).to_csv(outpath, sep='\t', index=None)


def find_graph_path(cluster_path="path_finder/clusters.txt",
                    inpgraph_path="path_finder/graph.txt",
                    centroids_path="path_finder/centroids.txt",
                    outputinfo_path="path_finder/graph_path.info.txt",
                    output_path="path_finder/graph_path.txt"):

    clusters = pd.read_table(cluster_path, usecols=['cdr3aa', 'epitope', 'gene', 'cid'])

    graph = pd.read_table(inpgraph_path)

    graph["cdr3aa"] = graph["from.cdr3aa"]

    graph = graph.merge(clusters, on=["epitope", "gene", "cdr3aa"], how='outer')

    centroids = pd.read_table(centroids_path)

    cluster_num = 0
    threshold_size = 5
    sorted_graph = pd.DataFrame()
    for cluster, cdrs in graph.groupby('cid'):
        if len(cdrs) < threshold_size:
            continue
        cluster_num += 1
        used_cdrs = set()
        k = 0
        work_cdrs = [centroids[centroids["cid"] == cluster]["cdr3aa"].iloc[0]]
        unused_cdrs = cdrs

        while k != 1:
            for work_cdr in work_cdrs:
                if work_cdr not in used_cdrs:
                    used_cdrs.add(work_cdr)
                    work_nodes = unused_cdrs[unused_cdrs["from.cdr3aa"] == work_cdr]
                    sorted_graph = sorted_graph.append(work_nodes)
                    work_cdrs += list(work_nodes["to.cdr3aa"])
                    unused_cdrs = unused_cdrs[(unused_cdrs["from.cdr3aa"] != work_cdr) & (~unused_cdrs["to.cdr3aa"].isin(work_cdrs))]
            if len(unused_cdrs) == 0:
                k = 1

    with open(outputinfo_path, "w") as out:
        out.write("Number of clusters with more than {} tcrs each: {}\n".format(threshold_size, cluster_num))
        out.write("Number of alpha clusters: {}\n".format(len(sorted_graph[sorted_graph["gene"]=="TRA"]["cid"].unique())))
        out.write("Number of beta clusters: {}\n".format(len(sorted_graph[sorted_graph["gene"] == "TRB"]["cid"].unique())))

    sorted_graph.to_csv(output_path, sep="\t", index=None, columns=["epitope", "gene",
                                                                                     "from.cdr3aa",	"to.cdr3aa",
                                                                                     "from.v", "from.j",
                                                                                     "to.v", "to.j",
                                                                                     "cid"])


def pipeline(folder="path_finder", overwrite=False):
    cluster_path = os.path.join(folder, "clusters.txt")
    inpgraph_path = os.path.join(folder, "graph.txt")
    centroids_path = os.path.join(folder, "centroids.txt")
    graph_info_path = os.path.join(folder, "graph_path.info.txt")
    graph_path = os.path.join(folder, "graph_path.txt")
    if overwrite is True or not os.path.exists(centroids_path):
        print("Search for centroids\n")
        find_centroids(cluster_path, centroids_path)
        print("Done\n")
    if overwrite is True or not os.path.exists(graph_path) or not os.path.exists(graph_info_path):
        print("Search for path in graph\n")
        find_graph_path(cluster_path, inpgraph_path, centroids_path, graph_info_path, graph_path)
        print("Done\n")
