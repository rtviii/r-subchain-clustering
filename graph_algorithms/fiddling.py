import networkx as nx
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import json
import os
from scipy.sparse import csr_matrix
import scipy.sparse as ssp
matplotlib.use("TkAgg")


def openjson(fullpath=str):
    with open(fullpath, 'r') as infile:
        jsonprofile = json.load(infile)
        return jsonprofile


def tree2tuplearr(nbrtree=dict):
    nbrtuples = []
    for kvpair in list(nbrtree.items()):
        for nbr in kvpair[1]:
            nbrtuples.append([kvpair[0], nbr])
    return nbrtuples


def test():

    pdbid = '5MYJ'
    clustersdatapath = './clusterdata/{}/'.format(pdbid)
    allavailable = os.listdir(clustersdatapath)
    clusters = openjson(clustersdatapath + allavailable[0])
    print("Opened [{}]\n\n\n".format(clustersdatapath+allavailable[0]))

    allchains = clusters['metadata']['allchains']
    nbrs = tree2tuplearr(clusters['nbrtree'])

    # print(clusters['metadata']['allchains'])
    # print("NBRS {}".format(nbrs))

    # ----------------------------------------------------------
    G = nx.Graph()
    G.add_nodes_from(allchains)
    G.add_edges_from(nbrs)

    # print(G.edges())
    # print(G['uS13']['bL31'])
    er = nx.erdos_renyi_graph(100, 0.15)
    ws = nx.watts_strogatz_graph(30, 3, 0.1) 
    K_5 = nx.complete_graph(5)
    pete =nx.petersen_graph()
    # nx.draw(G, with_labels=True,)

    A = nx.to_numpy_matrix(G)
    ScipA = csr_matrix(A)
    print( ScipA )
    perm = ssp.csgraph.reverse_cuthill_mckee(ScipA)

    
    arr_perm = A[perm, perm]
    # plt.spy(arr_perm )
    print(arr_perm.shape)
    # print(cm)
    # print(A)
    plt.matshow(arr_perm)
    plt.show()
