import networkx as nx
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import json
import os
from scipy.sparse import csr_matrix
from namespaces.total import total
import functools
import seaborn as sns
import scipy.sparse as ssp
from casadi import * 

matplotlib.use("TkAgg")
namespace = total

# TARGET GROUPS:

bacteria = ['5NJT',"5AFI", "5JVG", "5MYJ","5O60","5V7Q","5NRG","4Y4P","3J7Z", "5VP2"]
eukarya  = ["6EK0","5T2A", "3J79","5GAK","4V7E","5T5H","5XXB", "5XY3", "4UGO"]

#tight: 6eko, 5t5h

#large: 4v7e(e), 

def construct_adjacency_matrix(clusters=dict, nomenclature_namespace=dict, verbose=False):

    nbrpairs  = tree2tuplearr(clusters['nbrtree'])
    dim       = len(nomenclature_namespace.items())
    keys      = list(nomenclature_namespace.keys())
    substrate = np.zeros((dim, dim))

    for pair in nbrpairs:
        if (pair[0] not in keys or pair[1] not in keys):
            pass
        else:
            substrate[nomenclature_namespace[pair[0]],
                      nomenclature_namespace[pair[1]]] = 1

    if(verbose):
        plt.matshow(substrate)
        plt.show()
    return substrate


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



def get_laplacian_ndarray(nparray):
    csr = csr_matrix(nparray)
    L   = ssp.csgraph.laplacian(csr)
    return L.toarray()



def get_cuthill_mckee_ndarray(nparray):

    csr = csr_matrix(nparray)
    cuthillmckee   = ssp.csgraph.reverse_cuthill_mckee(csr)
    return cuthillmckee.toarray()


def get_simplemean(targetgroup=str):
    targetspath = './../clusterdata/targetgroups/{}/'.format(targetgroup)
    batch       = os.listdir(targetspath)
    print("BATCH >>>>", batch)
    adjmats     = [construct_adjacency_matrix(openjson(targetspath+x), namespace, False) for x in batch]
    
    # laplacians  = []
    reduced = functools.reduce(lambda x, y: np.add(x,y), adjmats)
    reduced = np.divide(reduced, len(adjmats))

    L = get_laplacian_ndarray(reduced)
    print(L)

    sns.clustermap(L)
    # plt.matshow(L)
    plt.show()
    return


get_simplemean('eukarya')
