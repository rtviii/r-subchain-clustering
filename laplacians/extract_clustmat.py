import numpy as np
import json
import os
from namespaces.total import total

# TODO: accessible neigbor-matrix explorer page.


# dynamic arguments:
# - NAMESPACE to contstruct a matrix from
# - BATCH of molecules to run mean on

#
namespace = total
#


def save_as_matrix(adjmatrix=np.array):
    np.savetxt('test.csv', np.around(adjmatrix, decimals=0), fmt='%1d')


def openjson(fullpath=str):
    with open(fullpath, 'r') as infile:
        jsonprofile = json.load(infile)
        return jsonprofile


def extract_adjacency_matrix(pdbid=str, radius=float):
    pdbid = str.upper(pdbid)
    clustersdatapath = './../clusterdata/{}/'.format(pdbid)
    allavailable = os.listdir(clustersdatapath)
    clusters = openjson(clustersdatapath + allavailable[0])

    construct_adjacency_matrix(clusters, namespace)

def tree2tuplearr(nbrtree=dict):
    nbrtuples = []
    for kvpair in list(nbrtree.items()):
        for nbr in kvpair[1]:
            nbrtuples.append([kvpair[0], nbr])
    return nbrtuples
        
    
def construct_adjacency_matrix(clusters=dict, nomenclature_namespace=dict):


    nbrpairs = tree2tuplearr(clusters['nbrtree'])

    # construct len(space) x len(space) matrix of 0s
    # 1s in all the proteins within a cluster
    # dim = len(nomenclature_namespace)
    data = np.array(nbrpairs)
    nodes = np.unique(data)
    noidx = {n: i for i, n in enumerate(nodes)}
    n = nodes.size
    numdata = np.vectorize(noidx.get)(data)
    A = np.zeros((n, n))
    for tail, head in numdata:
        A[tail, head] = 1
    save_as_matrix(A)


    # substrate = np.zeros((dim, dim))
    # save_as_matrix(substrate)

    # for cluster in clusters:
    #     for protein in cluster:
    #         inner_nbrs = [*cluster]
    #         inner_nbrs.remove(protein)

    #         index = get_protein_index(protein)
    #         print("\nproteins's {} index in space is {}".format(protein, index))
    #         for nbr in inner_nbrs:
    #             nbr_ind = get_protein_index(nbr)
    #             print('{} nbr\'s ind is {}'.format(nbr, nbr_ind))
    #             print("filling 1 at {},{}".format(index, nbr_ind))
    #             substrate[index-1][nbr_ind-1] = 1


extract_adjacency_matrix('5myj')
