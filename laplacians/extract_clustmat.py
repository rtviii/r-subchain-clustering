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


def construct_adjacency_matrix(clusters=dict, nomenclature_namespace=dict):
    a = list( nomenclature_namespace.keys() )
    no_dupes = [x for n, x in enumerate(a) if x not in a[:n]]
    print(no_dupes)  # [[1], [2], [3], [5]]

    dupes = [x for n, x in enumerate(a) if x in a[:n]]
    print(dupes)  # [[1], [3]]j

    # construct len(space) x len(space) matrix of 0s
    # 1s in all the proteins within a cluster
    # dim = len(nomenclature_namespace)
    allnames = list(nomenclature_namespace.keys())
    data = np.array(allnames)
    nodes = np.unique(data)
    noidx = {n: i for i, n in enumerate(nodes)}
    print(noidx)
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
