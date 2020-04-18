import numpy as np
import json
import os
import pandas as pd
from pprint import pprint as pp
from namespaces.total import total
import matplotlib.pyplot as plt
import collections
from scipy.sparse.csgraph import laplacian
import matplotlib
matplotlib.use("TkAgg")

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
    nom = list(nomenclature_namespace.keys())

    # data = pd.DataFrame(0, index=nom, columns=nom, dtype=int)
    # for x, y in nbrpairs:
    #     if (x not in nom or y not in nom):
    #         pass
    #     else:
    #         data.loc[x, y] = 1

    


    def degreeMatrix(adjmatrix):
        n = np.size(adjmatrix, 1)
        D = np.zeros((n, n))

        for row in adjmatrix:
            for col in adjmatrix:
                D[col, row] = D[col, row] + 1
        return D

    D = degreeMatrix(np.array(data))
    plt.matshow(D)
    plt.show()

    # # print("DATA",data)
    # substrate = np.unique(nom)
    # noidx = {n: i for i, n in enumerate(nom)}
    # print("Substrate", noidx)

    # def nbrIncidenceIndex(nbrpair=tuple):
    #     print("Inc index got", nbrpair)
    #     return [noidx.get(nbrpair[0]), noidx.get(nbrpair[1])]

    # print(numdata)
    # # print ("Numdata : \n", numdata)
    # A = np.zeros((n, n))
    # for tail, head in numdata:
    #     A[tail, head] = 1
    # print(A)

    # save_as_matrix(A)


# if not (os.path.isdir('./clusterdata/{}/'.format(str.upper(pdbid)))):
#     Path('./clusterdata/{}/'.format(pdbid)
#             ).mkdir(parents=True, exist_ok=True)
# with open(filename, 'w') as out:
#     json.dump(subchainclusters, out)
#     print('Saved successfully at \t [{}]'.format(filename))
# return subchainclusters

extract_adjacency_matrix('5myj')
