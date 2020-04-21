import numpy as np
import matplotlib.pyplot as plt
from matplotlib.pyplot import figure
import matplotlib
import json
import os
from scipy.sparse import csr_matrix
from total import total, totalArr
from helpers import openjson, construct_adjacency_matrix, tree2tuplearr, get_cuthill_mckee_ndarray, get_laplacian_ndarray
import subunits_namespace as names
import functools
import seaborn as sns
import scipy.sparse as ssp


def plot_subunit(clusterdatum=dict, largesmall=str):
    nbrpairs = tree2tuplearr(clusterdatum['nbrtree'])
    print(nbrpairs)
    print('got ls ', largesmall)
    namespace = names.largeSubunit if (
        largesmall == 'l') else names.smallSubunit
    namesArr = names.largeSubunitArr if (
        largesmall == 'l') else names.smallSubunitArr

    dim = len(namesArr)

    substrate = np.zeros((dim, dim))

    for pair in nbrpairs:
        if (pair[0] not in namesArr or pair[1] not in namesArr):
            pass
        else:
            indexself = namesArr.index(pair[0])
            indexnbr = namesArr.index(pair[1])
            substrate[indexself, indexnbr] = 1

    fig, ax = plt.subplots()
    im = ax.imshow(substrate)
    ax.set_xticks(np.arange(len(namesArr)))
    ax.set_yticks(np.arange(len(namesArr)))
    ax.set_yticklabels(namesArr)
    ax.set_xticklabels(namesArr)


    plt.setp(ax.get_yticklabels(),  ha="right",
             rotation_mode="anchor", fontsize=6)

    ax.xaxis.tick_top()
    plt.setp(ax.get_xticklabels(), rotation=90, ha="left",
             rotation_mode="anchor", fontsize=6)

    ax.set_title("Small Subunit")

    fig.tight_layout() 
    plt.show()


for i in range (0, 10):
    datum = openjson('./../clusterdata/5T5H/' +
                    os.listdir('./../clusterdata/5T5H/')[i])
    plot_subunit(datum, 'l')
