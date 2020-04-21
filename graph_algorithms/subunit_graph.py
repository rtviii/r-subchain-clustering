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


def ajdmat(nbrpairs, namespace):
    dim = len(namespace)
    substrate = np.zeros((dim, dim))

    for pair in nbrpairs:
        if (pair[0] not in namespace or pair[1] not in namespace):
            pass
        else:
            print(pair, "\n")
            indexself = namespace.index(pair[0])
            indexnbr = namespace.index(pair[1])
            substrate[indexself, indexnbr] = 1
    return substrate


def plt_subunit_cov(filename=str, clusterdatum=dict, largesmall=str, covariance=False):
    nbrpairs = tree2tuplearr(clusterdatum['nbrtree'])
    namesArr = names.largeSubunitArr if (
        largesmall == 'l') else names.smallSubunitArr

    substrate = ajdmat(nbrpairs, namesArr)

    fig, ax = plt.subplots()
    if covariance:
        substrate = np.cov(substrate)
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
    plt.suptitle(filename, fontsize=6)
    fig.tight_layout()
    plt.show()


def pltall_single(filename=str, clusterdatum=dict, largesmall=str, covariance=False):

    nbrpairs = tree2tuplearr(clusterdatum['nbrtree'])
    namesArr = names.largeSubunitArr if (
        largesmall == 'l') else names.smallSubunitArr

    substrate = ajdmat(nbrpairs, namesArr)

    fig, (ax1, ax2, ax3) = plt.subplots(1, 3)

    if covariance:
        substrate = np.cov(substrate)
    ax1.imshow(substrate)
    ax2.imshow(np.cov(substrate))
    ax3.imshow(get_cuthill_mckee_ndarray(substrate, namesArr)[0])
    ax1.set_xticks(np.arange(len(namesArr)))
    ax1.set_yticks(np.arange(len(namesArr)))
    ax1.set_yticklabels(namesArr)
    ax1.set_xticklabels(namesArr)

    plt.setp(ax1.get_yticklabels(),  ha="right",
             rotation_mode="anchor", fontsize=6)

    ax1.xaxis.tick_top()
    plt.setp(ax1.get_xticklabels(), rotation=90, ha="left",
             rotation_mode="anchor", fontsize=6)

    ax1.set_title("Small Subunit")
    plt.suptitle(filename, fontsize=6)
    fig.tight_layout()
    plt.show()


def plt_batch(directory, namespace):
    print("NAMSAPCE", namespace)
    batch = os.listdir(directory)
    nbrtrees = [tree2tuplearr(
        openjson(directory + x)['nbrtree']) for x in batch]

    adjmats = [ajdmat(tree, namespace) for tree in nbrtrees ]
    reduced = functools.reduce(lambda x, y: np.add(x, y), adjmats)
    reduced = np.divide(reduced, len(adjmats))

    fig, ax = plt.subplots()
    im = ax.imshow(np.cov( reduced ))
    ax.set_xticks(np.arange(len(namespace)))
    ax.set_yticks(np.arange(len(namespace)))
    ax.set_yticklabels(namespace)
    ax.set_xticklabels(namespace)

    plt.setp(ax.get_yticklabels(),  ha="right",
             rotation_mode="anchor", fontsize=6)

    ax.xaxis.tick_top()
    plt.setp(ax.get_xticklabels(), rotation=90, ha="left",
             rotation_mode="anchor", fontsize=6)

    # ax.set_title('Directory: ', directory)
    # plt.suptitle(filename, fontsize=6)
    # fig.tight_layout()
    sns.heatmap(np.cov(reduced))
    plt.show()

# for i in range (0, 10):
filename = os.listdir('./../clusterdata/5myj/')[40]
datum = openjson('./../clusterdata/5MYJ/' + filename)

# plot_subunit(filename, datum, 's', covariance=True)
# pltall_single(filename, datum, 'l', covariance=True)
plt_batch('./../clusterdata/targetgroups/all/', names.largeSubunitArr)
