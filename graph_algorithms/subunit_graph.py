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
import argparse
import scipy.sparse as ssp


matplotlib.style.use('classic')




def ajdmat(nbrpairs, namespace):
    dim = len(namespace)
    substrate = np.zeros((dim, dim))

    for pair in nbrpairs:
        if (pair[0] not in namespace or pair[1] not in namespace):
            pass
        else:
            indexself = namespace.index(pair[0])
            indexnbr = namespace.index(pair[1])
            substrate[indexself, indexnbr] = 1
    return substrate


# Plot adjacency matrix or covarinace matrix for the instance of a radius
# by SUBUNIT : largesmall = ['l', 's']
def plt_subunit(filename=str, clusterdatum=dict, largesmall=str, covariance=False):
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


# Plot adjacency, covariance, reverse cuthill-mckee decomp for an instance of a radius
# by SUBUNIT : largesmall = ['l', 's']
def pltall_single(filename=str, clusterdatum=dict, largesmall=str):

    nbrpairs = tree2tuplearr(clusterdatum['nbrtree'])
    namesArr = names.largeSubunitArr if (
        largesmall == 'l') else names.smallSubunitArr

    substrate = ajdmat(nbrpairs, namesArr)
    covmat = np.cov(substrate)
    ckm, perm = get_cuthill_mckee_ndarray(substrate, namesArr)
    fig, (ax1, ax2, ax3) = plt.subplots(1, 3)

    ax1.imshow(substrate)
    ax2.imshow(covmat)
    ax3.imshow(ckm)

    ax1.set_xticks(np.arange(len(namesArr)))
    ax1.set_yticks(np.arange(len(namesArr)))
    ax1.set_yticklabels(namesArr)
    ax1.set_xticklabels(namesArr)

    ax2.set_xticks(np.arange(len(namesArr)))
    ax2.set_yticks(np.arange(len(namesArr)))
    ax2.set_yticklabels(namesArr)
    ax2.set_xticklabels(namesArr)

    ax3.set_xticks(np.arange(len(perm)))
    ax3.set_yticks(np.arange(len(perm)))
    ax3.set_yticklabels(perm)
    ax3.set_xticklabels(perm)

    plt.setp(ax1.get_yticklabels(),  ha="right",
             rotation_mode="anchor", fontsize=6)
    plt.setp(ax2.get_yticklabels(),  ha="right",
             rotation_mode="anchor", fontsize=6)
    plt.setp(ax3.get_yticklabels(),  ha="right",
             rotation_mode="anchor", fontsize=6)

    ax1.xaxis.tick_top()
    ax2.xaxis.tick_top()
    ax3.xaxis.tick_top()
    plt.setp(ax1.get_xticklabels(), rotation=90, ha="left",
             rotation_mode="anchor", fontsize=5)
    plt.setp(ax2.get_xticklabels(), rotation=90, ha="left",
             rotation_mode="anchor", fontsize=5)
    plt.setp(ax3.get_xticklabels(), rotation=90, ha="left",
             rotation_mode="anchor", fontsize=5)

    # ax1.set_title("Small Subunit")
    plt.suptitle(filename, fontsize=6)
    fig.tight_layout()
    plt.show()


# plot covariance matrix of adjacency averages inside a given namespace
def plt_batch(directory, namespace):
    batch = os.listdir(directory)
    nbrtrees = [tree2tuplearr(openjson(directory + x)['nbrtree'])
                for x in batch]

    adjmats = [ajdmat(tree, namespace) for tree in nbrtrees]
    reduced = functools.reduce(lambda x, y: np.add(x, y), adjmats)
    reduced = np.divide(reduced, len(adjmats))

    fig, ax = plt.subplots()
    im = ax.imshow(np.cov(reduced))
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
    fig.tight_layout()
    plt.show()


def plt_inner_namespace(pathtodatum=str, covariance=False, pdbid=str):
    clusterdatum = openjson(pathtodatum)

    radius = round(clusterdatum['metadata']['radius'], 2)
    nbrpairs = tree2tuplearr(clusterdatum['nbrtree'])
    namesArr = clusterdatum['metadata']['allchains']
    substrate = ajdmat(nbrpairs, namesArr)


    fig, ax = plt.subplots()
    fig.patch.set_facecolor('xkcd:black')
    if covariance:
        substrate = np.cov(substrate)
    im = ax.imshow(substrate)
    ax.set_xticks(np.arange(len(namesArr)))
    ax.set_yticks(np.arange(len(namesArr)))
    ax.set_yticklabels(namesArr)
    ax.set_xticklabels(namesArr)

    plt.setp(ax.get_yticklabels(),  ha="right",
             rotation_mode="anchor", fontsize=8, c='white')

    ax.xaxis.tick_top()
    plt.setp(ax.get_xticklabels(), rotation=90, ha="left",
             rotation_mode="anchor", fontsize=8, c='white')

    ax.set_title( '{} | Neighbor radius : {} Å'.format(pdbid, radius), color='white')
    # plt.suptitle(filename, fontsize=4)
    fig.tight_layout()
    filename=' {}_r{}_covmat.png'.format(pdbid,radius)
    plt.savefig(filename, bbox_inches='tight', facecolor=fig.get_facecolor())
    # plt.show()


clustpath = './../clusterdata/'
pdbid = '3J7Z'
count = 20
# filename = os.listdir('./../clusterdata/{}/'.format(pdbid))[count]
# plt_inner_namespace(clustpath + '{}/'.format(pdbid) + filename, covariance=False, pdbid=pdbid)


filenum = len(os.listdir('./../clusterdata/{}/'.format(pdbid)))
for i in range(filenum):
    filename = os.listdir('./../clusterdata/{}/'.format(pdbid))[i]
    path = clustpath + '{}/'.format(pdbid) + filename
    plt_inner_namespace(path, True, pdbid)
