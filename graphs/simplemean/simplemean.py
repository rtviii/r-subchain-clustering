import json
from graphs.protspaceLU import proteinspace_LU 
from graphs.protspaceSU import proteinspace_SU 
from graphs.proteinspace import protein_space_initial 

from functools import reduce
import os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
matplotlib.use("TkAgg")

# TODO: more samples
# TODO: rerun the clusters with a lean on the smaller clustercount and larger clustersize: in the nbrhood of 4-6
# TODO: eukaryotes | bacteria | small su | large su

paths = {
    'bacteria': './../MedianClustersBacteria',
    'eukarya': './../MedianClustersEukarya',
    'total': './../MedianClustersTotal',
}

def openjsonfile(filename=str):

    if (filename):
        with open(filename) as infile:
            print("Infile: {}".format(infile))
            return json.load(infile)
    else:
        print("File {} not found. ".format(filename))
        raise FileNotFoundError

def open_clusters_file(pdbid=str):
    path_to_clusterfiles = './../MedianClustersTotal/'
    filepath = path_to_clusterfiles + pdbid.upper() + ".json"
    with open(filepath) as f:
        return json.load(f)


def save_as_matrix(matrix=np.array, name=str):
    np.savetxt(name+'_adjmat.csv', np.around(matrix, decimals=3), fmt='%0.4f')


def get_protein_index(proteinname=str, proteinspace=protein_space_initial):
    ischloroplast = proteinname[len(proteinname)-1] == 'c'
    if ischloroplast:
        print("Chloroplast! :{}".format(proteinname))
        proteinname = proteinname[0:len(proteinname)-1]
        print("Truncated: ", proteinname)

    if (proteinspace.get(proteinname, False)) != False:

        return proteinspace.get(proteinname)
    else:
        print("\n NOT FOUND IN THE PROTNAMESPACE:{}\n".format(proteinname))
        return False


def construct_adjacency_matrix(pdbid=str, proteinspace=protein_space_initial, saveflag=int):
    dim = len(proteinspace)
    substrate = np.zeros((dim, dim))
    clusters_dict = open_clusters_file(pdbid)
    for cluster in clusters_dict['clusters']:
        for protein in cluster:
            inner_nbrs = [*cluster]
            inner_nbrs.remove(protein)
            index = get_protein_index(protein, proteinspace)
            for nbr in inner_nbrs:
                nbr_ind = get_protein_index(nbr)
                print("got index ", nbr_ind)
                if (nbr_ind) != False:
                    substrate[index-1][nbr_ind-1] = 1
    if (saveflag > 0):
        save_as_matrix(substrate, pdbid)
    return substrate

def sum_over_matrices_in(proteinspace):
    adjamts = []
    for molecule in [x[0:len(x)-5] for x in os.listdir('./../MedianClustersTotal/')]:
        x = construct_adjacency_matrix(
            molecule, proteinspace, -1)
        adjamts.append(x)
    sigma = reduce(lambda x, y: np.add(x, y), adjamts)
    sigma = sigma.astype(float)
    sigma = sigma/len(adjamts)
    return sigma

substrate_total = sum_over_matrices_in(protein_space_initial)

colocated = []
protspace_indexed = list( protein_space_initial.keys() )
for i, _ in enumerate(protein_space_initial):
    for j, __ in enumerate(protein_space_initial):
       if substrate_total[i,j] >0.6:
           colocated.append(protspace_indexed[i])
           colocated.append(protspace_indexed[j])
colocated = list(dict.fromkeys(colocated))
print(colocated)

