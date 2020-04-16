import json
import numpy as np
import networkx as ns
from proteinspace import protein_space_initial


lunitproteins = json.load(open('./resources/lunitproteins.json'))
sunitproteins = json.load(open('./resources/sunitproteins.json'))


# FILENAME is a path to file
# Returns a json object

def openjsonfile(filename=str):
    if (filename):
        with open(filename) as infile:
            print("Infile: {}".format(infile))
            return json.load(infile)
    else:
        print("File {} not found. ".format(filename))
        raise FileNotFoundError


def open_clusters_file(pdbid=str):
    path_to_clusterfiles = './cluster_reports/'
    filepath = path_to_clusterfiles + pdbid.upper() + ".json"
    with open(filepath) as f:
        return json.load(f)

def get_protein_index(pdbid=str):
    data = protein_space_initial
    if (data[pdbid]):
        return data[pdbid]
    else:
        print("\nProtein not found in the space.\n")
        return IndexError


def construct_adjacency_matrix(pdbid=str):
    # construct len(space) x len(space) matrix of 0s
    # 1s in all the proteins within a cluster
    dim = len(protein_space_initial)
    print(dim)
    substrate = np.zeros((dim, dim))
    clusters_dict = open_clusters_file(pdbid)
    for cluster in clusters_dict['clusters']:
        for protein in cluster:
            inner_nbrs = [*cluster]
            inner_nbrs.remove(protein)
            index = get_protein_index(protein)
            print("\nproteins's {} index in space is {}".format(protein, index))
            for nbr in inner_nbrs:
                nbr_ind = get_protein_index(nbr)
                print('{} nbr\'s ind is {}'.format(nbr, nbr_ind))
                print("filling 1 at {},{}".format(index, nbr_ind))
                substrate[index-1][nbr_ind-1] = 1

    print(substrate)
    save_as_matrix(substrate)

def save_as_matrix(adjmatrix=np.array):
    np.savetxt('test.csv', np.around(adjmatrix, decimals=0), fmt='%1d')


construct_adjacency_matrix("5gak")
