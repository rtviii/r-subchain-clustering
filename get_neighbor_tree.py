from Bio import PDB
from get_clusters_of_neighbors import clustersets
import os
import numpy as np
import json

# this should be a view: all the calculations are to be extracted into a dedicated method


def get_neighbor_tree(filename, rnas, radius):
    filename = filename
    filepath = os.path.realpath("./../.cif_models/{}.cif".format(filename))

    requestedstruct = filename
    requestedradius = float(radius)
    rnas_to_exclude = rnas

    print("\n\n----------------------------------------------------------------------------------------")
    print("# Attempting to identify clusters of structure [{}] with the radius of [{}] Å       ".format(
        requestedstruct, radius))

    print("attempting to open {}".format(filepath))
    with open(filepath) as f:
        parser = PDB.FastMMCIFParser(QUIET=True)
        structure = parser.get_structure(requestedstruct, f)
        allchains = [chain.get_id() for chain in structure.get_chains()]

    # unfold to atoms
    atomwise_struct = PDB.Selection.unfold_entities(structure, "A")
    structwide_ns = PDB.NeighborSearch(atomwise_struct, bucket_size=5)
    all_nbr_pairs = structwide_ns.search_all(radius, "C")
    i = 0
    # Excluding RNAs from the protein neighbor pairs
    while i < len(all_nbr_pairs):
        if all_nbr_pairs[i][0].get_id() in rnas_to_exclude or all_nbr_pairs[i][1].get_id() in rnas_to_exclude:
            all_nbr_pairs.pop(i)
            pass
        else:
            i += 1

    # so, the neighbor tree is needed to provide a basis for clusters
    # only the chains, who have neighbors in common belong to the same cluster
    # Neverfcking mind.  T Hamelryck says this is fast enough.
    # So either rewrite in WAS or shut up and code on. Huh. Maybe i will.
    nbrtree = {}
    # if a key doesnt already exists in the tree, create it

    for nbrtuple in all_nbr_pairs:
        if not(nbrtuple[0].get_id() in nbrtree.keys()):
            nbrtree[nbrtuple[0].get_id()] = []
        if not(nbrtuple[1].get_id() in nbrtree.keys()):
            nbrtree[nbrtuple[1].get_id()] = []

    for nbrtuple in all_nbr_pairs:
        if nbrtuple[1].get_id() not in nbrtree[nbrtuple[0].get_id()]:
            nbrtree[nbrtuple[0].get_id()].append(nbrtuple[1].get_id())
        if nbrtuple[0].get_id() not in nbrtree[nbrtuple[1].get_id()]:
            nbrtree[nbrtuple[1].get_id()].append(nbrtuple[0].get_id())

    for key in nbrtree.keys():
        for polymer in rnas_to_exclude:
            if polymer in nbrtree[key]:
                print("RNA {} in protein-only neighbor tree!".format(polymer))
                raise ValueError
        if key not in allchains:
            print("Key {} is not struct's chains: {}!".format(key, allchains))
            raise ValueError

    # Set is not JSON Serializable. Parse to arrays for now, use pickle or yml later.
    nbr_clusters = [list(c)
                    for c in clustersets(nbrtree, rnas_to_exclude)]

    # --------------------------------------------TALLYING-----------------------------------------------------
    singularChains = []
    mergedClusters = set([])
    for cluster in nbr_clusters:
        mergedClusters.update(cluster)
    for chain in allchains:
        if chain not in mergedClusters and chain not in rnas_to_exclude:
            singularChains.append(chain)

    def totalchains(clusters=list):
        # returns the number of [chains, clusters]
        c = 0
        n = 0
        for cluster in clusters:
            n += 1
            for chain in cluster:
                c += 1
        return [c, n]
    total = totalchains(nbr_clusters)
    print("{} RNA polymers(excluded from clustering): {}".format(
        len(rnas_to_exclude), rnas_to_exclude))
    print("{} proteins in {} clusters: {}".format(
        total[0], total[1], nbr_clusters))
    print("{} singular proteins: {}".format(
        len(singularChains), singularChains))
    print("Nominal # of chains in the molecule: {}".format(len(allchains)))

    # --------------------------------------------------------------------------------------------------------

    return [total[0],   len(singularChains), total[1], {
            "clusters": nbr_clusters,
            "singular": singularChains
            }]


print(get_neighbor_tree('3j9m', ['A', 'B', 'u', 'AA'], 2.5))



# >>>>>Past tries, emailed Thomas Hamelryck, the author of NeighborSearch among other things.
# >>>> Still not sure if chain neighbors given a CHAIN could be gotten.
# >>>> center in search has to be np.shape = (3,)
# atomwise_K =        PDB.Selection.unfold_entities(subchain_dict['K'], 'A')

# np_atomwise_K =     np.array([atom.get_coord() for atom in list(atomwise_K)])
# struct_ns =         PDB.NeighborSearch(atomwise_struct, 10)
# print(np.shape(np_atomwise_K))                                   # (1451,3)

# chain_nbrs_single = struct_ns.search(np_atomwise_K[43], 10, "C") # works well, yields neighbors in 10 of 43'rd atom
# chain_nbrs = struct_ns.search(np_atomwise_K, 10, "C")            # Exepected a 3-dimensional Numpy array.
