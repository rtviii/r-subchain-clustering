from Bio import PDB
from get_clusters_of_neighbors import clustersets
import os
import numpy as np
import json


def extract_nbrtree_from_tuples(nbrtuples=list):
    nbrtree = {}
    # Create keys
    for seedpair in nbrtuples:
        if not(seedpair[0] in nbrtree.keys()):
            nbrtree[seedpair[0]] = []
        if not(seedpair[1] in nbrtree.keys()):
            nbrtree[seedpair[1]] = []

    for pair in nbrtuples:
        # if second chain of the two is not in the cluster of the first,
        # add it to the cluster (they are neigbor, after all)
        if pair[1] not in nbrtree[pair[0]]:
            nbrtree[pair[0]].append(pair[1])
        if pair[0] not in nbrtree[pair[1]]:
            nbrtree[pair[1]].append(pair[0])
    return nbrtree


def verify_no_rna(nbrtree=dict, allchains=list, rnas=list):
    for key in nbrtree.keys():
        for polymer in rnas:
            if polymer in nbrtree[key]:
                print("RNA {} in protein-only neighbor tree!".format(polymer))
                raise ValueError
        if key not in allchains:
            print("Key {} is not struct's chains: {}!".format(key, allchains))
            raise ValueError


def get_metadata(nbrclusters=list, allchains=list, rnas=list):
    # --------------------------------------------TALLYING-----------------------------------------------------
    singularChains = []
    mergedClusters = set([])
    # Merge everything
    for cluster in nbrclusters:
        mergedClusters.update(cluster)
    for chain in allchains:
        if chain not in mergedClusters and chain not in rnas:
            singularChains.append(chain)

    def countClustersChains(clusters=list):
        # returns the number of [chains, clusters]
        chaincount = 0
        clustercount = 0
        for cluster in clusters:
            clustercount += 1
            for chain in cluster:
                chaincount += 1
        return [chaincount, clustercount]

    totalcount = countClustersChains(nbrclusters)
    return {
        "n_clusters": totalcount[1],
        "n_clustered": totalcount[0],
        "n_singular": len(singularChains),
        "n_rnas": len(rnas),
        "n_total": len(allchains),
        "singular_chains": singularChains,
        "rnas": rnas,
        "allchains ": allchains
    }


def extract_clusters(pdbid=str, rnas=list, radius=int):
    filepath = os.path.realpath("./../.cif_models/{}.cif".format(pdbid))
    requestedstruct = pdbid
    requestedradius = radius
    rnas_to_exclude = rnas

    print("Opening {}".format(filepath))
    try:
        with open(filepath) as f:
            parser = PDB.FastMMCIFParser(QUIET=True)
            structure = parser.get_structure(requestedstruct, f)
    except:
        print("Failed to open {}".format(filepath))

    allchains = [chain.get_id() for chain in structure.get_chains()]
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
    nbridtuples = [(nbrtuple[0].get_id(), nbrtuple[1].get_id())
                   for nbrtuple in all_nbr_pairs]
    # so, the neighbor tree is needed to provide a basis for clusters
    nbrtree = extract_nbrtree_from_tuples(nbridtuples)
    verify_no_rna(nbrtree, allchains, rnas_to_exclude)

    # Extract clusters from the neighbortree
    nbrclusters = [list(c)for c in clustersets(nbrtree, rnas_to_exclude)]
    metadata = get_metadata(nbrclusters, allchains, rnas_to_exclude)

    return {
        "metadata": metadata,
        "clusters": nbrclusters
    }


# >>>>>Past tries, emailed Thomas Hamelryck, the author of NeighborSearch among other things.
# >>>> Still not sure if chain neighbors given a CHAIN could be gotten.
# >>>> center in search has to be np.shape = (3,)
# atomwise_K =        PDB.Selection.unfold_entities(subchain_dict['K'], 'A')
# np_atomwise_K =     np.array([atom.get_coord() for atom in list(atomwise_K)])
# struct_ns =         PDB.NeighborSearch(atomwise_struct, 10)
# print(np.shape(np_atomwise_K))                                   # (1451,3)
# chain_nbrs_single = struct_ns.search(np_atomwise_K[43], 10, "C") # works well, yields neighbors in 10 of 43'rd atom
# chain_nbrs = struct_ns.search(np_atomwise_K, 10, "C")            # Exepected a 3-dimensional Numpy array.
