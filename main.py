from os import listdir
from os.path import isfile, join
import os
from Bio import PDB
import get_neighbor_tree
import json


def main():
    thelist = list(molecule_rnas.keys())
    i = 13
    while i < len(thelist):
        pdbid = thelist[i]
        print("PDBID Is ", pdbid)
        r = 2.8
        while r < 10:
            master = {
                "molecule": pdbid,
                "rnas": molecule_rnas[pdbid],
                "radius": r,
                "clusters": {}
            }
            x = get_neighbor_tree(pdbid, molecule_rnas[pdbid], r)
            r += 0.05
            master["clusters"].update(x[3])
            filename = "./clusterdata/{}/{}_r{:.3f}__{}&{}__in{}.json".format(pdbid, pdbid, r, x[0], x[1], x[2])
            with open(filename, 'w') as out:
                json.dump(master, out)
        i += 1

    print("\n\n >>>>>>>>>>>> Finished. <<<<<<<<<<<<<")


molecule_rnas = {
    "5x8t": ["W", "A", "B"],
    "3j9m": ["A", "B", "u", "AA"],
    "5njt": ["A", "U", "V"],
    "5afi": ["a", "v", "w", "x", "y", "A", "B"],
    "5jvg": ["X", "Y"],
    "5myj": ["AA", "BA", "BB"],
    "5O60": ["A", "B", "2"],
    "5V7Q": ["A", "B"],
    "5nrg": ["X", "Y"],
    "4y4p": ["1A", "2A", "1B", "2B", "1a", "2a", "1v", "2v", "1w", "1y", "2w", "2y", "1x", "2x"],
    "4V9F": ["0", "9"],
    "4V6U": ["A1", "A2", "A0", "B1", "B3"],
    "6EK0": ["L5", "L7", "L8", "S2", "S6"],
    "5T2A": ["A", "B", "C", "D", "E", "2", "F", "G", "H"],
    "3J79": ["A", "B", "C"],
    "5GAK": ["1", "3", "4", "A", "B"],
    "4V7E": ["Ad", "Ae", "Af", "Aa", "Ac", "Ab"],
    "5T5H": ["A", "B", "C", "D", "E", "F", "G", "H"],
    # "5XXB": [],  # pdb unvailable
    "5XY3": ["1", "3", "4"],
    "3J7Z": ["B", "7", "A"],
    "5VP2": ["1A", "2A", "1B", "2B", "1a", "2a", "1v", "2v", "1w", "1y", "2w", "2y", "1x", "2x"],
    # "4UG0": []  # uniprot unavailable
}
main()
