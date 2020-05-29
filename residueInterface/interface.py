import numpy as np
from Bio import PDB
import pandas as pd
import argparse


# KMEANS on... what? residues?
# figure out how to pull prot seq, homology data for proteins
# ambient saturation field to go with


def open_struct(pdbid, filepath):
    """Returns an open structure"""
    return PDB.FastMMCIFParser(QUIET=True).get_structure(pdbid, filepath)


def residue_distance(residue_one, residue_two):
    """Returns the C-alpha distance between two residues"""

    # Taking the average of each atom's position in each of the two residues
    avg1 = np.average([atom.get_coord() for atom in residue_one])
    avg2 = np.average([atom.get_coord() for atom in residue_two])
    diff_vector = avg1 - avg2
    return np.sqrt(diff_vector * diff_vector)


def chain_distance_matrix(c1, c2):
    """Returns a matrix of alpha-carbon distances between two chains"""
    answer = np.zeros((len(c1), len(c2)), np.float)
    for i, residue_one in enumerate(c1):
        for j, residue_two in enumerate(c2):
            answer[i, j] = residue_distance(residue_one, residue_two)
    return answer


# pdbid = "3j7z"
# path = pdbid + ".cif"
# structure = open_struct(pdbid, path)

# # STRUCTURE --> MODEL -->  CHAINS  --> RESIDUES --> ATOMS
# model = structure[0]

# parser = argparse.ArgumentParser()
# parser.add_argument('chain1')
# parser.add_argument('chain2')

# args = parser.parse_args()
# c1 = model[args.chain1]
# c2 = model[args.chain2]

# distmat = chain_distance_matrix(c1, c2)
# print(np.shape(distmat))
# print("Chain {} has {} residues".format(c1.get_id(), len([* c1])))
# print("Chain {} has {} residues".format(c2.get_id(), len([* c2])))
