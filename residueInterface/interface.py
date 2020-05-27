import matplotlib
matplotlib.use('TKAgg')
import Bio.PDB
import numpy
# from tkinter import *
import matplotlib.pylab as plt





pdb_code = "5JVG"
pdb_filename = "5JVG.pdb" #not the full cage!


def calc_residue_dist(residue_one, residue_two) :
    """Returns the C-alpha distance between two residues"""
    diff_vector  = residue_one["CA"].coord - residue_two["CA"].coord
    return numpy.sqrt(numpy.sum(diff_vector * diff_vector))

def calc_dist_matrix(chain_one, chain_two) :
    """Returns a matrix of C-alpha distances between two chains"""
    answer = numpy.zeros((len(chain_one), len(chain_two)), numpy.float)
    for row, residue_one in enumerate(chain_one) :
        for col, residue_two in enumerate(chain_two) :
            answer[row, col] = calc_residue_dist(residue_one, residue_two)
    return answer


structure = Bio.PDB.PDBParser(QUIET=True).get_structure(pdb_code, pdb_filename)


model = structure[0]
coord1 = [* model['H'].get_residues() ]
print(coord1[0]['CA'].get_vector())
# calc_residue_dist(model['H'], model['I'])

dist_matrix = calc_dist_matrix(model["H"], model["H"])
contact_map = dist_matrix < 6.0

# print ( "Minimum distance", numpy.min(dist_matrix) )
# print ( "Maximum distance", numpy.max(dist_matrix) )


