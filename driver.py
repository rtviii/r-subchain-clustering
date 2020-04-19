from os import listdir
from os.path import isfile, join
from pathlib import Path
import os
import numpy as np
from Bio import PDB
import argparse
from clustering_scripts.save_clusters import save_clusters
from graph_algorithms.fiddling import test
import json

# DEFINE THE BATCH FOR NOW
with open('./../assets/kddtable.json', 'r') as infile:
    data = json.load(infile)
    batch = data.keys()

def openjson(fullpath=str):
    with open(fullpath, 'r') as infile:
        jsonprofile = json.load(infile)
        return jsonprofile



def cli():
    parser = argparse.ArgumentParser(description='Process pdb clustering')
    # parser.add_argument('-q',"--quantity",help='single or batch')
    parser.add_argument('-m','--mode', help='clustering or matrix operations')
    #[clust,graph]
    parser.add_argument('pdbid', help='PDB ID')
    parser.add_argument('radius', help='clustering radius', type=float)
    args = parser.parse_args()
    if (args.mode=='g'):
        test()
        return

    # if (args.mode=='clust'):
    save_clusters(args.pdbid, args.radius)



def cluster_all(ids, rad=float):
    print(batch)
    for key in batch:
        save_clusters(key, rad)


def openjson(fullpath=str):
    with open(fullpath, 'r') as infile:
        jsonprofile = json.load(infile)
        return jsonprofile





cli()

