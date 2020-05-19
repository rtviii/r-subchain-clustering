from os import listdir
from os.path import isfile, join
from pathlib import Path
import os
import numpy as np
from Bio import PDB
import argparse
# import seaborn as sns
from clustering_scripts.save_clusters import save_clusters
# from graph_algorithms.generalplot import plot_covmat, plot_simplemean
import json


# DEFINE THE BATCH FOR NOW


def openjson(fullpath=str):
    with open(fullpath, 'r') as infile:
        jsonprofile = json.load(infile)
        return jsonprofile


def cli():
    parser = argparse.ArgumentParser(description='Process pdb clustering')
    parser.add_argument('-m', '--mode', help='clustering or matrix operations')
    parser.add_argument('-a', '-all', help='run batch')
    parser.add_argument('-tgt', '--targets', help='batch of targets(clustererfiles) to plot')
    parser.add_argument('-plt', '--plottype', help='covariance matrix/ cuthillmcgee/ adjacency/ laplacians/ ')
    parser.add_argument('-t','--thread', help='single thread to process clusters for one id')

    parser.add_argument('--pdbid', help='PDB ID')
    parser.add_argument('--radius', help='clustering radius', type=float)

    args = parser.parse_args()

    print("GOT ARGS", args)
    if (args.mode == 'g'):
        print("GOT TARGETS: ", args.targets)
        if (args.plottype == 'covmat'):
            plot_covmat(args.targets)
            return
        if (args.plottype == 'cmk'):
            plot_simplemean(args.targets)
            return
        print("No specified plottyped.") 
        
    if (args.thread):
        cluster_single_thread(args.thread)

    # run clusters on batch
    if (args.mode == 'a'):
        cluster_all()
        return

    save_clusters(args.pdbid, args.radius)


def cluster_all():
    batch = []
    with open('./../assets/kddtable.json', 'r') as infile:
        data = json.load(infile)
        batch = data.keys()
    keys = list(batch)
    defective = ['5X8T', '3J9M', '5XXB', '4V7E', ]
    custom =[ '5MYJ', '5T2A']
    vp2 = [ '5VP2']
    for m in defective:
        keys.remove(m)

    for key in vp2:
        rad = 1.3
        while rad < 2.5:
            save_clusters(key, rad)
            rad += 0.1


def openjson(fullpath=str):
    with open(fullpath, 'r') as infile:
        jsonprofile = json.load(infile)
        return jsonprofile



def cluster_single_thread(pdbid):

    for key in [pdbid]:
        rad = 1
        while rad < 7:
            save_clusters(key, rad)
            rad += 0.1

cli()
