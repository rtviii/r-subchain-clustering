from os import listdir
from os.path import isfile, join
from pathlib import Path
import os
import numpy as np
from Bio import PDB
import argparse
from clustering_scripts.extract_clusters import extract_clusters
from nomenclature.nomenclature_map import nom_map_from_profile
import json

# DEFINE THE BATCH FOR NOW
with open('./../assets/kddtable.json', 'r') as infile:
    data = json.load(infile)
    batch = data.keys()


def openjson(fullpath=str):
    with open(fullpath, 'r') as infile:
        jsonprofile = json.load(infile)
        return jsonprofile


def save_clusters_of( pdbid=str, targetradius=float ):
    pdbid                  = str.upper(pdbid)
    #
    molecular_profile_path = './struct_profiles/{}.json'.format(str.upper(pdbid))
    profile                = openjson(molecular_profile_path)
    rnas                   = profile['metadata']['rnanames']
    nommap                 = nom_map_from_profile(profile)
    subchainclusters       = extract_clusters(pdbid, rnas, targetradius, nommap)
    filename               = "./clusterdata/{}/{}_withR={}_|_{}+{}__in{}_|_.json".format(
        pdbid,
        pdbid,
        f'{ targetradius :.3f}',
        subchainclusters['metadata']['n_clustered'],
        subchainclusters['metadata']['n_singular'],
        subchainclusters['metadata']['n_clusters'],
    )
    if not (os.path.isdir('./clusterdata/{}/'.format(str.upper(pdbid)))):
        Path('./clusterdata/{}/'.format(pdbid)
             ).mkdir(parents=True, exist_ok=True)
    with open(filename, 'w') as out:
        json.dump(subchainclusters, out)
        print('Saved successfully at \t [{}]'.format(filename))
    return subchainclusters



def main():
    step = 0.5  # angstrom
    temprad = 3.1

    print(batch)
    for key in batch:
        save_clusters_of(key, temprad)
        
        
def cli_get_clusters():
    parser = argparse.ArgumentParser(
        description='Process pdbid to find clusters of.')
    parser.add_argument('pdbid', metavar='id')
    parser.add_argument('radius', metavar='r')
    args = parser.parse_args()
    pdbid = args.pdbid
    radius = float(args.radius)
    save_clusters_of(pdbid, radius)

# runclusters_single()


# main()
cli_get_clusters()
