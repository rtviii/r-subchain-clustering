from os import listdir
from os.path import isfile, join
from pathlib import Path
import os
import numpy as np
from Bio import PDB
import argparse
from clustering_scripts.extract_clusters import extract_clusters
import json

# DEFINE THE BATCH FOR NOW
with open('./../assets/kddtable.json', 'r') as infile:
    data = json.load(infile)
    batch = data.keys()


def openjson(fullpath=str):
    with open(fullpath, 'r') as infile:
        jsonprofile = json.loads(infile)
        return jsonprofile


def download_cluster_configuration(molprofpath=str, targetradius=int):
    with open(molprofpath, 'r') as f:
        profile = json.load(f)
        pdbid = profile['metadata']['pdbid']
        rnas = profile['metadata']['rnanames']

    subchainclusters = extract_clusters(pdbid, rnas, targetradius)
    filename = "./clusterdata/{}/{}_withR={}_|_{}+{}__in{}_|_.json".format(
        pdbid,
        pdbid,
        f'{targetradius:.3f}',
        subchainclusters['metadata']['n_clustered'],
        subchainclusters['metadata']['n_singular'],
        subchainclusters['metadata']['n_clusters'],
    )

    if not (os.path.isdir('./clusterdata/{}/'.format(pdbid))):
        Path('./clusterdata/{}/'.format(pdbid)
             ).mkdir(parents=True, exist_ok=True)
    with open(filename, 'w') as out:
        json.dump(subchainclusters, out)
    return subchainclusters


def reconstructClustersInNomenclature(clusterdatum=json, pdbid=str):

    molprofile = openjson('./struct_profiles/{}.json'.format(pdbid))

    def chainid_to_nomenclature(jsonprofile=json, chainid=str):
        for polymer in jsonprofile.polymers:
            if chainid == polymer.chainid:
                return polymer.nomenclature[0]
        return -1

    print("clusterdatum", clusterdatum)


def driver():

    parser = argparse.ArgumentParser(description='Process pdbid.')
    parser.add_argument('pdbid', metavar='id')
    args = parser.parse_args()
    pdbid = args.pdbid

    molecular_profile_path = './struct_profiles/{}.json'.format(pdbid)
    cldatum = download_cluster_configuration(molecular_profile_path, 3)
    reconstructClustersInNomenclature(cldatum, pdbid)


def main():
    step = 0.5  # angstrom
    for key in batch:
        molecular_profile_path = './struct_profiles/{}.json'.format(key)
        print("Working with ", molecular_profile_path)
        for radius in np.arange(1.5, 8, step):
            x = download_cluster_configuration(molecular_profile_path, radius)


driver()

# main()
