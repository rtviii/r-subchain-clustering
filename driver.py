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


# def download_cluster_configuration(molprofpath=str, targetradius=int):
def download_cluster_configuration( pdbid, rnas, targetradius, nomenclatureMap=dict):

    pdbid = str.upper(pdbid)
    subchainclusters = extract_clusters(pdbid, rnas, targetradius, nomenclatureMap)
    filename = "./clusterdata/{}/{}_withR={}_|_{}+{}__in{}_|_.json".format(
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
    for key in batch:
        molecular_profile_path = './struct_profiles/{}.json'.format(key)
        print("Working with ", molecular_profile_path)
        for radius in np.arange(1.5, 8, step):
            x = download_cluster_configuration(molecular_profile_path, radius)


def runclusters_single():
    parser = argparse.ArgumentParser(
        description='Process pdbid to find clusters of.')
    parser.add_argument('pdbid', metavar='id')
    parser.add_argument('radius', metavar='r')
    args = parser.parse_args()
    pdbid = args.pdbid
    radius = float(args.radius)


    molecular_profile_path = './struct_profiles/{}.json'.format(str.upper(pdbid))
    profile                = openjson(molecular_profile_path)
    rnas                   = profile['metadata']['rnanames']

    nommap = nom_map_from_profile(profile)
    x = download_cluster_configuration(pdbid, rnas, radius, nommap)

def reconstructClustersInNomenclature(clusterdatum=json, pdbid=str):

    molprofile      = openjson('./struct_profiles/{}.json'.format(str.upper(pdbid)))
    nomenclaturemap = nom_map_from_profile(molprofile)

    print(molprofile)
    print('got molprofule')
    print('\n\n')
    print(clusterdatum)

    # def chainid_to_nomenclature(jsonprofile=json, chainid=str):
    #     for polymer in jsonprofile.polymers:
    #         if chainid == polymer.chainid:
    #             return polymer.nomenclature[0]
    #     return -1

    # print("clusterdatum", clusterdatum)



runclusters_single()


# main()
