from os import listdir
from os.path import isfile, join
from pathlib import Path
import os
import numpy as np
from Bio import PDB
from clustering_scripts.extract_clusters import extract_clusters
import json



#DEFINE THE BATCH FOR NOW
with open('./../assets/kddtable.json', 'r') as infile:
    data = json.load(infile)
    batch = data.keys()


def download_cluster_configuration(molprofpath=str, targetradius=int):
    with open(molprofpath, 'r') as f:
        profile = json.load(f)
        pdbid   = profile['metadata']['pdbid']
        rnas    = profile['metadata']['rnanames']

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


def main():
    step = 0.5  # angstrom
    for key in batch:
        molecular_profile_path = './struct_profiles/{}.json'.format(key)
        print("Working with ", molecular_profile_path)
        for radius in np.arange(1.5, 8, step):
            download_cluster_configuration(molecular_profile_path, radius)

main()
