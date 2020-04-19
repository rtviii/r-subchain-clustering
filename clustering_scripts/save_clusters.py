import json
from pathlib import Path
import numpy as np
import os
from clustering_scripts.extract_clusters import extract_clusters
from nomenclature.nomenclature_map import nom_map_from_profile

def openjson(fullpath=str):
    with open(fullpath, 'r') as infile:
        jsonprofile = json.load(infile)
        return jsonprofile



def save_clusters( pdbid=str, targetradius=float ):
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
