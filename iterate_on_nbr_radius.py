from pathlib import Path
from get_neighbor_tree import get_neighbor_tree
import os
import json


def save_clusters_for_radius_instance(pdbid=str, rnas=list, radius=int):
    r = radius
    master = {
        "pdbid": pdbid,
        "rnas": rnas,
        "radius": r,
        "clusters": {}
    }
    x = get_neighbor_tree(pdbid, rnas, r)
    master["clusters"].update(x[3])
    filename = "./clusterdata/{}/{}_r{}__{}&{}__in{}.json".format(
        pdbid, pdbid, round(r, 3), x[0], x[1], x[2])

    if not (os.path.isdir('./clusterdata/{}/'.format(pdbid))):
        Path('./clusterdata/{}/'.format(pdbid)
             ).mkdir(parents=True, exist_ok=True)
    with open(filename, 'w') as out:
        json.dump(master, out)


# DO YOUR UNI STUFF FIRST.

save_clusters_for_radius_instance('3j9m', ['A', 'u', 'B', 'AA'], 2.35)
