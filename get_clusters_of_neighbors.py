import os
import numpy as np


def clustersets(nbrtree, rnaIds=list):

    flatitems = []
    keys = nbrtree.keys()
    for key in keys:
        flatitems.append((str(key), [str(item) for item in nbrtree[key]]))

    clusters = []
    for kvpair in flatitems:
        newcluster = set(kvpair[1])
        # Make sure that .update() is passed a LIST. Strings get split into chars otherwise.
        newcluster.update([kvpair[0]])
        clusters.append(newcluster)

    def merge_sets(clusters=list):
        # Merge sets is a single pass of uniting non-disjoint sets in an array
        # Total disjoint is the recursive pass of mergesets
        print("\nMerging non-disjoint sets....")

        lead = 0
        counter = 1
        while lead < len(clusters) and len(clusters) > 1:
            # if lead is the only cluster then lead ===counter
            while counter < len(clusters):
                if counter == lead:
                    counter += 1
                    continue
                if not clusters[lead].isdisjoint(clusters[counter]):
                    clusters[lead].update(clusters[counter])
                    del clusters[counter]
                    if lead > counter:
                        lead -= 1
                    continue
                counter += 1
            lead = lead + 1
            counter = 0
        return clusters

    def total_disjoint(clusters=list):
        i = 0
        k = 0
        while i < len(clusters):
            while k < len(clusters):
                if k != i:
                    if clusters[k].isdisjoint(clusters[i]):
                        k = k + 1
                    else:
                        i = 0
                        k = 0
                        total_disjoint(merge_sets(clusters))
                else:
                    k = k+1
            k = 0
            i = i + 1
        print("Verifier has returned...")
        return clusters

    nbrclusters = total_disjoint(clusters)
    return nbrclusters
