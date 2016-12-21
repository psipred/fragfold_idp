import argparse
from scipy import stats
import os
import yaml

"""
   The script takes any 2 disorder profiles (PDB, FRAGFOLD, consensus, DynaMine)
   and outputs Spearman's rank correlation coefficient (Rs) between them
"""


def read_profile(fp):
    pr_rmsd = []
    with open(fp, 'r') as f:
        for line in f:
            value = line.split()
            if len(value) > 1:
                pr_rmsd.append(float(value[1]))
            else:
                pr_rmsd.append(float(value[0]))
    return pr_rmsd


def Rs(vecQ, vecT, pval=False):
    '''
    calculate Spearman's rank correlation coefficient (Rs) between 2 vectors

    vecQ, vecT: vectors (lists) of the same length with per-residue RMSD values
    pval (bool): display p-value

    returns:
    Rs value (float); *optionally: P-value (if pval==True)
    '''

    rs, P = stats.spearmanr(vecQ, vecT)

    if pval:
        return rs, P
    else:
        return rs


script_path = os.path.dirname(os.path.realpath(__file__))
paths_path = script_path+"/../paths.yml"
if os.path.isfile(paths_path):
    paths_yaml = open(paths_path)
    paths = yaml.load(paths_yaml)
else:
    print("Unable to find paths.yml.n\nShould be in top dir for fragfold_idp")
    exit()

parser = argparse.ArgumentParser(description='calculate Spearman\'s rank'
                                             'correlation coefficient (Rs)'
                                             'between 2 disorder profiles')
parser.add_argument('-i', '--input_query', dest='i1',
                    help="path to query disorder profile",
                    required=True)
parser.add_argument('-j', '--input_target', dest='i2',
                    help="path to target disorder profile",
                    required=True)

args = parser.parse_args()

profile1 = read_profile(args.i1)
profile2 = read_profile(args.i2)

print(Rs(profile1, profile2))
