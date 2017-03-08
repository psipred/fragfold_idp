import numpy as np
from itertools import combinations
import glob
from os.path import splitext
from manipulatePDB import ReadPDB
import os
import yaml
import argparse
from manipulatePDB import skip_comments, print_output

"""
    The script takes an ensemble of structures (PDB file) and performs sliding
    window superposition on all structures in an ensemble. It outputs
    a per-residue RMSD profile of the ensemble.
"""


def run_Sliding_Window(ens_fp, SW_size=10, to_file=False):

    # READ *ENS* FILES
    _ens = []
    if '*' in ens_fp:
        for inp_file in glob.glob(ens_fp):
            _ens.append(ReadPDB(inp_file))
    else:  # in CA mode
        _ens.append(ReadPDB(ens_fp))

    ens = np.array(_ens)
    PerRes_RMSD = sliding_window_RMSD(ens, sw=int(SW_size))

    if to_file is False:
        return PerRes_RMSD
    else:
        if '.' not in to_file:
            to_file = to_file + '.ffidp'
        np.savetxt(to_file, PerRes_RMSD, fmt='%.3f')


def sliding_window_RMSD(*args, **opts):
    '''
    Main Sliding_Window function

    args:

    opts:

    returns:
    np.array of per-residue RMSD values
    '''
    #print(args)
# **PRE-PROCESSING**
# prepare proteins np.array, and
# find out the length of the protein and number of conformations to analyse
    for arg in args:
        if type(arg) is not (np.ndarray or numpy.ndarray):
            raise TypeError('Numpy array required in PerResRMSD function!')
        if 'first_shape' in locals():
            if first_shape != arg.shape:
                raise ValueError('Shape mismatch in PerResRMSD')
        first_shape = arg.shape
    sw = 1
    if opts:
        if opts['sw']:
            sw = opts['sw']

# check input
    if len(args[0]) > 1:
        n_conf = len(args[0][0])
        prot1 = args[0][0]
        if len(prot1.shape) > 1:
            n_aa = prot1.shape[-2]
        else:
            n_aa = 1
        proteins = np.zeros((n_conf, n_aa, 3))
        for n in range(n_conf):
            proteins = args[0][n]
        del prot1
    elif len(args[0]) == 1:
        proteins = args[0][0]
        n_conf = proteins.shape[0]
        n_aa = proteins.shape[1]
    else:
        raise ValueError('No valid input provided in PerResRMSD')
##############

# **MAIN METHOD**
    SW_RMSD = {key: [] for key in range(n_aa)}
    tmp = np.zeros(sw)

    for aa in range(n_aa):
        begin_sw = int(aa-((sw-1)/2.0))
        if begin_sw < 0:
            begin_sw = 0
        end_sw = int(aa+((sw-1)/2.0))
        if end_sw > n_aa:
            end_sw = n_aa

        tmp = multi_RMSD(proteins[:, begin_sw:end_sw, :],
                         per_res=True,
                         max_iter=3,
                         thr_val=1.0)

        for t in range(len(tmp)):
            SW_RMSD[begin_sw+t].append(tmp[t])

    out = []
    for key in SW_RMSD:
        out.append(sum(SW_RMSD[key])/len(SW_RMSD[key]))

    return np.array(out)


def multi_RMSD(*args, **opts):
    '''
    calculates RMSD (total or per-residue) on a set of models for all
    pair-wise combinations

    input:
    a set of xyz arrays of coordinates (eg. prot1_xyz, prot2_xyz...) *OR*
    a multidimensional array with various conformations (eg. prot[conf,xyz]).
    NUMPY arrays required

    returns:
    a vector with per-residue RMSD values *OR* mean total RMSD
    '''

# *PRE-PROCESSING*
    per_res = False
    max_iter = False
    thr_val = False
    if opts:
        if opts['per_res']: per_res = True
        if opts['max_iter']: max_iter = opts['max_iter']
        if opts['thr_val']: thr_val = opts['thr_val']

# check input
    if len(args) > 1:
        n_conf = len(args)
        prot1 = args[0]
        if len(prot1.shape) > 1:
            n_aa = prot1.shape[-2]
        else:
            n_aa = 1
        proteins = np.zeros((n_conf,n_aa,3))

        for n in range(n_conf):
            prot = args[n] - Centroid(args[n])
            proteins[n,:,:] = prot
        del prot1
        del prot
    elif len(args) == 1:
        prot1 = args[0]
        n_conf = prot1.shape[0]
        n_aa = prot1.shape[1]
        proteins = np.zeros((n_conf,n_aa,3))
        # do a centroid flip for each protein
        for n in range(n_conf):
            prot = args[0][n] - Centroid(args[0][n])
            proteins[n,:,:] = prot
        del prot1
        del prot
    else:
        raise ValueError('No valid input provided in PerResRMSD')
#############

# calculate the list of all pairs
    combinations_list = list(combinations(range(n_conf), 2))
    n_combinations = len(combinations_list)

# calculate RMSD (tot OR per-res) for all pairs
    pairwise_RMSD = [None] * n_combinations
    for comb_no, comb in enumerate(combinations_list):
        pairwise_RMSD[comb_no] = RMSD_fit(proteins[comb[0],:,:],
                                          proteins[comb[1],:,:],
                                          per_res,
                                          thr_val,
                                          max_iter)

    pairwise_RMSD = np.array(pairwise_RMSD)**2

# calculate mean value accross all pairs
    mean_RMSD = np.zeros(pairwise_RMSD.shape[-1])
    for pair in pairwise_RMSD:
        mean_RMSD += pair

    if per_res:
        return np.sqrt(mean_RMSD / n_combinations)
    else:
        return np.sqrt(mean_RMSD[0] / n_combinations)


def Centroid(X):
    return sum(X)/len(X)


def RMSD_fit(P, Q, per_res=False, thr_val=1e-3, max_iter=20, mean_struc=False):
    """
    Varies the distance between P and Q, and optimizes rotation for each step
    until a minimum is found.

    ADAPTED FROM: https://github.com/charnley/rmsd
    """
# *PRE-PROCESSING*
    assert(P.shape == Q.shape)
###

    no_steps = 0
    step_size = P.max(0)
    threshold = step_size*thr_val
    rmsd_best = Kabsch(P, Q)
    P_prev = P
    while True:
        P = P_prev
        no_steps += 1
        for i in range(3):
            temp = np.zeros(3)
            temp[i] = step_size[i]
            rmsd_new = Kabsch(P+temp, Q)
            if rmsd_new < rmsd_best:
                rmsd_best = rmsd_new
                P_prev[:, i] += step_size[i]
            else:
                rmsd_new = Kabsch(P-temp, Q)
                if rmsd_new < rmsd_best:
                    rmsd_best = rmsd_new
                    P_prev[:, i] -= step_size[i]
                else:
                    step_size[i] /= 2
        if (step_size < threshold).all() or no_steps > max_iter:
            break

    if mean_struc:
        return RMSD_dumb(P, Q, True), np.array((P+Q)/2.0)
    else:
        if per_res:
            return Kabsch(P, Q, per_res=True)
        else:
            return rmsd_best


def Kabsch(arrA, arrB, output=False, per_res=False):
    '''
    Implementation of Kabsch algoritm for minimising the RMSD
    of 2 protein structures
    '''

    # covariance matrix
    C = np.dot(np.transpose(arrA), arrB)
    v, s, w = np.linalg.svd(C)

    # if the two proteins are reflections of each other
    d = (np.linalg.det(v) * np.linalg.det(w)) < 0.0
    if d:
        s[-1] = -s[-1]
        v[:, -1] = -v[:, -1]

    # rotation matrix
    U = np.dot(v, w)

    # rotate arrA
    arrA = np.dot(arrA, U)

    if output:
        return arrA, RMSD_dumb(arrA, arrB)
    else:
        if per_res:
            return RMSD_dumb(arrA, arrB, per_res)
        else:
            return RMSD_dumb(arrA, arrB)


def RMSD_dumb(arrA, arrB, per_res=False):
    '''
    naive RMSD calculation (no rotations)

    takes 2 arrays (arrA, arrB) and an argument if per-residue RMSD requested

    returns:
    list of RMSD values *OR*
    single RMSD value
    '''

# find number of atoms
    if len(arrA.shape) > 1:
        n_at = arrA.shape[-2]
        n_dim = arrA.shape[-1]
    else:
        n_at = 1
        n_dim = arrA.shape[0]

    rmsd = []
    for A, B in zip(arrA, arrB):
        rmsd.append([(A[i]-B[i])**2.0 for i in range(n_dim)])

    if per_res:
        return [np.sqrt(sum(r)) for r in rmsd]
    else:
        return np.sqrt(sum([sum(r) for r in zip(*rmsd)])/n_at)


script_path = os.path.dirname(os.path.realpath(__file__))
paths_path = script_path+"/../paths.yml"
if os.path.isfile(paths_path):
    paths_yaml = open(paths_path)
    paths = yaml.load(paths_yaml)
else:
    print("Unable to find paths.yml.n\nShould be in top dir for fragfold_idp")
    exit()

parser = argparse.ArgumentParser(description='Runs the SlidingWindow'
                                 'superposition code over an input NMR pdbfile')
parser.add_argument('--input_file', help="Input NMR PDB file", required=True)
parser.add_argument('--input_name', help="Results naming", required=True)
parser.add_argument('--outdir',
                    help="Where to put the output files",
                    default=script_path+"/../output/")
parser.add_argument('--window_size',
                    help="Where to put the output files",
                    default=10)
args = parser.parse_args()

# print(args.input_file)
# print(args.outdir+"/"+args.input_name+".pdb_ens")

sw_out = run_Sliding_Window(args.input_file,
                            SW_size=args.window_size)# value = run_Sliding_Window("/home/dbuchan/Code/fragfold_idp/output/a15a6b5e-9463-11e6-a62a-989096c13ee6_CLUSTER_001.pdb")


out_fp = args.outdir+"/"+args.input_name+".pdb_ens"
fasta_fp = args.indir+args.input_name+".fasta"

print_output(fasta_fp, sw_out, out_fp)
