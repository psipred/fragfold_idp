    #!/usr/bin/env python
# -*- coding: utf-8 -*-
#

import numpy as np
import pickle
import subprocess
import tempfile
import os
import shutil
import re
import yaml
import argparse
import glob

"""
    5th step in the process. Takes DynaMine and FFIDP RMSD profile and builds
    consensus output
"""


def run_network(ffidp_fp, dm_fp, ss_fp, aln_fp, net_fp, out_fp=None, scale=False):
    '''main script to run the neural network

    Parameters
    ----------
    ffidp_fp
        FRAGFOLD-IDP output file path
    dm_fp
        DynaMine output file path
    ss_fp
        PSIPRED SS2 file path
    aln_fp
        Multiple sequence alignment file path
    net_fp
        Neural network parameters file path (from data/)
    out_fp
        Output neural network results file path (default: None)
    scale: bool
        Should output network values be scaled for easier comparisons
        to other methods?

    Returns
    -------
    numpy.array
        predicted per-residue RMSD values
    '''

    net_fh = open(net_fp, 'r')
    net = pickle.load(net_fh)
    net_fh.close()

    inp = network_input(ffidp_fp, dm_fp, ss_fp, aln_fp)

    nn = []
    for res in inp:
        nn.append(net.activate(res)[0])
    nn = np.array(nn)

    if out_fp:
        np.savetxt(out_fp, nn, fmt='%.3f\n')

    if scale:
        nn = scale_out(nn)

    return nn


def network_input(ffidp_fp, dm_fp, ss_fp, aln_fp, win_size=9):
    '''Generate input for the neural network

    Parameters
    ----------
    ffidp_fp
        FRAGFOLD-IDP output file path
    dm_fp
        DynaMine output file path
    ss_fp
        PSIPRED SS2 file path
    aln_fp
        Multiple sequence alignment file path
    win_size: int
        Window size for sliding window size for neural network (default: 9)

    Returns
    -------
    list of lists
        list of per-residue features (246 features per residue if window size=9)

    Notes
    -----
    input is in the form:
    DM	sigmoid(FF-IDP)	PSIPRED(3)	AAcomposition(21)	log(length)	log(N_dist)	log(C_dist)
    '''

    # read DynaMine results
    dm = read_inp(dm_fp, columns=(1, 2))
    # read FRAGFOLD-IDP results and apply sigmoid transformation
    ff = list(map(lambda x: sigmoid(x, x0=1.432), read_inp(ffidp_fp)))
    # read secondary structure predictions
    # PSIPRED SS2 file
    ss = read_inp(ss_fp,  columns=(3, 6))

    # sanity checks
    try:
        assert len(dm) == len(ff)
        assert len(ss) == len(ff)
    except AssertionError:
        print('input file lengths do not match')

    # generate amino acid composition statistics
    aa = aa_composition(aln_fp)

    # generate list of input features
    inp_features = []
    for res, res_val in enumerate(ff):
        # ADD window features
        # DynaMine, FFIDP, Secondary structure (3) features
        _feat = [dm[res], ff[res], ss[res][0], ss[res][1], ss[res][2]]

        # ADD window features
        # AA composition (21)
        for i in range(21):
            _feat.append(aa[res][i])

        # ADD global fratures
        # log(len)
        n_res = len(ff)
        _feat.append(np.log10(n_res))
        # log(N_dist)
        _feat.append(np.log10(res+1))
        # log(C_dist)
        _feat.append(np.log10(n_res-res))

        inp_features.append(_feat)

    # perform sliding window on list of input features
    net_inp = SlidingWindow(inp_features, win_size)

    # list of lists
    return net_inp


def read_inp(finp, columns=False):
    '''read input data as a list of lists

    Parameters
    ----------
    finp: str
        input file path
    columns: tuple
        select input columns range to be read (default: False; all columns)

    Returns
    -------
    list of lists
        list of per-residue input features
    '''
    f = open(str(finp), 'r')
    lines = f.readlines()
    f.close()

    inp = []
    for line in lines:
        # ignore commented (#) and blank lines
        if line.startswith('#') or not line.strip():
            continue
        if not columns:
            inp.append(list(map(float, line.split())))
        else:
            inp.append(list(map(float, line.split()[columns[0]:columns[1]])))
    return np.squeeze(inp)


def aa_composition(msa_fp, hhfilter=None):
    '''calculate AA composition statistics

    Parameters
    ----------
    msa_fp: str
        file path to file with MSA without headers
        (e.g. `msa` file from example folder)
    hhfilter: str
        optional parameter to filter out similar sequences using HHfilter
        (default: None)

    Returns
    -------
    numpy.array
        2D numpy array of 21-element lists of residue frequencies
    '''

    if hhfilter:
    # perform HHfilt on input alignment
        hhfilter="$HHLIB/bin/hhfilter"
        tempdir = tempfile.mkdtemp(dir=os.getcwd())
        a3mfilt_fp = tempdir + '/a3mfilt'
        alnfilt_fp = tempdir + '/alnfilt'
        aln_fp = tempdir + '/aln'

        msa2fasta(msa_fp, aln_fp)

        hhfilter_cmd = '%s -i %s -id 90 -qid 30 -o %s' % (hhfilter,
                                                          aln_fp,
                                                          a3mfilt_fp)
        subprocess.call(hhfilter_cmd.split())

        grep = subprocess.Popen(' '.join(['grep -v "^>"', a3mfilt_fp,
                                "| sed '"'s/[a-z]//g'"' >", alnfilt_fp]),
                                shell=True, stdout=subprocess.PIPE)
        grep.communicate()

    aln_arr = []

    if hhfilter:
        f = open(alnfilt_fp, 'r')
    else:
        f = open(msa_fp, 'r')
    lines = f.readlines()
    f.close()
    for line in lines:
        aln_arr.append(list(line.rstrip()))

    aln = np.asarray(aln_arr)

    Naln = aln_to_num(aln)
    N = int(aln.shape[1])
    freq = np.zeros((N, 21))

    for pos, col in enumerate(Naln.transpose()):
        for aa in range(21):
            freq[pos, aa] = list(col).count(aa)/float(len(col))

    # clean up temp files
    if hhfilter:
        shutil.rmtree(tempdir)

    return freq


def SlidingWindow(inp_full, window_size=9):

    # assert window size is an odd number
    assert window_size % 2 == 1
    # check if `inp_full` is a list. If not, then make it a list
    if not type(inp_full) is list:
        inp_full = inp_full.tolist()

    inp_sw = []
    inp = []
    central_res = []
    winLR = int((window_size-1)/2)

    ## add missing residues feature
    # if missing: feature = 1
    # else: feature = 0
    for i in range(len(inp_full)):
    # chop out central residue features from the feature list and save them separately
        inp.append(inp_full[i][:-3])
        central_res.append(inp_full[i][-3:])
        inp[i].append(0)

    # num features - 3 central residue features [length of protein, N-term & C-term distances]
    n_feat = len(inp[0])-1

    # add blank features BEFORE and AFTER the list
    ## before
    inp.reverse()
    for i in range(winLR):
        empty = [0]*n_feat
        empty.append(1)
        inp.append(empty)
    inp.reverse()
    ## after
    for i in range(winLR):
        empty = [0]*n_feat
        empty.append(1)
        inp.append(empty)

    for res in range(len(inp)-window_size+1):
        tmp_inp = [i for sublist in inp[res:res+window_size] for i in sublist]
        # add central residue features
        tmp_inp = tmp_inp + central_res[res]
        inp_sw.append(tmp_inp)

    return inp_sw


def aln_to_num(aln):
    N_aln = np.zeros(aln.shape)
    for i, col in enumerate(aln.transpose()):
        N_aln[:, i] = AA_alphabet(col)
    return N_aln


def AA_alphabet(seq):
    '''
    input: list or string with 1-letter sequence
    output: list with numbered amino-acid code
    '''

    if type(seq) is list:
        seq = str(' '.join(seq))
    elif type(seq) is str and ',' in seq:
        seq = seq.replace(',', ' ')

    d = {'-': 0, 'C': 1, 'D': 2, 'S': 3, 'Q': 4, 'K': 5,
         'I': 6, 'P': 7, 'T': 8, 'F': 9, 'N': 10,
         'G': 11, 'H': 12, 'L': 13, 'R': 14, 'W': 15,
         'A': 16, 'V': 17, 'E': 18, 'Y': 19, 'M': 20,
         'X': 0, 'B': 10, 'U': 3}

    abc = np.zeros(len(seq))
    for i in range(len(seq)):
        abc[i] = d[seq[i]]
    return abc


def msa2fasta(in_fp, out_fp):
    f = open(in_fp, 'r')
    lines = f.readlines()
    f.close()
    out = open(out_fp, 'w')
    for line in lines:
        out.write(re.sub('^', '>\n', line))
    out.close()


def sigmoid(x, x0=0.0):
    y = 1.0/(1.0 + np.exp(-1.0*(x-x0)))
    return y

def scale_out(raw_out, mean=0.986):
    '''Scale neural network result to get results comparable with FRAGFOLD-IDP

    Parameters
    ----------
    raw_out: numpy.array
    mean: float
        mean value to scale to

    Returns
    -------
    numpy.array
        logistic-normalized output
    '''

    scaled_out = 1.0/(np.exp(-1.0*raw_out/(2+(2*mean))))
    return scaled_out

script_path = os.path.dirname(os.path.realpath(__file__))
paths_path = script_path+"/../paths.yml"
if os.path.isfile(paths_path):
    paths_yaml = open(paths_path)
    paths = yaml.load(paths_yaml)
else:
    print("Unable to find paths.yml.n\nShould be in top dir for fragfold_idp")
    exit()

parser = argparse.ArgumentParser(description='Builds a Consensus disorder ' \
                                 'profile from the FFIDP and Dynamine outputs')
parser.add_argument('--input_name', help="job name for the files generated in "
                                         "in the previous steps, usually "
                                         "a uuid", required=True)
parser.add_argument('--indir',
                    help="Where to find the .ens file. Defaults to the "
                         "runSeqAnalysis.py output directory",
                    default=script_path+"/../output/")
parser.add_argument('--outdir',
                    help="Where to put the output files",
                    default=script_path+"/../output/")
parser.add_argument('--ffidp_path',
                    help="The path for the FFIDP profile output",
                    default=script_path+"/../output/")
parser.add_argument('--dynamine_path',
                    help="The path for the FFIDP profile output",
                    default=script_path+"/../output/")
parser.add_argument('--psipred_path',
                    help="The path for the FFIDP profile output",
                    default=script_path+"/../output/")
parser.add_argument('--alignment_path',
                    help="The path for the FFIDP profile output",
                    default=script_path+"/../output/")

args = parser.parse_args()
# print(args.ffidp_path+args.input_name+".ffidp")
# print(args.ffidp_path+args.input_name+"/Dynamine_b_*")
# print(args.ffidp_path+args.input_name+".ss2")
# print(args.ffidp_path+args.input_name+".aln")
dynaResults = glob.glob(args.ffidp_path+args.input_name+"/Dynamine_b_*")[0]
dynamineOutput = ''
for dynaFile in glob.glob(dynaResults+"/*"):
    if ".pred" in dynaFile:
        dynamineOutput = dynaFile
# print(dynamineOutput)
# exit()
network_results = run_network(args.ffidp_path+args.input_name+".ffidp",
                              dynamineOutput,
                              args.ffidp_path+args.input_name+".ss2",
                              args.ffidp_path+args.input_name+".aln",
                              script_path+"/../data/network_params_example",
                              out_fp=args.outdir+args.input_name+".consensus"
                              )
