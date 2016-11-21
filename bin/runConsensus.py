#!/usr/bin/env python
# -*- coding: utf-8 -*-
#

# import pybrain
import numpy as np
import subprocess
import tempfile
import os
import shutil
import re


"""
    5th step in the process. Takes Dynamine and FFIDP RMSD profile and builds
    consensus output
"""


def run_network(ffidp_fp, dm_fp, ss_fp, aln_fp, net_fp, out_fp):
    '''
    main script to run the neural network

    returns:
    np.array with predicted per-residue RMSD values
    '''

    fout = open(out_fp, 'a')
    net = pybrain.tools.customxml.networkreader.NetworkReader.readFrom(net_fp)
    inp = network_input(ffidp_fp, dm_fp, ss_fp, aln_fp)

    nn = []
    for res in inp:
        nn.append(net.activate(res))

    return np.array(nn)


def network_input(ffidp_fp, dm_fp, ss_fp, aln_fp, win_size=9):
    '''
    input is in the form:
    DM	logistic(FF-IDP)	PSIPRED(3)	AAcomposition(21)	log(length)	log(N_dist)	log(C_dist)
    '''

    # read DynaMine results
    dm = read_inp(dm_fp, columns=(1, 2))
    # read FRAGFOLD-IDP results
    ff = list(map(lambda x: sigmoid(x), read_inp(ffidp_fp)))
    # read secondary structure predictions
    # PSIPRED SS file (NOT SS2!)
    ss = read_inp(ss_fp,  columns=(3, 6))

    # sanity checks
    try:
        assert len(dm) == len(ff)
        assert len(ss) == len(ff)
    except AssertionError:
        print('input file lengths do not match')

    # generate amino acid composition statistics
    hhfilt="/Users/tomasz/projects/MicroProt/tools/hh-suite/build/bin/hhfilter"
    aa = aa_composition(aln_fp, hhfilter=hhfilt)

    # generate list of input features
    inp_features = []
    for res in range(len(ff)):
        _feat = [dm[res], ff[res], ss[res][0], ss[res][1], ss[res][2]]
        for i in range(21):
            _feat.append(aa[res][i])
        inp_features.append(_feat)

    # ADD global fratures
    # log(len)
    # inp_features.append(np.log10(len(ff)))
    # log(N_dist)
    # log(C_dist)

    print(inp_features[1])

    # perform sliding window on list of input features
    net_inp = SlidingWindow(inp_features, win_size)

    # list of lists
    return net_inp


def read_inp(finp, columns=False):
    '''
    read input data as a list of lists
    '''
    f = open(str(finp), 'r')
    lines = f.readlines()
    f.close()

    inp = []
    for line in lines:
        if not columns:
            inp.append(list(map(float, line.split())))
        else:
            inp.append(list(map(float, line.split()[columns[0]:columns[1]])))
    return np.squeeze(inp)


def aa_composition(msa_fp, hhfilter="$HHLIB/bin/hhfilter"):
    '''
    calculate AA composition statistics

    input:
    file path to file with MSA without headers (e.g. ffaln file)

    returns:
    2D numpy array of 21-element lists of residue frequencies
    '''
    # perform HHfilt on input alignment
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
    f = open(alnfilt_fp, 'r')
    lines = f.readlines()
    f.close()
    for line in lines:
        aln_arr.append(list(line.rstrip()))

    aln = np.asarray(aln_arr)

    Naln = AlnToNum(aln)
    N = float(aln.shape[1])
    freq = np.zeros((N, 21))

    for pos, col in enumerate(Naln.transpose()):
        for aa in range(21):
            freq[pos, aa] = list(col).count(aa)/float(len(col))

    # clean up temp files
    shutil.rmtree(tempdir)

    return freq


def SlidingWindow(inp_full, window_size=9):
#####
#####	remove log(length), N_dist and C_dist features from non-central residues
#####

    # assert window size is an odd number
    assert window_size % 2 == 1

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


def AlnToNum(aln):
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


def sigmoid(x):
    y = 1.0/(1.0 + np.exp(-1.0*x))
    return y

def scale_out(raw_out):
    '''
    takes numpy.array

    returns:
    logistic-normalized output
    '''
    mean = 0.986
    scaled_out = 1.0/(np.exp(-1.0*raw_out/(2+(2*mean))))
    return scaled_out
