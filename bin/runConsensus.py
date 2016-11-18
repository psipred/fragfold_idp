#!/usr/bin/env python
# -*- coding: utf-8 -*-
#

import pybrain
from train_nn import *
import pybrain
import sys
import pickle


"""
    5th step in the process. Takes Dynamine and FFIDP RMSD profile and builds
    consensus output
"""


def run_network(ffidp_fp, dm_fp, net_fp, out_fp):
    '''
    main script to run the neural network

    returns:
    np.array with predicted per-residue RMSD values
    '''

    fout = open(out_fp, 'a')
    net = pybrain.tools.customxml.networkreader.NetworkReader.readFrom(net_fp)
    inp = network_input(ffidp_fp, dm_fp)

    nn = []
    for res in inp:
        nn.append(net.activate(res))

    return np.array(nn)


def network_input(ffidp_fp, dm_fp):
    return net_inp


def read_inp(finp):
    '''
    read input data as a list of lists
    '''
    f = open(str(finp), 'r')
    lines = f.readlines()
    f.close()

    inp = []
    for line in lines:
        inp.append(map(float, line.split()))
    return inp


def SlidingWindow(inp_full, window_size=9):
#####
#####	remove log(length), N_dist and C_dist features from non-central residues
#####

    # assert window size is an odd number
    assert window_size % 2 == 1

    inp_sw = []
    inp = []
    central_res = []
    winLR = (window_size-1)/2

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
