# fragfold_idp

## Introduction

FRAGFOLD-IDP is a complex multi stages processs with a large number of dependencies. These scripts are designed to help automate installation and running FF-IDP.

### Installation

We provide an ansible script for installation to any standard linux distro, see below. You must have ansible and python2 available. See https://github.com/ansible/ansible for further details.

### Runtime

1. Run basic sequences analysis
2. Build an ensemble of FRAGFOLD models
3. Cluster the models and produce a RMSD profile
4. Run DynaMine
5. Build Consensus output
6. Generate evaluation statistics

## How to Install

1. Checkout ansible https://github.com/ansible/ansible
2. OPTIONAL: switch to python2 virtualenv
3. install pyyaml

`pip install pyyaml`

4. Get a Dynamine API key from http://dynamine.ibsquare.be/download/
5. Edit the paths.yml to reflect where you want to install the various
   components or where they have already been installed.

6. Run setup script for ansible

`source ~/ansible/hacking/env-setup`

7. Change to the ansible script directory

`cd ~/fragfold_idp/ansible`

8. Run FFIDP installation script

`ansible-playbook -i hosts install.yml -f 1`

## How to run

1. Run the Sequence analysis script over each of your fasta files. This will
also produce the Fragfold input files. By default the results go in the
fragfold_idp/output dir. You can change the paths on the command line see --help

`python runSeqAnalysis.py --input example_data/2KJV.pdb`

## TODO

1. add msa and nfpar and ffaln output to runSeqAnalysis script

## NEXT UP

1. Write runFFIDP.py
2. Write runConsensus.py
3. Write RSEVAL.py
4. Write master control script (FFIDP.py)
5. Write Docs
6. script to run FF on SGE
7. PDB to Seq script
8. PDB sliding window thing
