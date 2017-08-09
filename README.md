# FRAGFOLD_idp

## Citation
This software and data related to the following publication:  
T. Kosciolek, D.W.A. Buchan & D.T. Jones (2017) "Predictions of Backbone Dynamics in Intrinsically Disordered Proteins Using De Novo Fragment-Based Protein Structure Predictions." Scientific Reports 7, 6999 (https://www.nature.com/articles/s41598-017-07156-1)

## Introduction

FRAGFOLD-IDP (FFIDP) calculate protein chain backbone dynamics, from which you can infer which residues in a protein may be intrinsically disordered. FFIDP is a complex multi-stage processs with a large number of dependencies. These scripts are designed to help automate both installation and running FFIDP.

###  The basic FFIDP process

1. Run basic sequences analysis
2. Build an ensemble of FRAGFOLD models
3. Cluster the models and produce a RMSD profile
4. Run DynaMine
5. Build Consensus output
6. Generate evaluation statistics

## Installation

We provide an ansible script for installation to any standard linux distro, see below.
You must have ansible and python2 available. See https://github.com/ansible/ansible for further details.

### Install Steps

1. The FFIDP code is implemented in python3. Some of the dependencies require other runtime environments such as the jvm and python2. You should ensure you have system installs of tar, python2, python3, pip2, pip3, java and C/C++ available.

2. Enter the src/ directory and run the following

    `make clean; make; make install`

3. Install python3 dependencies for your python3 installation or virtualenv

    `pip install -r requirements.txt`

    Note: pip install pybrain for python3 may not work with the structure module not correctly installing or runConsensus.py (below) may not run, instead try these 4 steps:

    `wget https://github.com/pybrain/pybrain/archive/master.zip`

    `unzip master.zip`

    `cd pybrain-master`

    `python setup.py install`

4. Checkout ansible

    `git clone https://github.com/ansible/ansible`

5. OPTIONAL: You may wish to use a python2 virtualenv for steps 6 to 10.  

6. Install python2 dependencies

    `pip2 install biopython`
    `pip2 install jinja2`

7. Get a DynaMine API key from http://dynamine.ibsquare.be/download/

8. Edit the paths.yml file to reflect where you want to install the various
   components or where they have already been installed. Ensure you have the
   version for blast and hhdb for your operating system. Note that you are
   free to update to more recent versions HHBlits or Blast+

9. Run the ansible intialisation script provided by the ansible package

    `source [path_to_ansible]/ansible/hacking/env-setup`

10. Change to the FFIDP ansible script directory

    `cd [ff_idp_root]/fragfold_idp/ansible`

11. Run FFIDP installation script

    `ansible-playbook -i hosts install.yml -f 1`

    After this point all the software listed in the paths.yml should be installed
    in the listed locations. If this fails or you wish to change anything the
    `ansible-playbook` command can be run repeatedly.

12. If you used a python2 virtualenv to run ansible then switch back to your python3 evironment once the ansible install has been successful and return to the root fragfold-idp directory.

## How to run FRAGFOLD-IDP

FRAGFOLD-IDP is composed of several scripts which automate the various steps
required to make the  prediction. There are two ways to run FRAGFOLD-IDP. You
can run each script in an independant step-by-step manner this also gives you
the option of moving some or all of these steps to the cluster/cloud computing
service of your choice. Alternatively we provide a "master" script which will
run all the scripts for you, though we note that as the FRAGFOLD step very
time consuming this master script is not ideal.

**Unless otherwise stated we assume Python3**

If you are only interested in making an FFIDP protein backbone dynamics prediction then
you can stop at the end of Step 3. If you wish to make a full FFIDP-dynamine consensus prediction you must continue until at least the end of step 6.

### Step-by-step Process

1. Run the Sequence analysis script over each of your fasta files. This will
also produce the FRAGFOLD input files. By default the results go in the
FRAGFOLD_idp/output dir. You can change the paths on the command line see
--help. If running this script on a cluster please ensure --num_threads is
set appropriately

    `python bin/runSeqAnalysis.py --input example_data/2KJV.pdb`

2. At this point we have generated the pdb file's MSA and FRAGFOLD input files.
A large ensemble of FRAGFOLD models (>200) needs to be generated. We include
here a python script which will generate a modest set of models (20) and
concatenate them together in to an ensemble of models (.ens file) for the later
step. --input_name should be take from the uuid value generated for the
runSeqAnalysis output files, yours will differ from the example below. If
running this script on a cluster please ensure --num_threads is set
appropriately

    `python bin/runFRAGFOLD.py --input_name a15a6b5e-9463-11e6-a62a-989096c13ee6`

  This step is non-ideal and you should attempt to generate at least 200models.
  FRAGFOLD is time consuming when generating large numbers of models, we note
  that users should use a computing cluster or cloud service to run many simultaneous FRAGFOLD runs. Arranging this we must leave as an exercise for
  the user.

  Fragfold requires as input the .nfpar and .ffaln file. In paths.ymlensure you set ffaln_dir a the directory on your cluster that thefragfold ptocess can see.  This is where you need to move the .ffaln and.nfpar files. Setting this
  path correctly ensures the runSeqAnalaysis will write correct .nfparfile.  

  If you have (sensibly) chosen to run FRAGFOLD on the cluster computing orcloud
  platform of your choice then the pdb files will now need to beconcatenated
  together. A cat command such as the following should suffice:

    `cat *.pdb > a15a6b5e-9463-11e6-a62a-989096c13ee6.ens`

  The .ens file should then be placed in the output directory you have been using
  (default is output/)

3. Step 3 takes the ensemble file which runFRAGFOLD outputted and runs PFclust
and the FRAGFOLD IDP superposition. It will output the FF_IDP RMSD profile

    `python bin/runFFIDP.py --input_name a15a6b5e-9463-11e6-a62a-989096c13ee6`

4. Step 1 will have generated a fasta file from your input pdb in the output
  directory. You can now use this with the DynaMine commandline script to
  make a protein dynamics prediction (note where you installed DynaMine in the
  paths.yml file). Note, *DynaMine is in python2 and requires biopython*
  If you are using a python2 virtualenv switch to it for this step

    `python2 /opt/dynamine/dynamine.py output/a15a6b5e-9463-11e6-a62a-989096c13ee6.fasta`

  NOTE: It is critical that dynamine completes, do not continue to steps 5-7
  unless dynamine has produced it's output

5. We would run the consensus predictor over the FF_IDP RMDS profile and the
Dynamine profile. This gives our final prediction values

    `python bin/runConsensus.py --input_name a15a6b5e-9463-11e6-a62a-989096c13ee6`

6. If your input PDB file is an NMR file with an ensemble of structures you
can now run the slidingwindow script. SKIP THIS STEP IF YOUR INPUT PDB FILE
IS NOT AN NMR ENSEMBLE

    `python bin/SlidingWindow.py --input_name a15a6b5e-9463-11e6-a62a-989096c13ee6 --input_file example_data/2KJV.pdb`

7. Finally we can calculate the correlation between the ensemble available
for benchmarking and diagnostic purposes. Assuming the outputs are in the
results directory the previous scripts have all been using you provide
the names of the output files. If you want to move these to another
directory then provide the --results_dir flag

    `python bin/RSEVAL.py -i a15a6b5e-9463-11e6-a62a-989096c13ee6.pdb_ens -j a15a6b5e-9463-11e6-a62a-989096c13ee6.consensus`

### Single Step Executable

1. The master script has two modes. Report mode and execution mode. In Report
mode the script only outputs the commands needed to run the steps but will
not actually execute each step. You can use this to understand what commands
you need to run should you plan to run each step on a cloud/cluster. In execution
mode the script will run each command.

    `python bin/FFIDP.py --input example/2KJV.pdb`
    `python bin/FFIDP.py --input 2KJV.pdb --mode execute`

    *This script has a great number of command line options allowing you to
    specifically configure each step.*

## Outputs

FFIDP produces a large number of intermediary files while executing each step.
The important files are the outputs of runFFIDP (`.ffidp`) and runConsensus files (`.consensus`). They represent predicted per-residue RMSD (Å) values. To compare FRAGFOLD-IDP and DynaMine results using e.g. `RSEVAL.py`, we recommend computing `1-S^2`. Sample FRAGFOLD-IDP output is in the form:
```
A	0.581
R	1.620
Y	1.722
E	1.535
```
where the first column indicates amino acid, as inferred from the input sequence, and the second column is the predicted per-residue RMSD. 

## NEXT UP
1. Write Docs
