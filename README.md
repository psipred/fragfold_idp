# FRAGFOLD_idp

## Introduction

FRAGFOLD-IDP is a complex multi-stage processs with a large number of dependencies. These scripts are designed to help automate both installation and running FF-IDP.

### Installation

We provide an ansible script for installation to any standard linux distro, see below. You must have ansible and python2 available. See https://github.com/ansible/ansible for further details.

###  The basit FF-IDP process

1. Run basic sequences analysis
2. Build an ensemble of FRAGFOLD models
3. Cluster the models and produce a RMSD profile
4. Run DynaMine
5. Build Consensus output
6. Generate evaluation statistics

## How to Install

1. The FFIDP code is implemented in python3. Some of the dependencies required
other runtime environments such as the jvm and python3. You should ensure you have system installs of python2, python3, pip2, pip3, java, C/C++ available.

2. Install python dependencies for your python3 installation or virtualenv

`pip install pyyaml`

`pip install biopython`

`pip install numpy`

`pip install pybrain`

  pip install pybrain for python3 may not work with the structure module
  not correctly installing, instead try these 4 steps:

`wget https://github.com/pybrain/pybrain/archive/master.zip`

`unzip master.zip`

`cd pybrain-master`

`python setup.py install`

3. Install python2 dependencies for your python2 installation ot virtualenv

`pip2 install biopython`

4. Checkout ansible https://github.com/ansible/ansible
5. OPTIONAL: switch to a python2 virtualenv
6. Get a Dynamine API key from http://dynamine.ibsquare.be/download/
7. Edit the paths.yml file to reflect where you want to install the various
   components or where they have already been installed.

8. Run intialisation script for ansible

`source ~/ansible/hacking/env-setup`

9. Change to the FF-idp ansible script directory

`cd ~/FRAGFOLD_idp/ansible`

10. Run FFIDP installation script

`ansible-playbook -i hosts install.yml -f 1`

After this point all the software listed in the paths.yml should be installed
in the listed locations. If this fails or you wish to change anything the
ansible-playbook command can be run repeatedly.

11. If you used a python2 virtualenv to run ansible then switch back to
your python3 evironment once the ansible install has been succesful

## How to run FRAGFOLD-IDP

FRAGFOLD-IDP is composed of several scripts which automate the various steps
required to make the  prediction. There are two ways to run FRAGFOLD-IDP. You
can run each script in an independant step-by-step manner this also gives you
the option of moving some or all of these steps to the cluster/cloud computing
service of your choice. Alternatively we provide a "master" script which will
run all the scripts for you, though we note that as the FRAGFOLD step very
time consuming this master script is not ideal.

Unless otherwise stated we assume python3

### Step-by-step Process

1. Run the Sequence analysis script over each of your fasta files. This will
also produce the FRAGFOLD input files. By default the results go in the
FRAGFOLD_idp/output dir. You can change the paths on the command line see
--help. If running this script on a cluster please ensure --num_threads is
set appropriately

`python runSeqAnalysis.py --input example_data/2KJV.pdb`

2. At this point we have generated the pdb file's MSA and FRAGFOLD input files.
A large ensemble of FRAGFOLD models (>200) needs to be generated. We include
here a python script which will generate a modest set of models (20) and
concatenate them together in to an ensemble of models (.ens file) for the later
step. --input_name should be take from the uuid value generated for the
runSeqAnalysis output files, yours will differ from the example below. If
running this script on a cluster please ensure --num_threads is set
appropriately

`python runFRAGFOLD.py --input_name a15a6b5e-9463-11e6-a62a-989096c13ee6`

This step is non-ideal and you should attempt to generate at least 200 models.
FRAGFOLD is time consuming when generating large numbers of models, we note
that users use a computing cluster or cloud service to run many simultaneous
FRAGFOLD runs. Executing this we must leave as an exercise for the reader.

If you have (sensibly) chosen to run FRAGFOLD on the cluster computing or cloud
platform of your choice then the pdb files will now need to be concatenated
together. A cat command such as the following should suffice:

`cat *.pdb > a15a6b5e-9463-11e6-a62a-989096c13ee6.ens`

3. Step 3 takes the ensemble file which runFRAGFOLD outputted and runs PFclust
and the FRAGFOLD IDP superposition. It will output the FF_IDP RMSD profile

`python runFFIDP.py --input_name a15a6b5e-9463-11e6-a62a-989096c13ee6`

4. Step 1 will have generated a fasta file from your input pdb in the output
directory. You can now use this with the dynamine commandline script to
make a protein dynamics prediction (note where you installed dynamine in the
paths.yml file). Note dynamine is in python2 and requires biopython

`python /opt/dynamine/dynamine.py output/a15a6b5e-9463-11e6-a62a-989096c13ee6.fasta`

5. We would run the consensus predictor over the FF_IDP RMDS profile and the
Dynamine profile. This gives our final prediction values

`python runConsensus.py --input_name a15a6b5e-9463-11e6-a62a-989096c13ee6`

6. If your input PDB file is an NMR file with an ensembl of structures you
can now run the slidingwindow script. SKIP THIS STEP IF YOUR INPUT PDB FILE
IS NOT AN NMR ENSEMBLE

`python SlidingWindow.py --input_name a15a6b5e-9463-11e6-a62a-989096c13ee6 --input_file example_data/2KJV.pdb`

7. Finally we can calculate the correlation between the ensemble available
for benchmarking and diagnositic purposes.

`python RSEVAL.py --input_name a15a6b5e-9463-11e6-a62a-989096c13ee6`

### Single Step Executable

1. The master script has two modes. Report mode and execution mode. In Report
mode the script only outputs the commands needed to run the steps but will
not actually execute each step. You can use this to understand what commands
you need to run should you plan to run each step on a cloud/cluster. In execution
mode the script will run each command.

`python FFIDP.py --input example/2KJV.pdb`
`python FFIDP.py --input 2KJV.pdb --mode execute`

This script has a great number of command line options allowing you to
specifically configure each step.

## TODO

4. Write master control script (FFIDP.py)
    done: added runSeqAnalysis
          added runFRAGFOLD
          added dynamine
    todo: add runFFIDP
          add runConsensus
          add RSEVAL
          add PDBfile slidingwindow superposition

## NEXT UP
5. Write Docs
6. script to run FF on SGE
