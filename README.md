# fragfold_idp

## Introduction

FRAGFOLD-IDP is a complex multi stages processs with a large number of dependencies. These scripts are designed to help automate installation and running FF-IDP.

### Installation

We provide an ansible script for installation to any standard linux distro, se below. You must have ansible and python2 available. See https://github.com/ansible/ansible for further details.

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

  ```pip install pyyaml```

4. Edit the paths.yml to reflect where you want to install the various
   components or where they have already been installed.
5. Run setup script for ansible

    ```source ~/ansible/hacking/env-setup```

6. Change to the ansible script directory

    ```cd ~/fragfold_idp/ansible```

7. Run FFIDP installation script

    ```ansible-playbook -i hosts install.yml -f 1```

## How to run

## TODO

1. Keep writing runFFIDP.py


## NEXT UP

2. Write runFFIDP.py
3. Write runConsensus.py
4. Write RSEVAL.py
5. Write master control script (FFIDP.py)
7. Write Docs
8. script to run FF on SGE
9. PDB to Seq script
10. PDB sliding window thing
