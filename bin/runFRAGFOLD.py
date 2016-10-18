from __future__ import print_function
import os
import yaml
import argparse
import sys
import subprocess
from multiprocessing import Pool
from functools import partial
import glob

"""
    2nd part of the process. Python script takes the job md5 generated in the
    previous step and runs several fragfold runs and then cats together
    the pdb files

    # python runFRAGFOLD.py --input_name a15a6b5e-9463-11e6-a62a-989096c13ee6
"""


def eprint(*args, **kwargs):
    print(*args, file=sys.stderr, **kwargs)


def run_fragfold(args, index):
    exe = args.fragfold_dir+"/bin/fragfold"
    fdata_dir = args.fragfold_dir+"/data/"
    ftdb_dir = args.fragfold_dir+"/data/"
    print(index)

    ff_args = ['FDATA_DIR='+fdata_dir,
               'FTDB_DIR='+ftdb_dir,
               args.fragfold_dir+"/bin/fragfold",
               args.indir+args.input_name+'.nfpar',
               args.outdir+args.input_name+"_ffmodels/"+str(index)+".pdb"]
    print(" ".join(ff_args))
    try:
        code = subprocess.call(' '.join(ff_args), shell=True)
    except Exception as e:
        print(str(e))
        sys.exit(1)
    if code != 0:
        print(name+" Non Zero Exit status: "+str(code))
        sys.exit(code)


def build_ensemble(args):
    pdb_files = sorted(glob.glob(args.outdir+args.input_name+"_ffmodels/*.pdb"))
    ens = open(args.outdir+args.input_name+".ens", 'w')

    for i in range(0, len(pdb_files)):
        ens.write("MODEL "+str(i)+"\n")
        with open(pdb_files[i], "r") as inputfile:
            contents = inputfile.read()
        ens.write(contents)
    ens.close()


script_path = os.path.dirname(os.path.realpath(__file__))
paths_path = script_path+"/../paths.yml"
if os.path.isfile(paths_path):
    paths_yaml = open(paths_path)
    paths = yaml.load(paths_yaml)
else:
    print("Unable to find paths.yml.n\nShould be in top dir for fragfold_idp")
    exit()

parser = argparse.ArgumentParser(description='Runs FRAGFOLD multiple times  '
                                             'and builds the pdb model '
                                             'ensemble')
parser.add_argument('--input_name', help="job name for the files generated in "
                                         "in the seq analysis step, usually "
                                         "an md5", required=True)
parser.add_argument('--indir',
                    help="Where to find the .nfpar file. Defaults to the "
                         "runSeqAnalysis.py output directory",
                    default=script_path+"/../output/")
parser.add_argument('--outdir',
                    help="Where to put the output files",
                    default=script_path+"/../output/")
parser.add_argument('--num_threads',
                    help="Number of concurrent fragfold processes to run",
                    default="10")
parser.add_argument('--num_models',
                    help="The number of fragfold runs to complete",
                    default="20")
parser.add_argument('--large_override',
                    help="Override 20 model limit",
                    default="FALSE")
parser.add_argument('--fragfold_dir',
                    help="installation location for FRAGFOLD",
                    default=paths["fragfold_dir"])

args = parser.parse_args()

eprint("WARNING: This script will only generate 20 fragfold models. An "
       "accurate fragfold_idp run requires that you generate >200 models. We "
       "advise users use an available cloud or cluster service to achieve this "
       "if you wish to use this script to generate more than 20 models "
       "you will need to set --large_override and --num_models\n")

if int(args.num_models) > 20 and 'TRUE' not in args.large_override:
    eprint("You have set --num_models greater than 20. These fragfold runs "
           "will take a very long time. If you are sure you wish to run this "
           "number of models set --large_override to TRUE")
    exit()

if not os.path.isdir(args.outdir+args.input_name+"_ffmodels"):
    os.makedirs(args.outdir+args.input_name+"_ffmodels")

seqs = {}
p = Pool(int(args.num_threads))
args_set = [args]
indices = range(0, int(args.num_models))
runff_partial = partial(run_fragfold, args)
p.map(runff_partial, indices)
p.close()
p.join()

print("Building model ensemble")
build_ensemble(args)
