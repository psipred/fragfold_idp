import os
import yaml
import argparse
import sys
import subprocess

"""
    3rd part of the process takes the fragfold ensemble and runs PFClust
    and sliding windwo thing. Outputs FFIDP RMSD profile
"""


def run_exe(args, name):
    """
        Function takes a list of command line args and executes the subprocess.
        Sensibly tries to catch some errors too!
    """
    code = 0
    print("Running "+name)
    try:
        code = subprocess.call(' '.join(args), shell=True)
    except Exception as e:
        print(str(e))
        sys.exit(1)
    if code != 0:
        print(name+" Non Zero Exit status: "+str(code))
        sys.exit(code)


script_path = os.path.dirname(os.path.realpath(__file__))
paths_path = script_path+"/../paths.yml"
if os.path.isfile(paths_path):
    paths_yaml = open(paths_path)
    paths = yaml.load(paths_yaml)
else:
    print("Unable to find paths.yml.n\nShould be in top dir for fragfold_idp")
    exit()

parser = argparse.ArgumentParser(description='Runs PFclust and the sliding '
                                             'window superposition code')
parser.add_argument('--input_name', help="job name for the files generated in "
                                         "in the previous steps, usually "
                                         "an uuid", required=True)
parser.add_argument('--indir',
                    help="Where to find the .ens file. Defaults to the "
                         "runSeqAnalysis.py output directory",
                    default=script_path+"/../output/")
parser.add_argument('--outdir',
                    help="Where to put the output files",
                    default=script_path+"/../output/")

args = parser.parse_args()

rmsdclust_args = [script_path+'/rmsdclust',
                  args.indir+args.input_name+".ens"
                  ]
run_exe(rmsdclust_args, "rmsd clust")


#run rmsdclust.c

#java -jar pfclust indir outdir comma,list,files,
#-cp? and needs some env vars
