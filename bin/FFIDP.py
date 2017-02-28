import os
import yaml
import argparse
import sys
import uuid
import glob
import itertools
from subprocess import call

"""
    Master scripts
    runs the following in order

    1. runSeqAnalysis.py *
    2. Fragfold
    3. runFFIDP.py
    4. DynaMine
    5. runConsensus.py
    6. RSEVAL.py
"""

script_path = os.path.dirname(os.path.realpath(__file__))
paths_path = script_path+"/../paths.yml"
if os.path.isfile(paths_path):
    paths_yaml = open(paths_path)
    paths = yaml.load(paths_yaml)
else:
    print("Unable to find paths.yml.n\nShould be in top dir for fragfold_idp")
    exit()

parser = argparse.ArgumentParser(description='Runs FRAGFOLD-IDP')
parser.add_argument('--input', help="PDB file for input", required=True)
parser.add_argument('--mode', help="Control whether the commands are printed"
                                   "or executed",
                                   default="print")

# Seq Analysis arguments
parser.add_argument('--chain',
                    help="The chain in the pdb file you want to analyse",
                    default="A")
parser.add_argument('--name', help="Optional job name, "
                                   "otherwise a UUID will be generated and"
                                   " used",
                    default=str(uuid.uuid1()))

# SeqAnalysis arguments
parser.add_argument('--blast_dir',
                    help="Default location of BLAST+ dir for runSeqAnalysis",
                    default=paths["blast_dir"])
parser.add_argument('--hhsuite_dir',
                    help="Default location of HHSuite for runSeqAnalysis",
                    default=paths["hhsuite_dir"])
parser.add_argument('--psipred_dir',
                    help="Default location of PSIPRED for runSeqAnalysis",
                    default=paths["psipred_dir"])
parser.add_argument('--uniref90',
                    help="Default location of BLAST+ db for runSeqAnalysis",
                    default=paths["uniref90_dir"]+"/uniref90.fasta")
parser.add_argument('--hhblits_uniref20',
                    help="Default location hhblits uniref20 for runSeqAnalysis",
                    default=paths["hhdb_dir"]+"/"+paths["hhdbversion"])
parser.add_argument('--uuid',
                    help="UUID from previous step or controlling script",
                    default=str(uuid.uuid1()))
parser.add_argument('--threads',
                    help="Number of threads to use for BLAST and HHBlits and "
                         "Fragfold",
                    default="10")
parser.add_argument('--outdir',
                    help="Where to put the output files",
                    default=script_path+"/../output/")

# FF arguments
parser.add_argument('--indir',
                    help="Where to find the .nfpar file. Defaults to the "
                         "--outdir directory, used by runFRAFOLD",
                    default=script_path+"/../output/")
# parser.add_argument('--outdir',
#                     help="Where to put the output files",
#                     default=script_path+"/../output/")
# parser.add_argument('--num_threads',
#                     help="Number of concurrent fragfold processes to run",
#                     default="10")
parser.add_argument('--num_models',
                    help="The number of fragfold runs to complete, used by "
                         "runFRAGFOLD",
                    default="20")
parser.add_argument('--large_override',
                    help="Override 20 model limit, used by runFRAGFOLD",
                    default="TRUE")
parser.add_argument('--fragfold_dir',
                    help="installation location for FRAGFOLD, used by "
                         "runFRAGFOLD",
                    default=paths["fragfold_dir"])

##Dynamine arguments
parser.add_argument('--dynamine_dir',
                    help="installation location for dynamine, used by ",
                    default=paths["dynamine_dir"])

##FFIDP arguments

##runConsensus arguments
parser.add_argument('--ffidp_profile_path',
                    help="The path for the FFIDP profile output for "
                         "runConsensus",
                    default=script_path+"/../output/")
parser.add_argument('--dynamine_profile_path',
                    help="The path for the dynamine profile output for "
                          "runConsensus",
                    default=script_path+"/../output/")
parser.add_argument('--psipred_result_path',
                    help="The path for the psipred ss output for runConsensus",
                    default=script_path+"/../output/")
parser.add_argument('--alignment_path',
                    help="The path for the MSA for runConsensus",
                    default=script_path+"/../output/")


##SlidingWindow arguments
parser.add_argument('--window_size',
                    help="size of the sliding window for SlidingWindow.py",
                    default=10)

args = parser.parse_args()

# Construct the runSeqAnalysis.py command 1st step
seq_args = ["python", script_path+"/runSeqAnalysis.py",
            '--input', args.input,
            '--uuid', args.name,
            ]

if "--blast_dir" in sys.argv:
    seq_args += ['--blast_dir', args.blast_dir]
if "--hhsuite_dir" in sys.argv:
    seq_args += ['--hhsuite_dir', args.hhsuite_dir]
if "--psipred_dir" in sys.argv:
    seq_args += ['--psipred_dir', args.psipred_dir]
if "--uniref90" in sys.argv:
    seq_args += ['--uniref90', args.uniref90]
if "--hhblits_uniref20" in sys.argv:
    seq_args += ['--hhblits_uniref20', args.hhblits_uniref20]
if "--threads" in sys.argv:
    seq_args += ['--threads', args.threads]
if "--outdir" in sys.argv:
    seq_args += ['--outdir', args.outdir]

print("#### runSeqAnalysis Command ####")
print(" ".join(seq_args)+"\n")

# Construct the runFRAGFOLD.py command 2nd step
ff_args = ["python", script_path+"/runFRAGFOLD.py",
           "--input_name", args.name,
           "--large_override", args.large_override,
           ]

if "--outdir" in sys.argv:
    ff_args += ['--outdir', args.outdir]
if "--indir" in sys.argv:
    ff_args += ['--indir', args.indir]
if "--threads" in sys.argv:
    ff_args += ['--num_threads', args.threads]
if "--num_models" in sys.argv:
    ff_args += ['--num_models', args.num_models]
if "--fragfold_dir" in sys.argv:
    ff_args += ['--fragfold_dir', args.fragfold_dir]

print("#### runFRAGFOLD Command ####")
print(" ".join(ff_args)+"\n")

# Construct the runFFIDP command 3rd step
ffidp_args = ["python", script_path+"/runFFIDP.py",
              "--input_name", args.name,
              ]
if "--outdir" in sys.argv:
    ffidp_args += ['--outdir', args.outdir]
if "--indir" in sys.argv:
    ffidp_args += ['--indir', args.indir]

print("#### runFFIDP Command ####")
print(" ".join(ffidp_args)+"\n")

# Construct the dynamine command 4th step
dynamine_args = ["python2", args.dynamine_dir+"/dynamine.py",
                 args.outdir+args.name+".fasta",
                 ]

print("#### Dynamine Command ####")
print(" ".join(dynamine_args)+"\n")

# Construct the runConsensus 5th step
consensus_args = ["python", script_path+"/runConsensus.py",
                  "--input_name", args.name,
                  ]
if "--outdir" in sys.argv:
    consensus_args += ['--outdir', args.outdir]
if "--indir" in sys.argv:
    consensus_args += ['--indir', args.indir]
if "--ffidp_profile_path" in sys.argv:
    consensus_args += ['--ffidp_path', args.ffidp_profile_path]
if "--dynamine_profile_path" in sys.argv:
    consensus_args += ['--dynamine_path', args.dynamine_profile_path]
if "--psipred_result_path" in sys.argv:
    consensus_args += ['--psipred_path', args.psipred_result_path]
if "--alignment_path" in sys.argv:
    consensus_args += ['--alignment_path', args.alignment_path]

print("#### runConsensus Command ####")
print(" ".join(consensus_args)+"\n")

# Construct the run sliding window command
sliding_args = ["python", script_path+"/SlidingWindow.py",
                "--input_name", args.name,
                "--input_file", args.input]
if "--outdir" in sys.argv:
    sliding_args += ['--outdir', args.outdir]
if "--window_size" in sys.argv:
    sliding_args += ['--window_size', args.window_size]

print("#### SlidingWindow Command ####")
print(" ".join(sliding_args)+"\n")


print("#### RSEVAL Commands ####")
profiles = [args.outdir+"/"+args.name+".ffidp",  # FFIDP profile
            args.outdir+"/"+args.name+".consensus",  # Consensus profile
            args.outdir+"/"+args.name+".pdb_ens",  # PDB profile
            ]
if glob.glob(args.outdir+"/"+args.name+"/Dynamine_b_*"):
    dynaResults = glob.glob(args.outdir+"/"+args.name+"/Dynamine_b_*")[0]
    dynamineOutput = ''
    for dynaFile in glob.glob(dynaResults+"/*"):
        if ".pred" in dynaFile:
            dynamineOutput = dynaFile
    profiles.append(dynamineOutput)

for pair in itertools.combinations(profiles, 2):
    if pair[0] == pair[1]:
        continue
    rseval_args = ['python', script_path+"/RSEVAL.py",
                   '-i', pair[0],
                   '-j', pair[1], ]
    print(" ".join(rseval_args)+"\n")

# Construct RSEVAL command
if args.mode == 'execute':
    call(seq_args)
    call(ff_args)
    call(ffidp_args)
    call(dynamine_args)
    call(consensus_args)
    call(sliding_args)
    print("RSEVAL calls are not executed. These should be run manually if "
          "needed")
