import sys
import argparse
import subprocess
import uuid
import yaml
import os

"""
    1st part of the process. Python script takes a sequence then
    runs PSIPRED and HHBlits and output MSA and SS for Fragfold

    # python runSeqAnalysis.py --input ../example_data/2KJV.seq --blast_dir /cs/research/bioinf/archive0/pfilt_test/ncbi-blast-2.2.31+ --hhsuite_dir /cs/research/bioinf/home1/green/dbuchan/Code/fragfold_idp/opt/hhsuite-2.0.16-linux-x86_64 --psipred_dir /cs/research/bioinf/home1/green/dbuchan/Code/fragfold_idp/opt/psipred --uniref90 /cs/research/bioinf/archive0/pfilt_test/uniref/test_db.fasta --hhblits_uniref20 /cs/research/bioinf/archive0/ffragfold_idp/hhsuite/uniprot20_2015_06
    # python runSeqAnalysis.py --input ../example_data/2KJV.seq --uniref90 /cs/research/bioinf/archive0/pfilt_test/uniref/test_db.fasta
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

# Here we grab the paths we set with the ansible install
script_path = os.path.dirname(os.path.realpath(__file__))
paths_path = script_path+"/../paths.yml"
if os.path.isfile(paths_path):
    paths_yaml = open(paths_path)
    paths = yaml.load(paths_yaml)
else:
    print("Unable to find paths.yml.n\nShould be in top dir for fragfold_idp")
    exit()

# Now grab anything from the commandline or set it's defaults
parser = argparse.ArgumentParser(description='Runs the preliminary sequence '
                                             'parser for FF-IDP')

parser.add_argument('--input', help="input pdb file", required=True)
parser.add_argument('--chain',
                    help="The chain in the pdb file you want to analyse",
                    default="A")
parser.add_argument('--blast_dir',
                    help="Default location of BLAST+ dir",
                    default=paths["blast_dir"])
parser.add_argument('--hhsuite_dir',
                    help="Default location of HHSuite",
                    default=paths["hhsuite_dir"])
parser.add_argument('--psipred_dir',
                    help="Default location of PSIPRED",
                    default=paths["psipred_dir"])
parser.add_argument('--uniref90',
                    help="Default location of BLAST+ db",
                    default=paths["uniref90_dir"]+"/uniref90.fasta")
parser.add_argument('--hhblits_uniref20',
                    help="Default location hhblits uniref20",
                    default=paths["hhdb_dir"]+"/"+paths["hhdbversion"])
parser.add_argument('--uuid',
                    help="UUID from previous step or controlling script",
                    default=str(uuid.uuid1()))
parser.add_argument('--threads',
                    help="Number of threads to use for BLAST and HHBlits",
                    default="10")
parser.add_argument('--outdir',
                    help="Where to put the output files",
                    default=script_path+"/../output/")

args = parser.parse_args()

fastapdbv3_args = [script_path+'/fasta_pdbv3.pl',
                   args.input,
                   args.outdir+args.uuid+".fasta",
                   args.chain,
                   ]
run_exe(fastapdbv3_args, "Converting PDB file")

# Going to run BLAST+
blast_args = [args.blast_dir+"/bin/psiblast",
              "-query", args.outdir+args.uuid+".fasta",
              "-db", args.uniref90,
              "-inclusion_ethresh", "0.001",
              "-out_pssm", args.outdir+args.uuid+".pssm",
              "-out", args.outdir+args.uuid+".bls",
              "-num_iterations", "3",
              "-num_alignments", "0",
              "-num_descriptions", "5000",
              "-num_threads", args.threads,
              ]
run_exe(blast_args, "BLAST+")

# Going to run PSIPRED components
chkparse_args = [args.psipred_dir+"/bin/chkparse",
                 args.outdir+args.uuid+".pssm",
                 ">",
                 args.outdir+args.uuid+".mtx",
                 ]
run_exe(chkparse_args, "chkparse")

psipred_args = [args.psipred_dir+"/bin/psipred",
                args.outdir+args.uuid+".mtx",
                args.psipred_dir+"/data/weights.dat",
                args.psipred_dir+"/data/weights.dat2",
                args.psipred_dir+"/data/weights.dat3",
                ">",
                args.outdir+args.uuid+".ss",
                ]
run_exe(psipred_args, "Psipred")

psipass_args = [args.psipred_dir+"/bin/psipass2",
                args.psipred_dir+"/data/weights_p2.dat",
                "1",
                "1.0",
                "1.0",
                args.outdir+args.uuid+".ss2",
                args.outdir+args.uuid+".ss",
                ">",
                args.outdir+args.uuid+".horiz",
                ]
run_exe(psipass_args, "Psipass2")

# Going to run HHBlits
hhblits_args = [args.hhsuite_dir+"/bin/hhblits",
                "-i", args.outdir+args.uuid+".fasta",
                "-o", args.outdir+args.uuid+".hh",
                "-oa3m", args.outdir+args.uuid+".a3m",
                "-n", "3",
                "-d", args.hhblits_uniref20+"/"+paths["hhdbversion"],
                "-cpu", args.threads,
                ]
run_exe(hhblits_args, "hhblits")

# Going to proces msa
msa_args = ['egrep', '-v', '"^>"', args.outdir+args.uuid+".a3m",
            '|',
            'sed', "'s/[a-z]//g'",
            '>',
            args.outdir+args.uuid+".msa",
            ]
# print(' '.join(msa_args))
run_exe(msa_args, "MSA_GEN")

# Here we basically replicate the old do_ffaln.sh script
createffaln_args = ['echo', '-e',
                    args.uuid,
                    '>>',
                    args.outdir+args.uuid+".ffaln",
                    ]
run_exe(createffaln_args, "Creating ffaln")

addlines_args = ['wc', '-l', args.outdir+args.uuid+".msa",
                 '|',
                 'awk', "'{print $1}'",
                 '>>',
                 args.outdir+args.uuid+".ffaln",
                 ]
run_exe(addlines_args, "Add msa count")

# TODO: this should be replaced with psfilt!!!!
processSS2_args = [script_path+'/psifilt',
                   '<',
                   args.outdir+args.uuid+".ss2",
                   ">>",
                   args.outdir+args.uuid+".ffaln",
                   ]
run_exe(processSS2_args, "Add SS")

seq_args = ['grep', '-v',
            '"^>"',
            args.outdir+args.uuid+".fasta",
            '|',
            'tr', '-d', "'\\n'",
            '|',
            'sed', "'s/$/\\n/'",
            '>>',
            args.outdir+args.uuid+".ffaln",
            ]
# print(' '.join(seq_args))
run_exe(seq_args, "Add Sequence")

addmsa_args = ['cat', args.outdir+args.uuid+".msa",
               '>>',
               args.outdir+args.uuid+".ffaln",
               ]
run_exe(addmsa_args, "Add MSA")

#Output the fragfold file
print("Printing nfpar file")
nfpar_string = "# FRAGFOLD Parameter File\n\n# Alignment file\n"
nfpar_string += "ALNFILE "+args.outdir+args.uuid+".ffaln\n\n# contacts file\n\n"
nfpar_string += "# Weighting mode\nWTMODE STDEV\n"
nfpar_string += "# Steric mode\nSTERMODE ALLATOM\n"
nfpar_string += "# Short range weighting\nSRWT 1.0\n"
nfpar_string += "# Long range weighting\nLRWT 1.0\n"
nfpar_string += "# Solvation weighting\nSOLVWT 1.0\n"
nfpar_string += "# Hydrogen bond weighting\nHBWT 1.0\n"
nfpar_string += "# Compactness weighting\nCOMPACTWT 1.0\n"
nfpar_string += "# Steric weighting\nSTERICWT 3.0\n"
nfpar_string += "# Disulphide weighting\nDSWT 0.0\n"
nfpar_string += "# RMSD-to-target weighting\nTARGWT 0.0\n"
nfpar_string += "# Residue-residue contact weighting\nRRWT 0.0\n"
nfpar_string += "# Max. annealing steps\nMAXSTEPS 10000000\n"
nfpar_string += "# Initial temperature\nINITEMP 0.6\n"
nfpar_string += "# MAXFRAGS\nMAXFRAGS 5\n"
nfpar_string += "# MAXFRAGS2\nMAXFRAGS2 25\n"
nfpar_string += "MODE FOLD\n"
nfpar_string += "POOLSIZE 10\n"
nfpar_string += "TRATIO 0.6\n"
nfpar_string += "OPTMETHOD REPMC\n"
nfpar = open(args.outdir+args.uuid+".nfpar", 'w')
nfpar.write(nfpar_string)
nfpar.close()
