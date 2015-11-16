import sys
import argparse
import subprocess
import uuid

"""
    1st part of the process. Python script takes a sequence then
    runs PSIPRED and HHBlits and output MSA and SS for Fragfold

    # python runSeqAnalysis.py --input ../example_data/2KJV.seq --blast_dir /cs/research/bioinf/archive0/pfilt_test/ncbi-blast-2.2.31+ --hhsuite_dir /cs/research/bioinf/home1/green/dbuchan/Code/fragfold_idp/opt/hhsuite-2.0.16-linux-x86_64 --psipred_dir /cs/research/bioinf/home1/green/dbuchan/Code/fragfold_idp/opt/psipred --uniref90 /cs/research/bioinf/archive0/pfilt_test/uniref/test_db.fasta --hhblits_uniref20 /cs/research/bioinf/archive0/ffragfold_idp/hhsuite/uniprot20_2015_06
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


parser = argparse.ArgumentParser(description='Runs the preliminary sequence '
                                             'parser for FF-IDP')

parser.add_argument('--input', help="input fasta file")
parser.add_argument('--blast_dir',
                    help="Default location of BLAST+ dir",
                    default="../opt/ncbi-blast-2.2.31+/")
parser.add_argument('--hhsuite_dir',
                    help="Default location of HHSuite",
                    default="../opt/hhsuite-2.0.16-linux-x86_64")
parser.add_argument('--psipred_dir',
                    help="Default location of PSIPRED",
                    default="../opt/psipred")
parser.add_argument('--uniref90',
                    help="Default location of BLAST+ db",
                    default="../opt/hhsuite_db/uniref90.fasta")
parser.add_argument('--hhblits_uniref20',
                    help="Default location hhblits uniref20",
                    default="../opt/uniprot20_2015_06/uniprot20_2015_06")
parser.add_argument('--uuid',
                    help="UUID from previous step or controlling script",
                    default=str(uuid.uuid1()))
args = parser.parse_args()
# Going to run BLAST+
blast_args = [args.blast_dir+"/bin/psiblast",
              "-query", args.input,
              "-db", args.uniref90,
              "-inclusion_ethresh", "0.001",
              "-out_pssm", "../output/"+args.uuid+".pssm",
              "-out", "../output/"+args.uuid+".bls",
              "-num_iterations", "3",
              "-num_alignments", "0",
              ]
run_exe(blast_args, "BLAST+")

# Going to run PSIPRED components
chkparse_args = [args.psipred_dir+"/bin/chkparse",
                 "../output/"+args.uuid+".pssm",
                 ">",
                 "../output/"+args.uuid+".mtx",
                 ]
run_exe(chkparse_args, "chkparse")

psipred_args = [args.psipred_dir+"/bin/psipred",
                "../output/"+args.uuid+".mtx",
                args.psipred_dir+"/data/weights.dat",
                args.psipred_dir+"/data/weights.dat2",
                args.psipred_dir+"/data/weights.dat3",
                ">",
                "../output/"+args.uuid+".ss",
                ]
run_exe(psipred_args, "Psipred")

psipass_args = [args.psipred_dir+"/bin/psipass2",
                args.psipred_dir+"/data/weights_p2.dat",
                "1",
                "1.0",
                "1.0",
                "../output/"+args.uuid+".ss2",
                "../output/"+args.uuid+".ss",
                ">",
                "../output/"+args.uuid+".horiz",
                ]
run_exe(psipass_args, "Psipass2")

# Going to run HHBlits
hhblits_args = [args.hhsuite_dir+"/bin/hhblits",
                "-i", args.input,
                "-o", args.uuid+".hh",
                "-n", "4",
                "-d", args.hhblits_uniref20,
                ]
run_exe(hhblits_args, "hhblits")
