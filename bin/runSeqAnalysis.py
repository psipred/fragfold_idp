import sys
import argparse

"""
    1st part of the process. Python script takes a sequence then
    runs PSIPRED and HHBlits and output MSA and SS for Fragfold

    # python runSeqAnalysis.py --input ../example_data/2KJV.seq --blast_dir /cs/research/bioinf/archive0/pfilt_test/ncbi-blast-2.2.31+ --hhsuite_dir /cs/research/bioinf/home1/green/dbuchan/Code/fragfold_idp/opt/hhsuite-2.0.16-linux-x86_64 --psipred_dir /cs/research/bioinf/home1/green/dbuchan/Code/fragfold_idp/opt/psipred --uniref90 /cs/research/bioinf/archive0/pfilt_test/uniref --hhblits_uniref20 /cs/research/bioinf/archive0/ffragfold_idp/hhsuite
"""
parser = argparse.ArgumentParser(description='Runs the preliminary sequence '
                                             'parser for FF-IDP')

parser.add_argument('--input', help="input fasta file")
parser.add_argument('--blast_dir', help="Default location of BLAST+ dir",
                    default="../opt/ncbi-blast-2.2.31+/")
parser.add_argument('--hhsuite_dir', help="Default location of HHSuite",
                    default="../opt/hhsuite-2.0.16-linux-x86_64")
parser.add_argument('--psipred_dir', help="Default location of PSIPRED",
                    default="../opt/psipred")
parser.add_argument('--uniref90', help="Default location of BLAST+ db",
                    default="../opt/hhsuite_db")
parser.add_argument('--hhblits_uniref20', help="Default location hhblits uniref20",
                    default="../opt/uniref")

args = parser.parse_args()

print(args.input)
