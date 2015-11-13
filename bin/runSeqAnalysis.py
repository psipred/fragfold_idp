import sys
import argparse

"""
    1st part of the process. Python script takes a sequence then
    runs PSIPRED and HHBlits and output MSA and SS for Fragfold
"""
parser = argparse.ArgumentParser(description='Runs the preliminary sequence '
                                             'parser for FF-IDP')

parser.add_argument('--input', help="input fasta file")

args = parser.parse_args()

try:
    fasta_file = sys.argv[1]
except:
    parser.print_help()
