#!/usr/bin/env python3

#
# 2024-04 JCF
#######################################################################
#
# Purpose:
#
# Check .fas1k files for diagnostic SNPs associated with major chr arm
#   inversions and output a table of frequencies for inversion-assoiciated 
#   SNPs
#
# Input:
# 
# "--fas1k" path to fas1k file (script checks fask1k file for all needed arms)
#   by transforming this string- so make sure all arms present in same folder
#
######################################################################

import argparse
import pandas as pd
from fas1k_utils import *
from itertools import compress
import re
from inv_snp_freq import *

####################################################################

parser = argparse.ArgumentParser(description="Calculate inv SNP frequency")  
parser.add_argument('--fas1k', help="Path to a fas1k file (assumes other Chr files are in directory)", 
            required=True, type=str)
parser.add_argument('--out', help="Path to directory where output should be written + outfile name", required=True, type=str)
args = parser.parse_args()

file_now = args.fas1k
out_file_name = args.out

######################################################################

# Check SNPs and write to out file
out = check_all_inv_snps(file_now)
#out_file_name = outdir + "/" + get_name(file_now) + "_inv_SNP_counts.tsv"
#print("Writing to file " + out_file_name)
out.to_csv(out_file_name, sep="\t", index = False) 








