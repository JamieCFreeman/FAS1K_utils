#!/usr/bin/env python3


##########################################################################################
# Purpose:
#  In Nexus pipeline, a headerless fasta file is written from the shifted vcf, and
#  then broken into lines of 1000 characters.  For troubleshooting purposes, want
#  a script to confirm correct length of these temporary fasta files.
#
# When run as main, takes an input file path 
# 
##########################################################################################

from validate_fas1k import get_fai_scaff_length
from fas1k_utils import get_chr_string, arm_to_int, read_in

##########################################################################################

def exp_length_arm(arm):
	arm = get_chr_string(file_in)
	exp_len = get_fai_scaff_length(arm_to_int(arm), fai_path)
	return exp_len

def check_fa_length(file_path, exp_len):
	'''
	For a .fasta file get length and match to an expected length
	'''
	# Read in fasta file, if no header all on line1
	seq   = read_in(file_in)
	line1 = seq[0]
	# Is there a header? The fasta files written by the DGN scripts do not 
	#	have a header, but they are standard in fasta files.
	header_present = line1[0] == '>'
	if (header_present):
		line1 = seq[1]
	obs_len = len(line1)
	# Does it match?
	return obs_len == exp_len


##########################################################################################

if __name__ == "__main__":
    import sys
    import os

    file_name = sys.argv[1]
    fai = sys.argv[2] if len(sys.argv) > 2 else "/home/jamie/dmel_ref/DmelRef.fasta.fai"
    
    check_val = 0
    # Check that file is not empty
    if( os.path.getsize(file_name) > 0 ):
        # If passes all three checks, sum=3
        check_val = check_fa_length(file_name, 
#    print(check_val)

    if( check_val != 3 ):
        print(file_name)
