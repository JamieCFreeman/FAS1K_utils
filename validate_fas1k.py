#!/usr/bin/env python3



##########################################################################################

# Usage:
# On the command line, run:
#	python validate_fas1k.py {file_path.fas1k} {optional reference genome fasta index (defaults to my dmel ref if not provided}
# With current defaults returns nothing if all checks pass, returns file name if any checks fail.

# I normally use on whole directory with parallel like so:
# parallel \
# "python validate_fas1k.py {}" \
# :::: <( find /raid10/jamie/EF_genomes/round2/fas1k -iname "*.fas1k" )

##########################################################################################

import os.path
import pandas as pd
from fas1k_utils import *

##########################################################################################

def validate_fas1k(fas1k_file, verbosity=1, ref_fai="", soft_mask=False):
    ''' Validate fas1k file. If a reference genome is provided, compares file length against expected.
    Return depends on verbosity: 
        if 0, returns sum of bool for 3 checks
        if 1, returns a list of booleans for 3 checks:
        if 2, returns descriptive string
    # 1. Check that each line is 1000 characters.
    # 2. Check no illegal characters are present. (By default, not expecting soft masked characters,
            if soft masking present, set soft_mask=True).
    # 3. (Optional) If ref genome provided, check against appropriate scaffold length.
    '''
    ploidy = get_fas1k_ploidy(fas1k_file, soft_mask)
    length = get_fas1k_length(fas1k_file)
    chrom = get_chr_string(fas1k_file)
    out_list = []
    # 1. Check that each line is 1000 characters.
    out_list.append( check_fas1k_wrap(fas1k_file) )
    # 2. Check no illegal characters are present.
    if ( (ploidy == 1) | (ploidy == 2) ):
        out_list.append(True)
    else:
        out_list.append(False)
    # 3. (Optional) If ref genome provided, check against appropriate scaffold length.
    if (len(ref_fai) > 0):
        # Get ref genome scaffold name
        ref_scaff = arm_to_int(chrom)
        # Get scaffold length from the fasta index
        lookup = pd.read_table(ref_fai, header=None)
        ref_length = lookup[lookup[0] == ref_scaff].iloc[:,1].item()
        length_match = ref_length == length
        if (length_match == True):
            out_list.append(True)
        else:
            out_list.append(False)
        #print( "Fas1k length matches ref genome? " + str(length) + "  " + str(length_match))
    if (verbosity == 0):
        return sum(out_list)
    elif (verbosity == 1):
        return out_list
    elif (verbosity == 2):
        s =  "ploidy: " + str(ploidy) + "length: " + str(length) + "correct length for chr " + str(chrom) + str(ref_length)
        return s
    # would be nice to add a more verbose version
    #return "ploidy: " + str(ploidy) + "length: " + str(length) + "chr: " + str(chrom)

def check_fas1k_wrap(fas1k_file):
    ''' Check file is make up of lines of 1000 characters (excluding last). Return True or False. '''
    fas1k = open(fas1k_file,'r')
    lines = fas1k.readlines()
    fas1k.close()
    all_line_len = [len(lines[i].strip()) for i in range(0, len(lines)-1)]
    bool_out =  all([ (all_line_len[i] == 1000) for i in range(0, len(all_line_len)) ])
    return bool_out

def get_fai_scaff_length(ref_scaff, ref_fai):
    lookup = pd.read_table(ref_fai, header=None)
    ref_length = lookup[lookup[0] == ref_scaff].iloc[:,1].item()
    return ref_length

if __name__ == "__main__":
    import sys
    import os

    file_name = sys.argv[1]
    fai = sys.argv[2] if len(sys.argv) > 2 else "/home/jamie/dmel_ref/DmelRef.fasta.fai"

    check_val = 0
    # Check that file is not empty
    if( os.path.getsize(file_name) > 0 ):
        # If passes all three checks, sum=3
        check_val = validate_fas1k(file_name, verbosity=0, ref_fai=fai, soft_mask=False)
#    print(check_val)

    if( check_val != 3 ):
        print(file_name)




