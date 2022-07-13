#!/usr/bin/env python3

import os.path
import pandas as pd

def subseq_from_lines(start, end, line_len, lines):
    ''' 
    Helper function for extract_strain_subseq, you probably don't want to call this.
    # Functions stolen from Chris to get nt of interest from fas1k files.
    # Specifically from "/raid10/chris/amplicon_design/primer_gen/inv_primer_gen.py"
    '''
    # Retrieve the sequence
    curr_line = start // line_len
    start_index = start %  line_len
    end_line = end // line_len
    end_index = end % line_len
    seq = ''
    if curr_line == end_line:
        seq += lines[curr_line].strip()[start_index:end_index]
    else:
        seq += lines[curr_line].strip()[start_index:]
        curr_line += 1
        while curr_line < end_line:
            seq += lines[curr_line].strip()
            curr_line += 1
            seq += lines[curr_line].strip()[:end_index]
    return(seq)

def extract_fas1k_subseq(start, end, fas1k):
    '''
    Extract subsequence from fas1k file. Coordinates in Python numbering space, e.g. from 0:len-1
    # Functions stolen from Chris to get nt of interest from fas1k files.
    # Specifically from "/raid10/chris/amplicon_design/primer_gen/inv_primer_gen.py"
    '''
    # Prep relevant values
    strain = open(strain_file,'r')
    lines = strain.readlines()
    strain.close()
    line_len = len(lines[0].strip())
    seq = subseq_from_lines(start,end,line_len,lines)
    return(seq)

def get_fas1k_length(fas1k_file):
    ''' Returns nt length of fas1k file..'''
    fas1k = open(fas1k_file,'r')
    lines = fas1k.readlines()
    fas1k.close()
    line_len = len(lines[0].strip())
    n_lines = len(lines) - 1
    last_line = len(lines[len(lines)-1].strip())
    return (n_lines*line_len) + last_line

def get_fas1k_ploidy(fas1k_file, soft_mask=False):
    ''' Compare all characters present in file to those expected based on ploidy, return ploidy as 1 or 2.
    If you know soft masking has been used, set soft_mask=True. '''
    # When using existing fas1k files, make sure you know what filtering has been performed
    # List legal characters for haploid and diploid
    diploid_set = set(('A', 'T', 'C', 'G', 'N', 'Y', 'R', 'W', 'S', 'K', 'M'))
    haploid_set = set(('A', 'T', 'C', 'G', 'N'))
    # Get set of all chars in file (excluding /n)
    all_char_set = get_fas1k_char(fas1k_file)
    # Compare set of characters from file to haploid, diploid, and return result
    if (all_char_set.union(haploid_set) == haploid_set):
        return 1
    elif (all_char_set.union(diploid_set) == diploid_set):
        return 2
    # Only look for soft-masked nts when argument soft_mask is set to true
    if (soft_mask == True):
        softmask_diploid_set = set(('a', 't', 'c', 'g', 'n', 'y', 'r', 'w', 's', 'k', 'm'))
        softmask_haploid_set = set(('a', 't', 'c', 'g', 'n'))
        if (all_char_set.union(haploid_set).union(softmask_haploid_set) == haploid_set.union(softmask_haploid_set)):
            return 1
        elif (all_char_set.union(diploid_set).union(softmask_diploid_set) == diploid_set.union(softmask_diploid_set)):
            return 2
    # If doesn't match either, return unexpected characters
    # (uncovered edge case- fas1k that doesn't have any occurance of 1 nt- 
    #   would only expect this to occur if you were trying to use it on a small
    #   slice of a fas1k file.)
    else:
        return "Unexpected characters encountered: " + str(all_char_set.difference(diploid_set))

def get_fas1k_char(fas1k_file):
    ''' Returns a python set of all characters present in file. '''
    fas1k = open(fas1k_file,'r')
    lines = fas1k.readlines()
    fas1k.close()
    all_char_set = set(lines[0].strip())
    for i in range(len(lines)-1):
        all_char_set.update(set(lines[i].strip()))
    return all_char_set

def validate_fas1k(fas1k_file, verbosity=1, ref_fai="", soft_mask=False):
    ''' Validate fas1k file. If a reference genome is provided, compares file length against expected. '''
    ploidy = get_fas1k_ploidy(fas1k_file, soft_mask)
    length = get_fas1k_length(fas1k_file)
    chrom = get_chr_string(fas1k_file)
    out_list = []
    # 1. Check that each line is 1000 characters.
    fas1k = open(fas1k_file,'r')
    lines = fas1k.readlines()
    fas1k.close()
    all_line_len = [len(lines[i].strip()) for i in range(0, len(lines)-1)]
    out_list.append(  all([ (all_line_len[i] == 1000) for i in range(0, len(all_line_len)) ]) )
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
        print("something")
    return out_list
    #return "ploidy: " + str(ploidy) + "length: " + str(length) + "chr: " + str(chrom)

def get_chr_string(fas1k_file):
    '''Get chromosome from name of.fas1k files. '''
    # Provide list of acceptable chr names
    chr_list = ['Chr2L', 'Chr2R', 'Chr3L', 'Chr3R', 'Chr4', 'ChrX', 'Yhet', 'ChrX']
    # string match file name against chr names
    bool_list = [chr in fas1k_file for chr in chr_list]
    # Get chr names contained w/i file name
    filtered_list = [i for (i, v) in zip(chr_list, bool_list) if v]
    if ((len(filtered_list) == 0) | (len(filtered_list) > 1)):
        return "Error, can't get chr string from file name: " + fas1k_file
    elif (len(filtered_list) == 1):
        return filtered_list[0]

def arm_to_int(chr_string, lookup_table="int_to_arm.txt"):
    ''' Reference genome has integer scaffold names, translate from arm names to ref names. '''
    # Some of the chr names start with "Chr"- need to strip this out without disturbing names of those 
    #   who don't. Split the string on "Chr"- if starts with "Chr" this outputs a list with 2 entries:
    #   0th empty string & 1st the name. If "Chr" not present, just 1 entry: str of name.
    split_list = chr_string.split("Chr")
    # bool("") = False and bool("anystring") = True
    bool_list = [bool(i) for i in split_list ]
    filtered_list = [i for (i, v) in zip(split_list, bool_list) if v]
    chr_string = filtered_list[0]
    # Compare to lookup table to get integer name from ref genome
    lookup = pd.read_table(lookup_table, header=None, sep=" ")
    int_out = lookup[lookup[1] == chr_string].iloc[0]
    return int_out.iloc[0]
