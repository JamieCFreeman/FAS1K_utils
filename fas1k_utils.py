#!/usr/bin/env python3

import os.path
import pandas as pd


# Functions stolen from Chris to get nt of interest from fas1k files.
# Specifically from "/raid10/chris/amplicon_design/primer_gen/inv_primer_gen.py"
# For retrieving a sequence from lines extracted from .fas1k
def subseq_from_lines(start,end,line_len,lines):
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

def extract_strain_subseq(start,end,strain_file):
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
    n_lines = len(lines)
    last_line = len(lines[len(lines)-1].strip())
    return (n_lines*line_len) + last_line

def get_fas1k_ploidy(fas1k_file):
    ''' Compare all characters present in file to those expected based on ploidy, return ploidy as 1 or 2.'''
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
    # If doesn't match either, return unexpected characters
    # (uncovered edge case- fas1k that doesn't have any occurance of 1 nt- 
    #   would only expect this to occur if you were trying to use it on a small
    #   slice of a fas1k file.)  
    else:
        return all_char_set.difference(diploid_set)

def get_fas1k_char(fas1k_file):
    ''' Returns a python set of all characters present in file. '''
    fas1k = open(fas1k_file,'r')
    lines = fas1k.readlines()
    fas1k.close()
    all_char_set = set(lines[0].strip())
    for i in range(len(lines)-1):
        all_char_set.update(set(lines[i].strip()))
    return all_char_set
