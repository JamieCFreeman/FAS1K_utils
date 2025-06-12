#!/usr/bin/env python3

# importing os.path module
import os.path
import pandas as pd


def calc_seq_diff(x,y):
    ''' Takes two arguments, each a string of nucleotides and outputs a distance score.
    Accepts ACGT and YRWSKM. This doesn't take into account Ns- the compare_fas1k function does that
    (many Ns in fas1k files, so only calling this when no Ns speeds things up significantly.'''
    if(x==y):
        return(0)
    # If x is homozygous and y is heterozygous
    elif(x=='A'):
        if(y=='R' or y=='W' or y=='M'):
            return(1./2)
        else:
            return(1)
    elif(x=='T'):
        if(y=='Y' or y=='W' or y=='K'):
            return(1./2)
        else:
            return(1)
    elif(x=='C'):
        if(y=='Y' or y=='S'or y=='M'):
            return(1./2)
        else:
            return(1)
    elif(x=='G'):
        if(y=='R' or y=='S' or y=='K'):
            return(1./2)
        else:
            return(1)
    # If x is heterozygous and y is homozygous
    elif(x=='W'):
        if(y=='A' or y=='T'):
            return(1./2)
        elif(y=='M' or y=='Y' or y=='R' or y=='K'):
            return(3./4)
        else:
            return(1)
    elif(x=='S'):
        if(y=='C' or y=='G'):
            return(1./2)
        elif(y=='Y' or y=='R' or y=='K' or y=='M'):
            return(3./4)
        else:
            return(1)
    elif(x=='M'):
        if(y=='C' or y=='A'):
            return(1./2)
        elif(y=='Y' or y=='R' or y=='W' or y=='S'):
            return(3./4)
        else:
            return(1)
    elif(x=='K'):
        if(y=='G' or y=='T'):
            return(1./2)
        elif(y=='Y' or y=='R' or y=='W' or y=='S'):
            return(3./4)
        else:
            return(1)
    elif(x=='R'):
        if(y=='A' or y=='G'):
            return(1./2)
        elif(y=='W' or y=='S' or y=='K' or y=='M'):
            return(3./4)
        else:
            return(1)
    elif(x=='Y'):
        if(y=='C' or y=='T'):
            return(1./2)
        elif(y=='W' or y=='S' or y=='K' or y=='M'):
            return(3./4)
        else:
            return(1)
    else:
        return(1)            

def compare_fas1k(path1, path2, comp_func=calc_seq_diff):
    ''' Expects two arguments, both paths to a .fas1k file. For a set of fas1k files, calculate pairwise seq distance
    Set function used to compare as either: calc_seq_diff (default), calc_seq_congruence
    Outputs a list of: [0]=path1, [1]=path2, [2]=number seq differences, [3]=total length seq (no comparison for any 
    nt where either seq has an N), [4]=differences/total length '''
    seqdiffs = 0
    total_length = 0
    with open(path1) as file1, open(path2) as file2:
        for line1, line2 in zip(file1,file2):
          cdiff = sum(comp_func(a,b) for a, b in zip(line1.strip(), line2.strip()) if a != 'N' and b != 'N')
          seqdiffs = seqdiffs + cdiff
          N_Line1 = [i for i, v in enumerate(line1) if v == 'N']
          N_Line2 = [i for i, v in enumerate(line2) if v == 'N']
          Shared_N = list(set(N_Line1) & set(N_Line2))
          total_length = total_length + len(line1) - len(N_Line1) - len(N_Line2) + len(Shared_N)
    # If total_length = 0, return NA
    def safe_div(x,y):
        if y == 0:
            return 'NA'
        return x / y
    prop_diff = safe_div(seqdiffs, total_length)
    out = [os.path.basename(path1), os.path.basename(path2), seqdiffs, total_length, prop_diff]
    return(out)

def read_file_list(file_list):
    ''' Read in list of files (one per line) and return as a list
    '''
    files = open(file_list, 'r', newline=None)
    lines = files.readlines()
    files.close()
    l     = [ x.strip() for x in lines ]
    return l

def compare_fas1k_list(path1, file_list, comp_func=calc_seq_diff):
    ''' Runs compare_fas1k against a list of fas1k files. Outputs pandas data frame.
    Set function used to compare as ethier: calc_seq_diff (default), calc_seq_congruence
    '''
    output_list = [ compare_fas1k(path1, x, comp_func) for x in file_list ]
    df          = pd.DataFrame(output_list, columns = ['File1', 'File2', 'NDIFF', 'NCOMP', 'DISTANCE']) 
    return(df)

def calc_seq_incongruence(x,y):
    ''' Takes two arguments, each a string of nucleotides and outputs a distance score.
    Differs from calc_seq_diff by only adding penalty for alleles that contradict, not just different
    eg: A <-> R(A+G) = 0, W(A+T) <-> M(C+A) = 0.5
    Testing for use in comparing between runs of the same line
    Accepts ACGT and YRWSKM. This doesn't take into account Ns- the compare_fas1k function does that
    (many Ns in fas1k files, so only calling this when no Ns speeds things up significantly.'''
    if(x==y):
        return(0)
    # If x is homozygous and y is heterozygous
    elif(x=='A'):
        if(y=='R' or y=='W' or y=='M'):
            return(0)
        else:
            return(1)
    elif(x=='T'):
        if(y=='Y' or y=='W' or y=='K'):
            return(0)
        else:
            return(1)
    elif(x=='C'):
        if(y=='Y' or y=='S'or y=='M'):
            return(0)
        else:
            return(1)
    elif(x=='G'):
        if(y=='R' or y=='S' or y=='K'):
            return(0)
        else:
            return(1)
    # If x is heterozygous and y is homozygous
    elif(x=='W'):
        if(y=='A' or y=='T'):
            return(0)
        elif(y=='M' or y=='Y' or y=='R' or y=='K'):
            return(1./2)
        else:
            return(1)
    elif(x=='S'):
        if(y=='C' or y=='G'):
            return(0)
        elif(y=='Y' or y=='R' or y=='K' or y=='M'):
            return(1./2)
        else:
            return(1)
    elif(x=='M'):
        if(y=='C' or y=='A'):
            return(0)
        elif(y=='Y' or y=='R' or y=='W' or y=='S'):
            return(1./2)
        else:
            return(1)
    elif(x=='K'):
        if(y=='G' or y=='T'):
            return(0)
        elif(y=='Y' or y=='R' or y=='K' or y=='M'):
            return(1./2)
        else:
            return(1)
    elif(x=='R'):
        if(y=='A' or y=='G'):
            return(0)
        elif(y=='W' or y=='S' or y=='K' or y=='M'):
            return(1./2)
        else:
            return(1)
    elif(x=='Y'):
        if(y=='C' or y=='T'):
            return(0)
        elif(y=='W' or y=='S' or y=='K' or y=='M'):
            return(1./2)
        else:
            return(1)
    else:
        return(1)            


def calc_directional_seq_incong(x,y):
    ''' Takes two arguments, each a string of nucleotides and outputs a distance score.
    Assumes a directional relationship over time between x and y where x is an older seq than y
    Or a directional relationship in inbreeding status x < y
    So in this case A -> R(A+G) = 0.5, but R -> A = 0
    Accepts ACGT and YRWSKM. This doesn't take into account Ns- the compare_fas1k function does that
    (many Ns in fas1k files, so only calling this when no Ns speeds things up significantly.'''
    if(x==y):
        return(0)
    # If x is homozygous and y is heterozygous
    elif(x=='A'):
        if(y=='R' or y=='W' or y=='M'):
            return(1./2)
        else:
            return(1)
    elif(x=='T'):
        if(y=='Y' or y=='W' or y=='K'):
            return(1./2)
        else:
            return(1)
    elif(x=='C'):
        if(y=='Y' or y=='S'or y=='M'):
            return(1./2)
        else:
            return(1)
    elif(x=='G'):
        if(y=='R' or y=='S' or y=='K'):
            return(1./2)
        else:
            return(1)
    # If x is heterozygous and y is homozygous
    elif(x=='W'):
        if(y=='A' or y=='T'):
            return(0)
        elif(y=='M' or y=='Y' or y=='R' or y=='K'):
            return(3./4)
        else:
            return(1)
    elif(x=='S'):
        if(y=='C' or y=='G'):
            return(0)
        elif(y=='Y' or y=='R' or y=='K' or y=='M'):
            return(3./4)
        else:
            return(1)
    elif(x=='M'):
        if(y=='C' or y=='A'):
            return(0)
        elif(y=='Y' or y=='R' or y=='W' or y=='S'):
            return(3./4)
        else:
            return(1)
    elif(x=='K'):
        if(y=='G' or y=='T'):
            return(0)
        elif(y=='Y' or y=='R' or y=='K' or y=='M'):
            return(3./4)
        else:
            return(1)
    elif(x=='R'):
        if(y=='A' or y=='G'):
            return(0)
        elif(y=='W' or y=='S' or y=='K' or y=='M'):
            return(3./4)
        else:
            return(1)
    elif(x=='Y'):
        if(y=='C' or y=='T'):
            return(0)
        elif(y=='W' or y=='S' or y=='K' or y=='M'):
            return(3./4)
        else:
            return(1)
    else:
        return(1)            



#path1 = "/home/jamie/DGN_compatible/EF_testing/round2/fas1k/16Mar22-56_EF10N_S109_round2_Chr2L_diploid.fas1k"
#path2 = "/raid10/octavius_SECOND_6TB_RAID/mac_mini/EF_genomes_Tiago/diploid_calls/gwas/fas1k/EF1N_Chr2L_diploid.fas1k"
#compare_fas1k(path1, path2)
