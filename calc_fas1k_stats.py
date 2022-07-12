#usr/bin/python

import sys
import os
import glob

fstats = open("Fas1k_Stats", 'w')
fstats.write('filename' + '\t' + 'no_genotype' + '\t' + 'heterozygous' + '\t' + 'total_sites' + '\n')
for files in glob.glob("*fas1k"):
    files_vect = files.split('.')
    filename_root = files_vect[0]
    vfile = open(files, 'r')
    n_calls = 0
    het_calls = 0
    total_sites = 0
    for line in vfile:
        line.strip('\n')
        all_line = len(line)
        ncall_line = line.count('N')
        n_calls = n_calls + ncall_line
        # sample single variant heterozygous allele
        het_W_line = line.count('W')
        het_S_line = line.count('S')
        het_M_line = line.count('M')
        het_K_line = line.count('K')
        het_R_line = line.count('R')
        het_Y_line = line.count('Y')
        het_call_line = het_W_line + het_S_line + het_K_line + het_Y_line + het_R_line + het_M_line
        het_calls = het_calls + het_call_line
        total_sites = total_sites + all_line
    # For Y fas1k files have 0 sites that aren't N, so if globbing everything, scripts needs to be able to handle 
    #   cases where total_sites = 0
    def safe_div(x,y):
        if y == 0:
            return 'NA'
        return x / y
    prop_n = str(safe_div(n_calls, total_sites))
    prop_het = str(safe_div(het_calls, total_sites))

    fstats.write(filename_root + '_' + '\t' + str(n_calls) + '\t' + str(het_calls) + '\t' + str(total_sites) + '\n')
