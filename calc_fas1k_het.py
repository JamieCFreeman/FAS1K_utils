#usr/bin/python

import os.path
import pandas as pd

# Function to calc per-site heterozygosity
def calc_seq_het(file_name, skip_line1=False):
    ''' 
    Count heterozygous sites for fas1k file (or fasta file).
    If fasta header is present, set skip_line1=True to ignore header line.
    Output is list containing: file name, count of N calls, count of het calls, total sites (including Ns).
    '''
    # Open input file , initialize empty variables
    vfile = open(file_name, 'r')
    # If fasta file w/ header, skip first line
    if skip_line1==True:
        line1 = vfile.readline().rstrip()
    n_calls = 0
    het_calls = 0
    total_sites = 0
    # Line by line count of het & N sites
    for line in vfile:
        line = line.strip('\n')
        all_line = len(line)
        ncall_line = line.count('N')
        n_calls = n_calls + ncall_line
        # count heterozygous alleles
        het_W_line = line.count('W')
        het_S_line = line.count('S')
        het_M_line = line.count('M')
        het_K_line = line.count('K')
        het_R_line = line.count('R')
        het_Y_line = line.count('Y')
        het_call_line = het_W_line + het_S_line + het_K_line + het_Y_line + het_R_line + het_M_line
        het_calls = het_calls + het_call_line
        total_sites = total_sites + all_line
    
    return [os.path.basename(file_name), n_calls, het_calls, total_sites]

def per_site_het(file_name, skip_lin1=False):
    '''
    Returns per-site heterozygosity for .fas1k file.
    '''
    # Get counts
    counts = calc_seq_het(file_name, skip_line1=False)
    # Get proportion of het sites/ called sites
    het = counts[2] / (counts[3] - counts[1])
    return het

def get_hom_tract_lengths(file_name):
    ''' 
    Get lengths of homozygous tracts from fas1k file- currently ignores Ns, 
    (which probably inflates tract lengths). Consider in-progress. 
    '''
    # Read in file
    files = open(file_name, 'r', newline=None)
    lines = files.readlines()
    files.close()
    # Initialize variables
    tract_current = 0
    tracts_out = []
    het_list = ['Y', 'W', 'M', 'R', 'S', 'K']
    # Process line by line
    for line in lines:
        now = line.strip()
        # Are there any heterozygous sites in this line?
        het_current_line = sum([now.count(i) for i in het_list])
        # If no, add whole line to current tract length (minus Ns)
        if (het_current_line == 0):
            tract_current = tract_current + len(now) - now.count('N')
        # If yes, replace all heterozygous characters with "h", and split the line on all instances
        #   1. the len of the 0th element is added to prev tract length, and that tract is added to output
        #   2. any tracts contained within the line are added to output
        #   3. the last tract is used as the working length to add to next line 
        elif (het_current_line >= 1):
        # replace all het char with h, for simpler string splitting
            for char in het_list:
                now = now.replace(char, 'h')
            # remove Ns so they aren't counted
            now = now.replace('N', '')
            now_split = now.split('h')
            tract_current = tract_current + len( now_split[0] )
            tracts_out.append(tract_current)
            tract_current = 0
            tracts_out.extend( [len(now_split[i]) for i in range(1, len(now_split)-1) ] )
            tract_current = len(now_split[len(now_split)-1])
    # Return list of tract lenghts
    return tracts_out

###############################################################

if __name__ == "__main__":
    import sys
    fas1k_in = sys.argv[1]
    print( per_site_het(fas1k_in, skip_lin1=False) )







