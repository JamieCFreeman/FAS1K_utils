#!/usr/bin/env python3

############################################################
#
# 2024 JCF
#
# Purpose: For PCA genotyping of In(3L)Ok need Zambia diploid sequences 
#	   for clustering. From orginal DGN have
#          many haploid, so will use to make synthetic diploid.
#
############################################################

import random

############################################################

het_dict= {
    "R": {"A", "G"},
    "Y": {"C", "T"},
    "S": {"G", "C"},
    "W": {"A", "T"},
    "K": {"G", "T"},
    "M": {"A", "C"}
    }

# Taken from gt_mat (fix there!)
def expand_set_het(s):
    '''
    For a set of nucleotides, expand heterozygous codes into 
    component nucleotides, and return expanded set of nucleotides
    '''
    het_sites = set(het_dict.keys()).intersection(s)
    if len(het_sites) > 0:
        for x in list(het_sites):
            s = s.union(het_dict[x])
        return s - het_sites
    elif len(het_sites) == 0:
        return s

def check_het(l, verbosity=0):
    '''
    Are any sites heterozygous? Return ploidy of SNP set
    When verbosity is 0, return integer result of ploidy
    When verbosity is 1, return het sites
    '''
    homoz_set = set(('A', 'T', 'C', 'G', 'N'))
    het_set   = set(('Y', 'R', 'W', 'S', 'K', 'M'))
    snp_set   = set(l)
    
    # Does the snp_set intersect with the heterozygous codes?
    if (verbosity == 0):
        if ( len(snp_set.intersection(het_set)) > 0):
            return 2
        elif ( len(snp_set.intersection(het_set)) == 0):
            return 1
    if (verbosity == 1):
        return list( snp_set.intersection(het_set) )


seq1 = "ATGCWRA"
seq2 = "AACTMYN"
poss = ["AWSYMWA", "AWSYMSA", "AWSYMKA", "AWSYMMA", 
        "AWSYWWA", "AWSYWSA", "AWSYWKA", "AWSYWMA" ]


def add_nt(x, y):
    ''' 
    Make heterozygote of two sequences. If no heterozygosity, should be
    determenistic. If heterozygosity, sample alleles per parent.
    '''
    # Expand the nucleotide codes to all included 
    z = expand_set_het(set([x,y]))
    # If there are >= 2 nt, check whether each parental allele is het  
    #   if not, add the allele to out list
    #   if it is, sample one nt and add to out list
    if(len(z)>=2):
        out = []
        for a in [x,y]:
            if(   (check_het(a) == 1) & (a != 'N')):
                out += [a]
            elif( check_het(a) == 2):
                out += random.sample(list(het_dict[a]), 1)
        z = set(out)
    # If there is 1 nt, return that nt
    if(len(z)==1):
        return(x)
    for key, value in het_dict.items():
        if(value == z):
            return(key)

def sample_het(x):
    '''
    Sample a diploid sequence down to one allele
    '''    
    # Expand the nucleotide codes to all included 
    z = expand_set_het(set([x]))
    out = x
    # If heterozygous, sample 1 (returns list w 1 element, so we take the 0th) 
    if(len(z)>=2):
        out = random.sample(list(het_dict[x]), 1)[0]
    return out


