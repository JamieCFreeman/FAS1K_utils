
#!/usr/bin/env python3

# 2024-06-06 JCF

# Purpose:
# Functions to generate a genotype matrix from a set of fas1k files


###############################################################

from fas1k_utils import *

###############################################################

def flatten_list(l):
	return [item for sublist in l for item in sublist]

def allele_count(l):
	# Each entry string is of length sample + 1 for ref
	# Return a count of non-N alleles
	# First, get non-N alleles
	non_n = set([*l]).difference(set('N'))
	# If no heteozygous sites, allele count is just number of non-N alleles
	if (check_het(l, 0) == 1):
		return len(non_n)
	# If het sites, expand them out to their components and merge with homo sites
	elif (check_het(l, 0) == 2):
		# Next, check for heterozygous sites
		het  = check_het(non_n, 1) 
		exp  = flatten_list([ het_dict[x] for x in het ] )
		
		homoz_set = set(('A', 'T', 'C', 'G', 'N'))
		hom = non_n.intersection(homoz_set)
		
		return len( hom.union( set(exp) ))

def all_n(l):
	# Check whether any samples have called sites, 
	# Returns 0 if no called sites
	m = l[1:len(l)]
	return len( set(m).difference(set('N')) )


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

het_dict= {
	"R": ["A", "G"],
	"Y": ["C", "T"],
	"S": ["G", "C"],
	"W": ["A", "T"],
	"K": ["G", "T"],
	"M": ["A", "C"]
	}

def score_invariant(r, x):
	# If site is invariant, there are two possible calls:
	# N
	if (x == 'N'):
		return '9'
	# homozygous ref (0)
	elif (x == r):
		return '0'

def score_biallelic(r, x):
	if   ( x == 'N'):
		return '9'
	# homozygous ref (0)
	elif ( x == r):
		return '0'
	# het (1)
	elif ( check_het(x) == 2):
		return '1'
	# homozy alt (2)
	elif ( check_het(x) == 1):
		return '2'

def gt_site(l):
	# First take each string, and split into indiv char
	v = [*l]
	# Remove ref allele from front
	r = v.pop(0)
	s = ''
	if ( allele_count(l) == 1 ):
		for i in v:
			s = s + score_invariant(r, i)
	if ( allele_count(l) == 2 ):
		for i in v:
			s = s + score_biallelic(r, i)
	return s

def convert_site(l):
	# If all sites are n, return a 9 for each sample
	if ( all_n(l) == 0 ):
		s = '9' * ( len(l) -1 )
		return s
	# Now consider completely homozygous sites
	elif ( allele_count(l) <= 2):
		s = gt_site(l)
		return s
	elif ( allele_count(l) >= 3):
		s = '3' * ( len(l) -1 )
		return s

###############################################################

