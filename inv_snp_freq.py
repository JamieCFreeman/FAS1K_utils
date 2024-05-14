#!/usr/bin/env python3

# #!/usr/bin/env python3
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




######################################################################

def transform_arm(file_path, arm):
	file_arm           = file_path.split("Chr")[1][0:2]
	return re.sub(file_arm, arm, file_path)



def get_name(file_path):
	line_name = file_path.split("/")[len(file_path.split("/"))-1].split("_")[0]
	return line_name


def get_inv_snp_files(dir_path="/home/jamie/FAS1K_utils/inv_fixeddiff_kapun/"):
	# From the directory, get all inv snp files
	dir_all = os.listdir(dir_path)
	subs_files = [ "fixeddiff_kapun.txt" in x for x in dir_all ]
	file_list = list(compress(dir_all, subs_files))
	file_list = [ dir_path + x for x in file_list ]
	return file_list


#def stop_if_fas1k_bad(fas1k_file):
#	try:
#		if sum( validate_fas1k(fas1k_file, ref_fai="/home/jamie/dmel_ref/DmelRef.fasta.fai") ) != 3:
#			raise Exception
#	except:
#		#sys.exit("File names are not unique.")
#		print(".fas1k not according to spec.")
#	else:
#		return "good"


def check_inv_snps(fas1k_file, snp_file, controls=False):
	'''
	For troubleshooting, set controls to True
	'''
	
	# 1. Test fas1k file
	# If the fas1k files are invalid, throw an error
	#stop_if_fas1k_bad(fas1k_file)
	
	# Get current arm
	snp_table     = pd.read_table(snp_file, header=None)
	current_arm   = snp_table.iloc[1][1]
	strain_files  = [ fas1k_file ]
	current_files = [ transform_arm(x, current_arm) for x in strain_files ]
	
	[ get_fas1k_ploidy(x) for x in strain_files ]
	
	if (controls == True):
		# Add inversion control and no inversion control to the run
		inv_control = pd.read_table("inv_control_files.txt", header=None)
		inv_control_file = inv_control[ inv_control[0] == snp_table[0][1] ].iloc[0][1]
		noinv_control_file = transform_arm("/raid10/backups/genepool/DPGP2plus/wrap1kb/Chr3R/FR153N_Chr3R.fas1k", current_arm)
		current_files = [inv_control_file, noinv_control_file] + current_files
	
	pos = snp_table[2]
	for y in current_files:
		short_name = y.split("/")[len(y.split("/"))-1]
		tmp = []
		for x in pos:
			tmp.append( extract_fas1k_subseq(x-1, x, y) )
		snp_table[short_name] = tmp
	
	return snp_table

def count_snp_table(snp_table):
	# Number of Ns
	names        = list(snp_table.columns)[4:snp_table.shape[1]]
	n_counts     = [ sum( snp_table.iloc[:, x] == "N" ) for x in range(4,snp_table.shape[1]) ]
	homoz_match  = [ sum( snp_table.iloc[:, x] == snp_table[3] ) for x in range(4,snp_table.shape[1]) ]
	site_counts  = [ sum( snp_table.iloc[:, x] != "N" ) for x in range(4,snp_table.shape[1]) ]
	inv          = [ snp_table[0][1] for x in range(4,snp_table.shape[1]) ]
	
	# Are there any heterozygous sites? If there are go through and look for het matches, otherwise set het matches to 0
	if (check_het(snp_table)==2):
		het_match = [ count_het_table(snp_table, x) for x in range(4,snp_table.shape[1]) ]
	elif (check_het(snp_table)==1):
		het_match = 0
	
	data_out = pd.DataFrame({'names': names, 'inv': inv, 'N': n_counts, 'inv_match': homoz_match, 'het_match': het_match,'called_sites': site_counts})
	data_out["perc_inv_snps"] = (data_out["inv_match"] + (data_out["het_match"]*0.5) )/ data_out["called_sites"]
	return data_out
	

def check_het(snp_table):
	'''
	Are any sites heterozygous? Return ploidy of SNP set
	'''
	homoz_set = set(('A', 'T', 'C', 'G', 'N'))
	het_set = set(('Y', 'R', 'W', 'S', 'K', 'M'))
	snp_set = set((snp_table.iloc[:, 4].str.cat()))
	
	# Does the snp_set intersect with the heterozygous codes?
	if ( len(snp_set.intersection(het_set)) > 0):
		return 2
	elif ( len(snp_set.intersection(het_set)) == 0):
		return 1

het_dict= {
	"R": ["A", "G"],
	"Y": ["C", "T"],
	"S": ["G", "C"],
	"W": ["A", "T"],
	"K": ["G", "T"],
	"M": ["A", "C"]
	}

def count_het_table(snp_table, x):
	het_sites = [ snp_table.iloc[y, x] in het_dict.keys() for y in range(0,snp_table.shape[0]) ]
	het_idx = list(compress(range(len(het_sites)), het_sites))
	return( sum([ snp_table.iloc[i, 3] in het_dict[snp_table.iloc[i, x]] for i in het_idx ]) )

def check_all_inv_snps(fas1k, file_list=get_inv_snp_files()):
	collect_all = pd.DataFrame()
	
	for z in file_list:
		print(z)
		tmp_out = count_snp_table(check_inv_snps(fas1k, z))
		collect_all = pd.concat([collect_all, tmp_out], ignore_index=True)
	
	return collect_all










