
Python functions to handle the Pool lab's .fas1k file spec.

.fas1k files are headerless fasta files with line breaks every 1000 chars to make multifile comparisons easy and to allow indexed acces to files.
The chromosome arm is indicated in the file name. In the Nexus release, heterozygous regions were soft-masked (lowercase instead of uppercase).


## Scripts for command-line use

##  Functions


**extract_fas1k_subseq**(start, end, fas1k_file)
Returns a specified subsequence of a fas1k file in 0-based coordinates (so if you are using reference genome coordinates make sure to subtract 1!). 

**get_chr_string**(fas1k_file)
Extracts the chromosome identifier from the file name of a fas1k file. Should return one of ['Chr2L', 'Chr2R', 'Chr3L', 'Chr3R', 'Chr4', 'ChrX', 'Yhet', 'mtDNA']

**arm_to_int**(chr_string, lookup_table=arm_lookup)
For a chromsome arm string formatted as in fas1k file names, return the integer scaffold in the Drosophila reference genome.

**get_fas1k_length**(fas1k_file)
Returns the nucleotide length of a fas1k file. Does not check length of lines.

**get_fas1k_ploidy**(fas1k_file, soft_mask=False)
Returns the ploidy (1 or 2) of a fas1k file. Checks that the character set includes only legal characters 
(nucleotides and ambiguity codes) and will return error if unexpected characters are present.

**validate_fas1k(** file_name,verbosity=1, ref_fai, soft_mask=False **)**
Checks a fas1k file for conformation to the fas1k specifications. Performs three checks (unless reference genome fasta index file is not provided, then 2).
When run as main, only prints output (failing file name) if tests fail, otherwise completes silently. When sourced, returns either number of tests passing (if verbosity=0) or list of boolean results for each test (when verbosity=1). 
Tests are: 1. Is each line 1000 characters? 2. Are only allowable characters present? set soft_mask to True if the file has been soft masked (contains lowercase nucleotide characters). 3. Does the length of the file match the length of the appropriate reference genome scaffold?
To call from the command line, run: python fas1k_utils.py fas1k_file.fas1k

**calc_fas1k_het(** file_name, skip_line1=False **)** 
