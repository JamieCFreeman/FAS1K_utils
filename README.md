
Python functions to handle the Pool lab's .fas1k file spec.

.fas1k files are headerless fasta files with line breaks every 1000 chars to make multifile comparisons easy and to allow indexed acces to files.

##  Functions

**validate_fas1k(** file_name,verbosity=1, ref_fai, soft_mask=False **)**
Checks a fas1k file for conformation to the fas1k specifications. Performs three checks (unless reference genome fasta index file is not provided, then 2).
When run as main, only prints output (failing file name) if tests fail, otherwise completes silently. When sourced, returns either number of tests passing (if verbosity=0) or list of boolean results for each test (when verbosity=1). 
Tests are: 1. Is each line 1000 characters? 2. Are only allowable characters present? set soft_mask to True if the file has been soft masked (contains lowercase nucleotide characters). 3. Does the length of the file match the length of the appropriate reference genome scaffold?
To call from the command line, run: python fas1k_utils.py fas1k_file.fas1k

**calc_fas1k_het(** file_name, skip_line1=False **)** 
