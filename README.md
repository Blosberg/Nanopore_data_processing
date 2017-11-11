# pigx_NP
 Processing nanopore data into human-readable output.

 This script takes fastq data as raw input and outputs a report in html format.
 To run it, first edit the config file: pigx_NP_config.json
and supply the Transcriptome file (with full path) against which you wish to compare the nanopore RNA data, as well as the genome location, genome version and path to input and output folders.

From that point, run the following:

$ snakemake --snakefile pigx_NP.py  

Requires: GLIB_2.14, otherwise returns the following error:
minimap2/minimap2: /lib64/libc.so.6: version `GLIBC_2.14' not found (required by minimap2/minimap2)

the number of mismatches is determined by taking the NM take from minimap2's output and subtracting the "gap" values: thus, mismatches = NM-I-D  , where I and D are taken from the CIGAR string of the read.

