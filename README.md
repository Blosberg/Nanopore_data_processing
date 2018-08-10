# Nanopore pipeline:

Processing nanopore data into human-readable output.

This script takes fastq data as raw input and outputs a report in html format.
To run it, first edit the config file: `config.json` and supply the Transcriptome file (with full path) against which you wish to compare the nanopore RNA data, as well as the genome location, genome version and path to input and output folders.

## config file
in the config file there are three variables to define:
target_out can be one of three values: "bam, "GR", or "report", in that order of hierarchy.
"bam": the script will stop when it has generated bam files
"GR": the script will create the bam file and GRanges object with alignments (including current statistics)
"report": both of the previous + the final html report 

From that point, run the following:

`$ snakemake `  

Requires: GLIB_2.14, otherwise returns the following error:
minimap2/minimap2: /lib64/libc.so.6: version `GLIBC_2.14' not found (required by minimap2/minimap2)

the number of mismatches is determined by taking the NM take from minimap2's output and subtracting the "gap" values: thus, mismatches = NM-I-D  , where I and D are taken from the CIGAR string of the read.

Remember that you have to install it from the local machine, so if you get errors like "cannot execute binary", then try `$make clean` and then `$make`.

you also have to `$ unset PYTHONPATH` to make sure that you're doing everything in python 3 (in which snakemake is written). this needs to be done if you have an environment that uses python2 packages (such as ont-tombo).


