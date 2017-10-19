# pigx_NP
 Processing nanopore data into human-readable output.

 This script takes fastq data as raw input and outputs a report in html format.
 To run it, first edit the config file: pigx_NP_config.json
and supply the Transcriptome file (with full path) against which you wish to compare the nanopore RNA data, as well as the genome location, genome version and path to input and output folders.

From that point, run the following:

$ snakemake --snakefile pigx_NP.py  

# Outstanding problems:
This pipeline uses minimap2; the output of this has been observed to yield output peculiar in the following ways:
[1] Reads that do not map to the genome are still carried over to the output with a seqname of "O" (instead of "chr1", "chr2", "chrX", etc.) 
-- A step has been added to filter these entries out to avoid downstream errors.

[2] Additional alignment entries have been observed that do not trace back clearly to any particular read from the input data. These tend to have entries such as the following:
 [1] SA:Z:chr7,128210465,-,261S269M29D640S,1,54;                                         
 [2] SA:Z:chr1,145969268,+,109S166M1D791S,0,24;  

-- These entries have also been filtered out, however their appearance in the first place is somewhat suspicious and warrants caution in the use of this pipeline, until their meaning can be more clearly ascertained.

