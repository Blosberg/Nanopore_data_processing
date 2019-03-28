 _   _                         _
| \ | |                       (_)
|  \| | __ _ _ __   ___  _ __  _ _ __   ___ _ __
| . ` |/ _` | '_ \ / _ \| '_ \| | '_ \ / _ \ '__|
| |\  | (_| | | | | (_) | |_) | | |_) |  __/ |
\_| \_/\__,_|_| |_|\___/| .__/|_| .__/ \___|_|
                        | |     | |
                        |_|     |_|              

# ==============================================

This pipeline proocesses nanopore date into reports, and aligned reads in bam or 
GRangeslist format. 

Two forms of raw input are accepted, single base-called fastq files, and raw 
data directly output from the minION device. In the latter case, reads are 
typically binned into files of ~4000 reads each, and raw-current data is 
preserved; analysis of the raw current values is performed. Incorporation of 
this analysis is yet to be included in the final report (the contents of
which are under continuing development).

To run the program, first edit the config file `config.json`, as described
below, and supply the Transcriptome file (with full path) against which you
wish to compare the nanopore reads (if none is supplied, this section is simply
skipped). You must also provide  the genome location, genome version and path
to input and output folders.

The second step is to configure your executables. The assumption here is that
you will be using guix; if you are working from the MDC in Berlin then
everything is already setup for this purpose --otherwise, you may have to adapt
the following steps to your system.  

Navigate to the repository's subfolder `dev/guix` and create your environment from
the manifest file that is provided using the following command:

`{guix} package -p $PWD/.guix-profile -m manifest.scm` 

or, alternatively, just use:

`source generate_env_from_manifest.sh`

Since this file simply contains the above command. Guix will then go about
configuring the executables that you will need to run the pipeline into a
prepared environment which can then be used whenever you want. This will take
some time, but only needs to be done once; go for lunch.  

Once that's finished, for any future terminal session in which you want to run
the pipeline, you should navigate to the `guix` subfolder and type: 

`source load_env.sh` 

to load your environment. Assuming your config file is setup (see below), you
should then simply be able to type `Snakemake` and everything will run. 

This pipeline is a work in progress, but feedback is welcome; if you encounter
any problems, please email Brendan.Osberg@mdc-berlin.de

## config file 

=====

| Variable name | description |
| ------------- |:-----------:|
| PATHIN        | string: Required: location of source data. 
| PATHOUT       | string: Required: location of output path to be written to. 
| samplelist    | list of samples with identifying name (RUN_ID), number of bins (MAXSAMPLEi), and hashID (fastq_prefix), and file extension (fastq_suffix).
| target_out | desired output format; the pipeline will proceed until it has accomplished the desired output; either a "bam", an RData file with "GR"anges, or a full html "report" (the latter being the default).
| ref        | Various data related to the reference genome: "Transcriptome" -absolute path to the file with transcripts to compare to; "Genome_DIR": absolute path leading to the reference genome (fasta); and "Genome_version" (string), the name of the genome
| scripts   | paths to various executables needed for the pipeline ( it is recommended to not change these.)
| execution: nice | the degree to which your process is made "nice" to reduce congestion on your computation nodes (integer, max. 19.)


=====

in the config file there are three variables to define:
target_out can be one of three values: "bam, "GR", or "report", in that order of hierarchy.
"bam": the script will stop when it has generated bam files
"GR": the script will create the bam file and GRanges object with alignments (including current statistics)
"report": both of the previous + the final html report 

From that point, run the following:

`$ snakemake `  
