# Nanopore pipeline:

This pipeline proocesses raw minION nanopore data directly to reports and 
figures, as well as bam-file and .RData in GRanges format.

Two forms of raw input are accepted, and can be used to produce a report in
html format.  If the variable `input_data_type`, in the config file, is set to "fastq",
then we assume a single base-called fastq file for each sample, which is then
processed into a report. If this variable is set to "raw_minION", then it is
assumed that raw data directly output from the minION device is being used (in
which case, the fastq data are typically binned into files of ~4000 reads each,
and raw-current data is preserved.) In the latter case, analysis of the raw
current values (before base-calling by Albacore) is performed, although most of
this analysis is not yet included in the final report (the contents of the
report are continually under development).

To run the program, first edit the config file `config.json`, as described
below, and supply the Transcriptome file (with full path) against which you
wish to compare the nanopore reads (if none is supplied, this section is simply
skipped). You must also provide  the genome location, genome version and path
to input and output folders.

The second step is to configure your executables. The assumption here is that
you will be using guix; if you are working from the MDC in Berlin then
everything is already setup for this purpose --otherwise, you may have to adapt
the following steps to your system. 

Navigate to the repository's subfolder `guix` and create your environment from
the manifest file that is provided with the following command:

`/gnu/remote/bin/guixr package -p $PWD/.guix-profile -m manifest.scm` 

Note that this text is the contents of the file
`generate_env_from_manifest.txt`, so you should just be able to type `source
generate_env_from_manifest.txt`. Guix will then go about configuring the
executables that you will need to run the pipeline into a prepared environment
which can then be used whenever you want. This will take some time; go for
lunch.  

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
| input_data_type | string: instructs the pipeline what format of data to expect (either "fastq", or "raw_minION").
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

Requires: GLIB_2.14, otherwise returns the following error:
minimap2/minimap2: /lib64/libc.so.6: version `GLIBC_2.14' not found (required by minimap2/minimap2)

the number of mismatches is determined by taking the NM take from minimap2's output and subtracting the "gap" values: thus, mismatches = NM-I-D  , where I and D are taken from the CIGAR string of the read.

Remember that you have to install it from the local machine, so if you get errors like "cannot execute binary", then try `$make clean` and then `$make`.

you also have to `$ unset PYTHONPATH` to make sure that you're doing everything in python 3 (in which snakemake is written). this needs to be done if you have an environment that uses python2 packages (such as ont-tombo).


