<!-- language: lang-none -->
     _   _                         _
    | \ | |                       (_)
    |  \| | __ _ _ __   ___  _ __  _ _ __   ___ _ __
    | . ` |/ _` | '_ \ / _ \| '_ \| | '_ \ / _ \ '__|
    | |\  | (_| | | | | (_) | |_) | | |_) |  __/ |
    \_| \_/\__,_|_| |_|\___/| .__/|_| .__/ \___|_|
                            | |     | |
                            |_|     |_|              

--- 

This pipeline proocesses nanopore date, aligns reads in bam and 
GRangeslist format, and profiles raw current data along regions of interest. 

Generally, the pipeline is designed to work on raw data directly as it comes
out of the minION device, although a branch for base-called reads is in
progress. Reads are typically binned into "chunks" of ~4000 reads each, and
analysis of the raw current values is performed in conjunction with alignment.
A consolidated final report (the contents of which are under continuing
development) is the ultimate goal.

To run the program, first edit the config file `config.json`, as described
below, and supply a reference file (with full path) populated with regions of
interest against which you wish to analyse the current of the reads (if none is
supplied, this section is simply skipped). You must also provide  the genome
location, genome version and path to input and output folders.

Before running the pipeline, you will need to configure your executables. The
assumption here is that you will be using guix; if you are working from the MDC
in Berlin then everything is already setup for this purpose --otherwise, you
may have to adapt the following steps to your system.  

Navigate to the repository's subfolder `dev/guix` and create your environment from
the manifest file that is provided using the following command:

`{guix} package -p $PWD/.guix-profile -m manifest.scm` 

or, alternatively, just use:

`source generate_env_from_manifest.sh`

Since this file simply contains the above command. Guix will then go about
preparing the dependencies that you will need to run the pipeline in a
bespoke environment which can then be used whenever you want. This will take
some time, but only needs to be done once; go for lunch while it's running. 

Once that's finished, for any future terminal session in which you want to run
the pipeline, you should navigate to the `guix` subfolder and type: 

`source load_env.sh` 

to load your environment. Assuming your config file is setup (see below), you
should then  be able to run the pipeline (see below). 

This pipeline is a work in progress, (note the version number) and feedback is
welcome; if you encounter any problems, please post an issue.

## config file 

=====

| Variable name | description |
| ------------- |:-----------:|
| PATHIN        | string: Required: location of source data. 
| PATHOUT       | string: Required: location of output path to be written to. 
| ReadType    | Type of read under analysis. For now, this should always be "ONT_Direct_RNA", but other strand types will be added in time. 
| samplelist    | list of samples with identifying name (RUN_ID), and the "sampledir" which provides the name of the subdirectory for this dataset within "PATHIN".
| ref           | data regarding the reference --e.g. "Genome_Dir": Path to the reference genome; "Genome_version": name thereof, "RsoI_abspath": path to the directory containing regions of interest. And finally, "RsoI": RDS file containing regions of interest (in a specialized format).
| Execution     | target_out Desired output of the pipeline; e.g. a "bam"; the read data with alignment, in Granges List format "reads_GRL", among various other possibilities. "jobs": max number of rule executions to be run simultaneously. "clustersub": should this task be submitted to an SGE cluster (true/false).

The remaining options can generally be left to their default values in `dev/config_defaults.json`

## Temporary interpreter solution
One last preparatory step must be taken before running the pipeline --it's an unfortunate, kludge solution that will hopefully be removed soon. 

The first line in the file `nanopiper.py` reads `#!...` and what follows  must be replaced by your python interpreter. To find out what that is, type `which python`, and paste the result into the first line in nanopiper.py in place of what is currently there. 

This is an ugly and temporary solution for which I apologize, and it will be removed soon.


## Execution

Once the above have been completed, you should be able to navigate to the folder in which you've saved this repository, and enter: 

`./nanopiper.py -c config_batch_all.json -n`

The `-n` flag signifies a dry-run to confirm that everything is configured in place to run correctly. If no red errors appear when you do you this, then you are ready to run the real thing. 

