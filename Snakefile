#!/usr/bin/env python3.5
import os, sys, json, csv, yaml
import argparse

# import IPython;

#-----------------------------------------------------------------
#--- define settings for config file -----------------------------
config_defaults="dev/config_defaults.json"
config = yaml.safe_load(open(config_defaults, 'r'))

config_userin="config.json"
update_config( config, yaml.safe_load( open( config_userin, 'r')))

config_log="configlog_out.json"

# look in the path provided for each sample and build a list of chunks in each case.
for sample in config["samplelist"]:
    # set path within this sample's data set to look for chunks
    DATPATH = os.path.join(config["PATHIN"], config["samplelist"][sample]["subdir"], "fast5", "pass" )

    # define a list of these subentries:
    config["samplelist"][sample]["chunkdirlist"] = [entry.name for entry in os.scandir( DATPATH ) if entry.is_dir()]

# store a log of the config file we just built:
config_log="configlog_out.json"

with open(config_log, 'w') as outfile:
    dumps = json.dumps(config,
                       indent=4, sort_keys=True,
                       separators=(",",": "), ensure_ascii=True)
    outfile.write(dumps)


# ------------------------------------------------------
#--- Dependencies:

include: os.path.join( config["scripts"]["script_folder"], config["scripts"]["pyfunc_defs"] )

MM2        = config["progs"]["minimap"]
SAMTOOLS   = config["progs"]["SAMTOOLS"]
nanopolish = config["progs"]["nanopolish"]

RefTranscriptome   = config["ref"]["Transcriptome"]
GENOME_VERSION     = config["ref"]["Genome_version"]
RmdReportScript    = os.path.join(config["scripts"]["script_folder"],"final_report","Nanopore_report.Rmd")
input_data_type    = config["input_data_type"]

R_tables2GR_main     = os.path.join( config["scripts"]["script_folder"], config["scripts"]["Rmain_table2GRconv"] )
R_tables2GR_funcs    = os.path.join( config[ "scripts"]["script_folder"], config[ "scripts"]["Rfuncs_table2GRconv"] )

R_flattenreads_main  = os.path.join( config["scripts"]["script_folder"], config["scripts"]["Rmain_flattenreads"] )
R_flattenreads_funcs = os.path.join( config["scripts"]["script_folder"], config["scripts"]["Rfuncs_flattenreads"] )

R_build_histlist_main   = os.path.join( config["scripts"]["script_folder"], config["scripts"]["Rmain_build_histlist"] )
#------------------------------------------------------
#--- Define output directories

SUBDIR_SYMLINKS           = "01_symlinks_fastqindex_chunks/"
os.makedirs( os.path.join( config["PATHOUT"], config["samplelist"][sample]["subdir"], SUBDIR_SYMLINKS),  exist_ok=True)
# -- all the rest will be created by snakemake automatically

SUBDIR_ALIGNED_MINIMAP    = "02_MM_aligned_chunks/"
SUBDIR_FILTERED_MINIMAP   = "03_MM_filtered_chunks/"
SUBDIR_SORTED_MINIMAPPED  = "04_MM_sortedbam/"
SUBDIR_EVENTALIGN         = "05_eventalign/"
SUBDIR_GR                 = "06_GRobjects"
SUBDIR_REPORT             = "Final_report/"
DIR_REFGENOME             = config['ref']['Genome_DIR']

#------------------------------------------------------
#--- check that the pipeline can be executed:

# check for write access to refgenome dir
if ( not os.access(DIR_REFGENOME, os.W_OK) ):
   print("Write access to refgenome folder is denied. Checking if necessary indexing files already exist: ... ")

   if( not os.path.isfile(DIR_REFGENOME , config['ref']['Genome_version']+ ".mmi")):
      bail("minimap index files not found, and cannot be created. Aborting")

   else:
      print("Refgenome index files are present. Continuing... ")

#---  Create symbolic links to PATHIN so that indexing/etc can be performed in
# written pathout

if ( input_data_type == "raw_minION"):
   # the snakemake rules defined in the following included script assume
   # that the minION output has been "chunked" into sequential files, and
   # will have to be reassembled.
   include   : os.path.join( config["scripts"]["script_folder"], config["scripts"]["rules_chunks"] )

   for sample in config["samplelist"]:

       # get subdir for this sample:
       samplePATH = os.path.join(config["PATHIN"], config["samplelist"][sample]["subdir"] )

       # linke to each chunk directly:
       for linkindex in config["samplelist"][sample]["chunkdirlist"]:
          linkname = sample + "_" + str(linkindex) + config['samplelist'][sample]["fastq_suffix"]


          source   = getPathCase( samplePATH, "fastq", "pass",config["samplelist"][sample]["fastq_prefix"] + str(linkindex) +  config["samplelist"][sample]["fastq_suffix"], "raw_minION")

          makelink(  source, os.path.join( config["PATHOUT"], config["samplelist"][sample]["subdir"], SUBDIR_SYMLINKS, linkname) )


elif( input_data_type == "fastq"):
   # Do some other stuff
   include   : os.path.join( config["scripts"]["script_folder"], config["scripts"]["rules_wholefastq"] )

   for sample in config["samplelist"]:

       # get subdir for this sample:
       samplePATH = os.path.join(config["PATHIN"], config["samplelist"][sample]["subdir"] )

       linkname = sample + config['samplelist'][sample]["fastq_suffix"]
       makelink( os.path.join(samplePATH, config["samplelist"][sample]["fastq_prefix"] + config["samplelist"][sample]["fastq_suffix"] ), os.path.join( config["PATHOUT"], config["samplelist"][sample]["subdir"], SUBDIR_SYMLINKS, linkname))

else:
   print("Unrecognized input data format. Terminating.")
   exit(1)


#------------------------------------------------------
#---  Define output (target) files:

if ( config["target_out"] == "report" ):
   OUTPUT_FILES=  [ os.path.join( config["PATHOUT"], config["samplelist"][sample]["subdir"], SUBDIR_REPORT, ""+ sample +"_report.html") for sample in config["samplelist"] ]
elif ( config["target_out"] == "histlist" ):
   OUTPUT_FILES=  [
                  os.path.join( config["PATHOUT"], config["samplelist"][sample]["subdir"], SUBDIR_GR, sample+"_kmer_histlist.rds")  for sample in config["samplelist"]
                  ]
elif ( config["target_out"] == "flatreads_GRL" ):
   OUTPUT_FILES=  [
                  os.path.join( config["PATHOUT"], config["samplelist"][sample]["subdir"], SUBDIR_GR, sample+"_reads_flat_GRL.rds")  for sample in config["samplelist"]
                  ]
elif ( config["target_out"] == "reads_GRL" ):
   OUTPUT_FILES=  [
                  os.path.join( config["PATHOUT"], config["samplelist"][sample]["subdir"], SUBDIR_GR, sample+"_reads_GRL.rds")  for sample in config["samplelist"]
                  ]
elif ( config["target_out"] == "bam" ):
   OUTPUT_FILES=  [
                  os.path.join( config["PATHOUT"], config["samplelist"][sample]["subdir"], SUBDIR_SORTED_MINIMAPPED, "run_"+ sample+".sorted.bam") for sample in config["samplelist"]
                  ]
else:
   print("Unrecognized target output file format: ", config["target_out"], " ... Terminating.")
   exit(1)

#---  DEBUGGING:
#------------------------------------------------------
# print("input_data_type = " + config["input_data_type"])
# print("target out = " + config["target_out"])
print("OUTPUT_FILES=")
for x in OUTPUT_FILES:
  print(x)
print("\n finished outputting output files \n\n ")
# IPython.embed()
#
# ========================================================================
#
#   BEGIN RULES
#
# ========================================================================

rule all:
    input:
        [ OUTPUT_FILES ]

#------------------------------------------------------

# build the final output report in html format
rule make_report:
    input:
        aligned_reads_bam = os.path.join( config["PATHOUT"], config["samplelist"][sample]["subdir"], SUBDIR_SORTED_MINIMAPPED, "run_{sample}.sorted.bam"),
        transcriptome     = RefTranscriptome,
        GRobj             = os.path.join( config["PATHOUT"], config["samplelist"][sample]["subdir"], SUBDIR_GR, "{sample}_reads_GRL.rds")
    output:
        os.path.join( config["PATHOUT"], config["samplelist"][sample]["subdir"], SUBDIR_REPORT, "{sample}_report.html")
    params:
        " readcov_THRESH = 10;   ",
        " yplotmax = 10000; "
    log:
        logfile = os.path.join( config["PATHOUT"], config["samplelist"][sample]["subdir"], SUBDIR_REPORT, "final_report_{sample}.log")
    message:
        fmt("producing final report")
    shell:
        " Rscript -e  '{params} "
        " fin_readalignment_bam = \"{input.aligned_reads_bam}\"; "
        " fin_Transcript        = \"{input.transcriptome}\";"
        " fin_GRobj             = \"{input.GRobj}\";"
        " Genome_version        =\"{GENOME_VERSION}\"; "
        " rmarkdown::render(\"{RmdReportScript}\", output_file = \"{output}\" ) ' "

#------------------------------------------------------

# Align the events to the reference genome.
# The wildcard "chunk" can simply be "full", in cases
# where there are no chunks
rule np_event_align:
    input:
        sortedbam             = os.path.join( config["PATHOUT"], config["samplelist"][sample]["subdir"], SUBDIR_SORTED_MINIMAPPED, "bam_chunks", "run_{sample}_{chunk}.sorted.bam"),
        NOTCALLED_indexedbam  = os.path.join( config["PATHOUT"], config["samplelist"][sample]["subdir"], SUBDIR_SORTED_MINIMAPPED, "bam_chunks", "run_{sample}_{chunk}.sorted.bam.bai"),
        fastq_file            = os.path.join( config["PATHOUT"], config["samplelist"][sample]["subdir"], SUBDIR_SYMLINKS,  "{sample}_{chunk}" +config["samplelist"][sample]["fastq_suffix"]),
        fastq_npi             = os.path.join( config["PATHOUT"], config["samplelist"][sample]["subdir"], SUBDIR_SYMLINKS,  "{sample}_{chunk}" + config["samplelist"][sample]["fastq_suffix"] + ".index"),
        refgenome_fasta       = os.path.join( DIR_REFGENOME, config['ref']['Genome_version']+ ".fa" ),
        NOTCALLED_bwt         = os.path.join( DIR_REFGENOME, config['ref']['Genome_version']+ ".fa.bwt"),
        NOTCALLED_pac         = os.path.join( DIR_REFGENOME, config['ref']['Genome_version']+ ".fa.pac")
    output:
        Ealigned         = os.path.join( config["PATHOUT"], config["samplelist"][sample]["subdir"], SUBDIR_EVENTALIGN, "csv_chunks", 'Ealign_{sample}_{chunk}.csv' )
    log:
        logfile  = os.path.join( config["PATHOUT"], config["samplelist"][sample]["subdir"], SUBDIR_EVENTALIGN, "csv_chunks", 'Ealign_{sample}_{chunk}.log')
    message: """---- align events from sample {wildcards.sample}, chunk {wildcards.chunk} to the genome ----"""
    shell:
        " {nanopolish} eventalign --reads {input.fastq_file} --bam {input.sortedbam} --genome {input.refgenome_fasta} --scale-events  > {output}  2> {log.logfile} "

#------------------------------------------------------
# rule quickcheck: (TODO)
#------------------------------------------------------

# Index the sorted bam file with samtools
rule index_sortedbam:
    input:
        sortedbam  = os.path.join( config["PATHOUT"], config["samplelist"][sample]["subdir"], SUBDIR_SORTED_MINIMAPPED, "bam_chunks", "run_{sample}_{chunk}.sorted.bam")
    output:
        indexedbam = os.path.join( config["PATHOUT"], config["samplelist"][sample]["subdir"], SUBDIR_SORTED_MINIMAPPED, "bam_chunks", "run_{sample}_{chunk}.sorted.bam.bai")
    log:
        logfile    = os.path.join( config["PATHOUT"], config["samplelist"][sample]["subdir"], SUBDIR_SORTED_MINIMAPPED, "bam_chunks", 'run_{sample}_{chunk}_samtoolsindex.log')
    message: """---- index the bam files for {wildcards.sample} chunk {wildcards.chunk} ----"""
    shell:
        " {SAMTOOLS} index  {input.sortedbam}  2> {log.logfile} "

#------------------------------------------------------
# Index the reads and the fast5 files for nanopolish
rule np_index:
    input:
        fast5_folder = lambda wc: getPathCase( os.path.join( config["PATHIN"], config["samplelist"][sample]["subdir"] ), 'fast5', 'pass', wc.chunk, input_data_type ),
        fastq_file   = lambda wc: os.path.join( config["PATHOUT"], config["samplelist"][sample]["subdir"], SUBDIR_SYMLINKS, wc.sample + "_" + str(wc.chunk) + config['samplelist'][wc.sample]["fastq_suffix"] )
    output:
        npi    = os.path.join( config["PATHOUT"], config["samplelist"][sample]["subdir"], SUBDIR_SYMLINKS, "{sample}_{chunk}"+config["samplelist"][sample]["fastq_suffix"]+".index"  ),
        fai    = os.path.join( config["PATHOUT"], config["samplelist"][sample]["subdir"], SUBDIR_SYMLINKS, "{sample}_{chunk}"+config["samplelist"][sample]["fastq_suffix"]+".index.fai" ),
        gzi    = os.path.join( config["PATHOUT"], config["samplelist"][sample]["subdir"], SUBDIR_SYMLINKS, "{sample}_{chunk}"+config["samplelist"][sample]["fastq_suffix"]+".index.gzi" ),
        readdb = os.path.join( config["PATHOUT"], config["samplelist"][sample]["subdir"], SUBDIR_SYMLINKS, "{sample}_{chunk}"+config["samplelist"][sample]["fastq_suffix"]+".index.readdb")
    params:
        options    = " index -d "
    log:
        logfile  = os.path.join( config["PATHOUT"], config["samplelist"][sample]["subdir"], SUBDIR_SYMLINKS,  "{sample}_{chunk}_npi.log" )
    message: """---- index the reads from chunk {wildcards.chunk} against the fast5 files from the same. ----"""
    shell:
        " nice -19 {nanopolish} {params.options} {input.fast5_folder} {input.fastq_file} 2> {log.logfile} "

