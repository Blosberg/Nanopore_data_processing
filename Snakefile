#!/usr/bin/env python3.5
import os

# import IPython; 

# set config file
configfile: "./config.json"
include   : os.path.join( config["scripts"]["script_folder"], config["scripts"]["pyfunc_defs"] )

#------------------------------------------------------
# --- Dependencies:

MM2        = config["progs"]["minimap"]
SAMTOOLS   = config["progs"]["SAMTOOLS"]
nanopolish = config["progs"]["nanopolish"]

RefTranscriptome = config["ref"]["Transcriptome"]
GENOME_VERSION   = config["ref"]["Genome_version"]
RmdReportScript  = os.path.join(config["scripts"]["script_folder"],"final_report","Nanopore_report.Rmd")
input_data_type           = config["input_data_type"]

#------------------------------------------------------
# --- define output directories

DIR_SYMLINKS           = config["PATHOUT"]+"01_symlinks/"
os.makedirs( DIR_SYMLINKS,  exist_ok=True)
# -- all the rest will be created by snakemake automatically

DIR_ALIGNED_MINIMAP    = config["PATHOUT"]+"02_MM_aligned/"
DIR_FILTERED_MINIMAP   = config["PATHOUT"]+"03_MM_filtered/"
DIR_SORTED_MINIMAPPED  = config["PATHOUT"]+"04_MM_sortedbam/"
DIR_EVENTALIGN         = config["PATHOUT"]+"05_eventalign/"
DIR_GR                 = config["PATHOUT"]+"06_GRobjects"
DIR_REPORT             = config["PATHOUT"]+"Final_report/"
DIR_REFGENOME          = config['ref']['Genome_DIR']

#------------------------------------------------------
# check that the pipeline can be executed: 

if ( not os.access(DIR_REFGENOME, os.W_OK) ):
   print("Write access to refgenome folder is denied. Checking if necessary indexing files already exist: ... ")

   if( not os.path.isfile(os.path.join(DIR_REFGENOME , config['ref']['Genome_version']+ ".mmi")) ):  
      bail("minimap index files not found, and cannot be created. Aborting")

   else:
      print("Refgenome index files are present. Continuing... ")

# Create symbolic links to PATHIN so that indexing/etc can be performed in
# written pathout 

if ( input_data_type == "raw_minION"): 
   # the snakemake rules defined in the following included script assume
   # that the minION output has been "chunked" into sequential files, and
   # will have to be reassembled. 
   include   : os.path.join( config["scripts"]["script_folder"], config["scripts"]["rules_chunks"] )

   for sample in config["samplelist"]:
       for linkindex in range(0, config['samplelist'][sample]["MAXSAMPLEi"] + 1):
          linkname = config['samplelist'][sample]["RUN_ID"] + "_" + str(linkindex) + config['samplelist'][sample]["fastq_suffix"]

          source   = getPathCase( config["PATHIN"], "fastq", "pass",config["samplelist"][sample]["fastq_prefix"] + str(linkindex) +  config["samplelist"][sample]["fastq_suffix"], "raw_minION") 

          makelink(  source, os.path.join( DIR_SYMLINKS, linkname) )

 
elif( input_data_type == "fastq"):
   # Do some other stuff
   include   : os.path.join( config["scripts"]["script_folder"], config["scripts"]["rules_wholefastq"] )

   for sample in config["samplelist"]:

       linkname = config['samplelist'][sample]["RUN_ID"] + config['samplelist'][sample]["fastq_suffix"]
       makelink( os.path.join(config["PATHIN"], config["samplelist"][sample]["fastq_prefix"] + config["samplelist"][sample]["fastq_suffix"] ),
                 os.path.join( DIR_SYMLINKS, linkname))

else: 
   print("Unrecognized input data format. Terminating.")
   exit(1)

#------------------------------------------------------
# Define output (target) files:

if ( config["target_out"] == "report" ):
   OUTPUT_FILES=  [
                  os.path.join( DIR_REPORT, ""+config["samplelist"][sample]["RUN_ID"]+"_report.html") for sample in config["samplelist"]
                  ]
elif ( config["target_out"] == "GR" ):
   OUTPUT_FILES=  [
                  os.path.join( DIR_GR, config["samplelist"][sample]["RUN_ID"]+"_GR.RData")  for sample in config["samplelist"]
                  ]
elif ( config["target_out"] == "bam" ):
   OUTPUT_FILES=  [
                  os.path.join( DIR_SORTED_MINIMAPPED, "run_"+ config["samplelist"][sample]["RUN_ID"]+".sorted.bam") for sample in config["samplelist"]
                  ]
else:
   print("Unrecognized target output file format: ", config["target_out"], " ... Terminating.")
   exit(1)

# DEBUGGING:
#------------------------------------------------------
# print("input_data_type = " + config["input_data_type"])
# print("target out = " + config["target_out"])
# print("OUTPUT_FILES=")
# for x in OUTPUT_FILES: 
#   print(x)
# print("\n finished outputting output files \n\n ")
# IPython.embed()
# 
#=========================================================================
#
#   BEGIN RULES    
#
#=========================================================================

rule all:
    input:
        [ OUTPUT_FILES ]

#------------------------------------------------------

# build the final output report in html format
rule make_report:
    input:
        aligned_reads = os.path.join( DIR_SORTED_MINIMAPPED, "run_{sample}.sorted.bam"),
        transcriptome = RefTranscriptome,
        GRobj         = os.path.join( DIR_GR, "{sample}_GR.RData")
    output:
        os.path.join( DIR_REPORT, "{sample}_report.html")
    params:
        " readcov_THRESH = 10;   ",
        " yplotmax = 10000; "
    log:
        logfile = os.path.join( DIR_REPORT, "final_report_{sample}.log")
    message: """--- producing final report."""

    shell: 
        " Rscript -e  '{params} "
        " fin_readalignment = \"{input.aligned_reads}\"; "
        " fin_Transcript = \"{input.transcriptome}\";"
        " fin_GRobj      = \"{input.GRobj}\";"
        " Genome_version=\"{GENOME_VERSION}\"; " 
        " rmarkdown::render(\"{RmdReportScript}\", output_file = \"{output}\" ) ' "  

#------------------------------------------------------

# Align the events to the reference genome. 
# The wildcard "chunk" can simply be "full", in cases 
# where there are no chunks
rule np_event_align:
    input:
        sortedbam             = os.path.join( DIR_SORTED_MINIMAPPED, "bam_chunks", "run_{sample}_{chunk}.sorted.bam"),
        NOTCALLED_indexedbam  = os.path.join( DIR_SORTED_MINIMAPPED, "bam_chunks", "run_{sample}_{chunk}.sorted.bam.bai"),
        fastq_file            = os.path.join( DIR_SYMLINKS,  "{sample}_{chunk}" +config["samplelist"][sample]["fastq_suffix"]),
        fastq_npi             = os.path.join( DIR_SYMLINKS,  "{sample}_{chunk}" + config["samplelist"][sample]["fastq_suffix"] + ".index"),
        refgenome_fasta       = os.path.join( DIR_REFGENOME, config['ref']['Genome_version']+ ".fa" ),
        NOTCALLED_bwt         = os.path.join( DIR_REFGENOME, config['ref']['Genome_version']+ ".fa.bwt"),
        NOTCALLED_pac         = os.path.join( DIR_REFGENOME, config['ref']['Genome_version']+ ".fa.pac")
    output:
        Ealigned         = os.path.join( DIR_EVENTALIGN, "csv_chunks", 'Ealign_{sample}_{chunk}.csv' )
    log:
        logfile  = os.path.join( DIR_EVENTALIGN, "csv_chunks", 'Ealign_{sample}_{chunk}.log')
    message: """---- align events from sample {wildcards.sample}, chunk {wildcards.chunk} to the genome ----"""
    shell:
        " {nanopolish} eventalign --reads {input.fastq_file} --bam {input.sortedbam} --genome {input.refgenome_fasta} --scale-events  > {output}  2> {log.logfile} "

#------------------------------------------------------
# rule quickcheck: (TODO)
#------------------------------------------------------

# Index the sorted bam file with samtools
rule index_sortedbam:
    input:
        sortedbam  = os.path.join( DIR_SORTED_MINIMAPPED, "bam_chunks", "run_{sample}_{chunk}.sorted.bam")
    output:
        indexedbam = os.path.join( DIR_SORTED_MINIMAPPED, "bam_chunks", "run_{sample}_{chunk}.sorted.bam.bai")
    log:
        logfile    = os.path.join( DIR_SORTED_MINIMAPPED, "bam_chunks", 'run_{sample}_{chunk}_samtoolsindex.log')
    message: """---- index the bam files for {wildcards.sample} chunk {wildcards.chunk} ----"""
    shell:
        " {SAMTOOLS} index  {input.sortedbam}  2> {log.logfile} "

#------------------------------------------------------

# Index the reads and the fast5 files for nanopolish
rule np_index:
    input:
        fast5_folder = lambda wc: getPathCase( config['PATHIN'], 'fast5', 'pass', wc.chunk, input_data_type ),
        fastq_file   = lambda wc: os.path.join( DIR_SYMLINKS, config['samplelist'][wc.sample]["RUN_ID"] + "_" + str(wc.chunk) + config['samplelist'][wc.sample]["fastq_suffix"] ) 
    output:
        npi    = os.path.join( DIR_SYMLINKS, "{sample}_{chunk}"+config["samplelist"][sample]["fastq_suffix"]+".index"  ), 
        fai    = os.path.join( DIR_SYMLINKS, "{sample}_{chunk}"+config["samplelist"][sample]["fastq_suffix"]+".index.fai" ),
        gzi    = os.path.join( DIR_SYMLINKS, "{sample}_{chunk}"+config["samplelist"][sample]["fastq_suffix"]+".index.gzi" ),  
        readdb = os.path.join( DIR_SYMLINKS, "{sample}_{chunk}"+config["samplelist"][sample]["fastq_suffix"]+".index.readdb")
    params:
        options    = " index -d "
    log:
        logfile  = os.path.join( DIR_SYMLINKS,  "{sample}_{chunk}_npi.log" )
    message: """---- index the reads from chunk {wildcards.chunk} against the fast5 files from the same. ----"""
    shell:
        " nice -19 {nanopolish} {params.options} {input.fast5_folder} {input.fastq_file} 2> {log.logfile} "

# The following commented-out rules were used for the BWA branch 
# -----------------------------------------------------
# Create BWA-indexed version of reference genome for fast
#  alignment with bwa later:
# rule bwa_index:
#     input:
#         refgenome_fasta  = os.path.join(DIR_REFGENOME , config['ref']['Genome_version']+ ".fa" )
#     output:
#         bwt      = os.path.join(DIR_REFGENOME , config['ref']['Genome_version']+ ".fa.bwt"),
#         pac      = os.path.join(DIR_REFGENOME , config['ref']['Genome_version']+ ".fa.pac")
#     params:
#         options  = " index  "
#     log:
#         logfile  = os.path.join( DIR_REFGENOME, config['ref']['Genome_version'], "_bwa_indexing.log")
#     message: 
#         """---- creating bwa index of the reference genome. ----"""
#     shell:
#         "{BWA} {params.options}  {input} > {log.logfile}"
#------------------------------------------------------
# Align the reads to the reference with BWA
# rule align_bwa_mem_ont2d:
#     input:
#         refg_fasta = os.path.join(DIR_REFGENOME , config['ref']['Genome_version']+ ".fa" ),
#         refg_bwt   = os.path.join(DIR_REFGENOME , config['ref']['Genome_version']+ ".fa.bwt"),
#         reads      = os.path.join( DIR_SYMLINKS,  "{sample}_{chunk}" + config["samplelist"][sample]["fastq_suffix"] ),
#         npi        = os.path.join( DIR_SYMLINKS,  "{sample}_{chunk}" + config["samplelist"][sample]["fastq_suffix"] + ".index")
#     output:
#         sortedbam  = os.path.join( DIR_SORTED_ALIGNED_BWA, "chunks", "fastq_run_{sample}_{chunk}.bwaligned.sorted.bam")
#     params:
#         options    = " mem -x ont2d ",
#         tempfile   = os.path.join( DIR_SORTED_ALIGNED_BWA, "chunks", "fastq_run_{sample}_{chunk}.bwaligniment.log")
#     log:
#         logfile  = os.path.join( DIR_SORTED_ALIGNED_BWA, "chunks", 'alignment_{sample}_{chunk}_bwaMemOnt2d.log')
#     message: """---- Align the reads from chunk {wildcards.chunk} to the reference ----"""
#     shell:
#         " {BWA} {params.options} {input.refg_fasta} {input.reads} | samtools sort -o {output.sortedbam} -T {params.tempfile}  > {log.logfile} 2>&1 "
# 
# 
