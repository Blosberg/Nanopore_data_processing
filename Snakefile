#!/usr/bin/env python3.5
import os
# import IPython; 

# set config file
configfile: "./config.json"
include   : os.path.join( config["scripts"]["script_folder"], config["scripts"]["pyfunc_defs"] )

#------------------------------------------------------
# --- Dependencies:

BWA        = config["progs"]["BWA"]
MM2        = config["progs"]["minimap"]
SAMTOOLS   = config["progs"]["SAMTOOLS"]
nanopolish = config["progs"]["nanopolish"]

RefTranscriptome = config["ref"]["Transcriptome"]
GENOME_VERSION   = config["ref"]["Genome_version"]
RmdReportScript  = "scripts/Nanopore_report.Rmd"
intype = config["intype"]

#------------------------------------------------------
# --- define output directories

DIR_ALIGNED_MINIMAP    = config["PATHOUT"]+"01_MM_aligned/"
DIR_FILTERED_MINIMAP   = config["PATHOUT"]+"02_MM_filtered/"
DIR_SORTED_MINIMAPPED  = config["PATHOUT"]+"03_MM_sortedbam/"
DIR_SORTED_ALIGNED_BWA = config["PATHOUT"]+"04_BWA_sortedbam/"
DIR_EVENTALIGN         = config["PATHOUT"]+"05_BWA_eventalign/"
DIR_GR                 = config["PATHOUT"]+"06_GRobjects"
DIR_REPORT             = config["PATHOUT"]+"07_report/"


#------------------------------------------------------
# define file lists based on input formatting.

if ( intype == "RAW_minION"): 
   # the snakemake rules defined in the following included script assume that the minION
   # output has been "chunked" into sequential files, and will have to be reassembled. 
   include   : os.path.join( config["scripts"]["script_folder"], config["scripts"]["rules_chunks"] )
  
elif( intype == "fastq"):
   # Do some other stuff
   include   : os.path.join( config["scripts"]["script_folder"], config["scripts"]["rules_wholefastq"] )

else: 
   print("Unrecognized input data format. Terminating.")
   exit(1)

#------------------------------------------------------
# Define output (target) files:

if ( config["target_out"] == "report" ):
   OUTPUT_FILES=  [
                  # os.path.join( DIR_EVENTALIGN, 'E_aligned_all.cvs'),
                  os.path.join( DIR_REPORT, ""+config["samplelist"][sample]["RUN_ID"]+"_report.html") for sample in config["samplelist"]
                  ]
elif ( config["target_out"] == "GR" ):
   OUTPUT_FILES=  [
                  # os.path.join( DIR_EVENTALIGN, 'E_aligned_all.cvs'),
                  os.path.join( DIR_GR, config["samplelist"][sample]["RUN_ID"]+"_GR.RData")  for sample in config["samplelist"]
                  ]
elif ( config["target_out"] == "bam" ):
   OUTPUT_FILES=  [
                  os.path.join( DIR_SORTED_MINIMAPPED, config["samplelist"][sample]["RUN_ID"]+".sorted.bam") for sample in config["samplelist"]
                  ]
else:
   print("Unrecognized target output file format: ", config["target_out"], " ... Terminating.")
   exit(1)

#------------------------------------------------------
print("intype = " + config["intype"])
print("target out = " + config["target_out"])
print("OUTPUT_FILES=")
for x in OUTPUT_FILES: 
  print(x)
print("\n finished outputting output files \n\n ")
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

rule make_report:
# build the final output report in html format
    input:
        aligned_reads = os.path.join( DIR_SORTED_MINIMAPPED, "{sample}.sorted.bam"),
        transcriptome = RefTranscriptome,
        GRobj         = os.path.join( DIR_GR, "{sample}_GR.RData")
    output:
        os.path.join( DIR_REPORT, "{sample}_report.html")
    params:
        " readcov_THRESH = 10;   ",
        " yplotmax = 10000; "
    log:
        logfile = os.path.join( DIR_REPORT, "finale_report_{sample}.log")
    message: """--- producing final report."""

    shell: 
        " Rscript -e  '{params} "
        " fin_readalignment = \"{input.aligned_reads}\"; "
        " fin_Transcript = \"{input.transcriptome}\";"
        " fin_GRobj      = \"{input.GRobj}\";"
        " Genome_version=\"{GENOME_VERSION}\"; " 
        " rmarkdown::render(\"{RmdReportScript}\", output_file = \"{output}\" ) ' "  

#------------------------------------------------------

rule np_event_align:
# Align the events to the reference genome the wildcard "chunk" can simply be "full", in cases where there are no chunks
    input:
        sortedbam             = os.path.join( DIR_SORTED_ALIGNED_BWA, "chunks", "fastq_runid_{sample}_{chunk}.bwaligned.sorted.bam"),
        NOTCALLED_indexedbam  = os.path.join( DIR_SORTED_ALIGNED_BWA, "chunks", "fastq_runid_{sample}_{chunk}.bwaligned.sorted.bam.bai"),
        fastq_file            = os.path.join( config['PATHIN'], 'fastq', 'pass',"fastq_runid_{sample}_{chunk}.fastq"),
        NOTCALLED_fastq_npi   = os.path.join( config['PATHIN'], 'fastq', 'pass',"fastq_runid_{sample}_{chunk}.fastq.index"),
        refgenome_fasta  = os.path.join(config['ref']['Genome_DIR'] , config['ref']['Genome_version']+ ".fa" ),
        NOTCALLED_bwt    = os.path.join(config['ref']['Genome_DIR'] , config['ref']['Genome_version']+ ".fa.bwt"),
        NOTCALLED_pac    = os.path.join(config['ref']['Genome_DIR'] , config['ref']['Genome_version']+ ".fa.pac")
    output:
        Ealigned         = os.path.join( DIR_EVENTALIGN, 'Ealign_{sample}_{chunk}.cvs' )
    log:
        logfile  = os.path.join( DIR_EVENTALIGN, 'Ealign_{sample}_{chunk}.log')
    message: """---- align events from sample {wildcards.sample}, chunk {wildcards.chunk} to the genome ----"""
    shell:
        " {nanopolish} eventalign --reads {input.fastq_file} --bam {input.sortedbam} --genome {input.refgenome_fasta} --scale-events  > {output}  2> {log.logfile} "


#------------------------------------------------------
# rule quickcheck: (TODO)
#------------------------------------------------------

rule index_sortedbam:
# Index the sorted bam file
    input:
        sortedbam  = os.path.join( DIR_SORTED_ALIGNED_BWA, "chunks", "fastq_runid_{sample}_{chunk}.bwaligned.sorted.bam")
    output:
        indexedbam = os.path.join( DIR_SORTED_ALIGNED_BWA, "chunks", "fastq_runid_{sample}_{chunk}.bwaligned.sorted.bam.bai")
    log:
        logfile  = os.path.join( DIR_SORTED_ALIGNED_BWA, "chunks", 'index_{sample}_{chunk}_bwaMemOnt2d.log')
    message: """---- index the bam files for {wildcards.sample} chunk {wildcards.chunk} ----"""
    shell:
        " {SAMTOOLS} index  {input.sortedbam}  2> {log.logfile} "

#------------------------------------------------------

rule align_bwa_mem_ont2d:
# Align the reads to the reference
    input:
        refg_fasta = os.path.join(config['ref']['Genome_DIR'] , config['ref']['Genome_version']+ ".fa" ),
        refg_bwt   = os.path.join(config['ref']['Genome_DIR'] , config['ref']['Genome_version']+ ".fa.bwt"),
        reads      = os.path.join( config['PATHIN'], 'fastq', 'pass', "fastq_runid_{sample}_{chunk}.fastq"),
        npi        = os.path.join( config['PATHIN'], 'fastq', 'pass', "fastq_runid_{sample}_{chunk}.fastq.index")
    output:
        sortedbam  = os.path.join( DIR_SORTED_ALIGNED_BWA, "chunks", "fastq_runid_{sample}_{chunk}.bwaligned.sorted.bam")
    params:
        options    = " mem -x ont2d ",
        tempfile   = os.path.join( DIR_SORTED_ALIGNED_BWA, "chunks", "fastq_runid_{sample}_{chunk}.bwaligniment.log")
    log:
        logfile  = os.path.join( DIR_SORTED_ALIGNED_BWA, "chunks", 'alignment_{sample}_{chunk}_bwaMemOnt2d.log')
    message: """---- Align the reads from chunk {wildcards.chunk} to the reference ----"""
    shell:
        " {BWA} {params.options} {input.refg_fasta} {input.reads} | samtools sort -o {output.sortedbam} -T {params.tempfile}  > {log.logfile} 2>&1 "

#------------------------------------------------------

rule np_index:
# Index the reads and the fast5 files themselves
    input:
        fast5_folder = getPathCase( config['PATHIN'], 'fast5', 'pass', '{chunk}', intype ),
        fastq_file   = getPathCase( config['PATHIN'], 'fastq', 'pass', '{sample}_pass.fq.gz', intype ) 
    output:
        npi    = getPathCase( config['PATHIN'], 'fastq', 'pass', "fastq_runid_{sample}_{chunk}.fastq.index", intype ), 
        fai    = getPathCase( config['PATHIN'], 'fastq', 'pass', "fastq_runid_{sample}_{chunk}.fastq.index.fai", intype ),
        gzi    = getPathCase( config['PATHIN'], 'fastq', 'pass', "fastq_runid_{sample}_{chunk}.fastq.index.gzi", intype ),  
        readdb = getPathCase( config['PATHIN'], 'fastq', 'pass', "fastq_runid_{sample}_{chunk}.fastq.index.readdb", intype )
    params:
        options    = " index -d "
    log:
        logfile  = getPathCase( config['PATHIN'], 'fastq', 'pass', "fastq_runid_"+'{sample}'+'_{chunk}_npi.log', intype  )
    message: """---- index the reads from chunk {wildcards.chunk} against the fast5 files from the same. ----"""
    shell:
        " nice -19 {nanopolish} {params.options} {input.fast5_folder} {input.fastq_file} 2> {log.logfile} "
# OUTPUT_FILES  = getPathCase( config['PATHIN'], 'fastq', 'pass', "fastq_runid_"+singleID'_full.fastq.index.fai'),
#------------------------------------------------------

rule bwa_index:
# Create indexed version of reference genome for fast alignment with bwa later:
    input:
        refgenome_fasta  = os.path.join(config['ref']['Genome_DIR'] , config['ref']['Genome_version']+ ".fa" )
    output:
        bwt  = os.path.join(config['ref']['Genome_DIR'] , config['ref']['Genome_version']+ ".fa.bwt"),
        pac  = os.path.join(config['ref']['Genome_DIR'] , config['ref']['Genome_version']+ ".fa.pac")
    params:
        options  = " index  "
    log:
        logfile  = os.path.join( config['ref']['Genome_DIR'], config['ref']['Genome_version'], "_bwa_indexing.log")
    message: """---- creating bwa index of the reference genome. ----"""
    shell:
        "{BWA} {params.options}  {input} > {log.logfile}"



