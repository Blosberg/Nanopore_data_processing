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

DIR_SYMLINKS           = config["PATHOUT"]+"01_symlinks/"
os.makedirs( DIR_SYMLINKS,  exist_ok=True)
# -- all the rest will be created by snakemake automatically

DIR_ALIGNED_MINIMAP    = config["PATHOUT"]+"02_MM_aligned/"
DIR_FILTERED_MINIMAP   = config["PATHOUT"]+"03_MM_filtered/"
DIR_SORTED_MINIMAPPED  = config["PATHOUT"]+"04_MM_sortedbam/"
DIR_SORTED_ALIGNED_BWA = config["PATHOUT"]+"05_BWA_sortedbam/"
DIR_EVENTALIGN         = config["PATHOUT"]+"06_BWA_eventalign/"
DIR_GR                 = config["PATHOUT"]+"07_GRobjects"
DIR_REPORT             = config["PATHOUT"]+"08_report/"
DIR_REFGEMONE          = config['ref']['Genome_DIR']


#------------------------------------------------------
# Create symbolic links to PATHIN so that indexing/etc can be performed in written pathout 

if ( intype == "raw_minION"): 
   # the snakemake rules defined in the following included script assume that the minION
   # output has been "chunked" into sequential files, and will have to be reassembled. 
   include   : os.path.join( config["scripts"]["script_folder"], config["scripts"]["rules_chunks"] )

   for sample in config["samplelist"]:
       for index in range(0, config['samplelist'][sample]["MAXSAMPLEi"] + 1):
          linkname = config['samplelist'][sample]["RUN_ID"] + "_" + str(index) + config['samplelist'][sample]["fastq_suffix"]

          source   = getPathCase( config["PATHIN"], "fastq", "pass",config["samplelist"][sample]["fastq_prefix"] + str(index) +  config["samplelist"][sample]["fastq_suffix"], "raw_minION") 

          makelink(  source, os.path.join( DIR_SYMLINKS, linkname) )

 
elif( intype == "fastq"):
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

rule np_event_align:
# Align the events to the reference genome the wildcard "chunk" can simply be "full", in cases where there are no chunks
    input:
        sortedbam             = os.path.join( DIR_SORTED_ALIGNED_BWA, "chunks", "fastq_run_{sample}_{chunk}.bwaligned.sorted.bam"),
        NOTCALLED_indexedbam  = os.path.join( DIR_SORTED_ALIGNED_BWA, "chunks", "fastq_run_{sample}_{chunk}.bwaligned.sorted.bam.bai"),
        fastq_file            = os.path.join( DIR_SYMLINKS,  "{sample}_{chunk}" +config["samplelist"][sample]["fastq_suffix"]),
        NOTCALLED_fastq_npi   = os.path.join( DIR_SYMLINKS,  "{sample}_{chunk}" + config["samplelist"][sample]["fastq_suffix"]),
        refgenome_fasta       = os.path.join( DIR_REFGEMONE, config['ref']['Genome_version']+ ".fa" ),
        NOTCALLED_bwt         = os.path.join( DIR_REFGEMONE, config['ref']['Genome_version']+ ".fa.bwt"),
        NOTCALLED_pac         = os.path.join( DIR_REFGEMONE, config['ref']['Genome_version']+ ".fa.pac")
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

rule index_sortedbam:
# Index the sorted bam file
    input:
        sortedbam  = os.path.join( DIR_SORTED_ALIGNED_BWA, "chunks", "fastq_run_{sample}_{chunk}.bwaligned.sorted.bam")
    output:
        indexedbam = os.path.join( DIR_SORTED_ALIGNED_BWA, "chunks", "fastq_run_{sample}_{chunk}.bwaligned.sorted.bam.bai")
    log:
        logfile  = os.path.join( DIR_SORTED_ALIGNED_BWA, "chunks", 'index_{sample}_{chunk}_bwaMemOnt2d.log')
    message: """---- index the bam files for {wildcards.sample} chunk {wildcards.chunk} ----"""
    shell:
        " {SAMTOOLS} index  {input.sortedbam}  2> {log.logfile} "

#------------------------------------------------------

rule align_bwa_mem_ont2d:
# Align the reads to the reference
    input:
        refg_fasta = os.path.join(DIR_REFGEMONE , config['ref']['Genome_version']+ ".fa" ),
        refg_bwt   = os.path.join(DIR_REFGEMONE , config['ref']['Genome_version']+ ".fa.bwt"),
        reads      = os.path.join( DIR_SYMLINKS,  "{sample}_{chunk}" + config["samplelist"][sample]["fastq_suffix"] ),
        npi        = os.path.join( DIR_SYMLINKS,  "{sample}_{chunk}" + config["samplelist"][sample]["fastq_suffix"] + ".index")
    output:
        sortedbam  = os.path.join( DIR_SORTED_ALIGNED_BWA, "chunks", "fastq_run_{sample}_{chunk}.bwaligned.sorted.bam")
    params:
        options    = " mem -x ont2d ",
        tempfile   = os.path.join( DIR_SORTED_ALIGNED_BWA, "chunks", "fastq_run_{sample}_{chunk}.bwaligniment.log")
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
        fastq_file   = os.path.join( DIR_SYMLINKS, config['samplelist'][sample]["RUN_ID"] + "_" + str(index) + config['samplelist'][sample]["fastq_suffix"] ) 
    output:
        npi    = os.path.join( DIR_SYMLINKS, "{sample}_{chunk}"+config["samplelist"][sample]["fastq_suffix"]+".index"  ), 
        fai    = os.path.join( DIR_SYMLINKS, "{sample}_{chunk}"+config["samplelist"][sample]["fastq_suffix"]+".index.fai" ),
        gzi    = os.path.join( DIR_SYMLINKS, "{sample}_{chunk}"+config["samplelist"][sample]["fastq_suffix"]+".index.gzi" ),  
        readdb = os.path.join( DIR_SYMLINKS, "{sample}_{chunk}"+config["samplelist"][sample]["fastq_suffix"]+".index.readdb")
    params:
        options    = " index -d "
    log:
        logfile  = os.path.join( DIR_SYMLINKS,  "{sample}_{chunk}_npi.log", intype  )
    message: """---- index the reads from chunk {wildcards.chunk} against the fast5 files from the same. ----"""
    shell:
        " nice -19 {nanopolish} {params.options} {input.fast5_folder} {input.fastq_file} 2> {log.logfile} "
 

rule bwa_index:
# Create indexed version of reference genome for fast alignment with bwa later:
    input:
        refgenome_fasta  = os.path.join(DIR_REFGEMONE , config['ref']['Genome_version']+ ".fa" )
    output:
        bwt      = os.path.join(DIR_REFGEMONE , config['ref']['Genome_version']+ ".fa.bwt"),
        pac      = os.path.join(DIR_REFGEMONE , config['ref']['Genome_version']+ ".fa.pac")
    params:
        options  = " index  "
    log:
        logfile  = os.path.join( DIR_REFGEMONE, config['ref']['Genome_version'], "_bwa_indexing.log")
    message: 
        """---- creating bwa index of the reference genome. ----"""
    shell:
        "{BWA} {params.options}  {input} > {log.logfile}"
