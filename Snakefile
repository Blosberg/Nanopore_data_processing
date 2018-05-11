#!/usr/bin/env python3.5
import os

# set config file
configfile: "./config.json"

BWA=config["progs"]["BWA"]
MM2=config["progs"]["minimap"]
SAMTOOLS=config["progs"]["SAMTOOLS"]
nanopolish=config["progs"]["nanopolish"]

RefTranscriptome = config["ref"]["Transcriptome"]
GENOME_VERSION = config["ref"]["Genome_version"]


# DIR_ALIGNED  = config["PATHOUT"]+"01_aligned/"
# DIR_FILTERED = config["PATHOUT"]+"02_filtered/"
DIR_SORTED     = config["PATHOUT"]+"03_sortedbam/"
DIR_EVENTALIGN = config["PATHOUT"]+"04_eventalign/"
# DIR_REPORT   = config["PATHOUT"]+"04_report/"

OUTPUT_FILES = [ 
                expand ( os.path.join( DIR_EVENTALIGN, 'Ealign_'+sample+'.cvs' ), ) for sample in config["SAMPLES"] 
               ] 
		
# os.path.join( config['PATHIN'], 'fastq', 'pass', config['RUN_ID']+'_'+sample+'.fastq.index'), )
# print("OUTPUT_FILES=")
# for x in OUTPUT_FILES: 
#    print( x )
# 
#=========================================================================================================
#
#                                         BEGIN RULES    
#
#=========================================================================================================
rule all:
    input:
        [ OUTPUT_FILES ]
#------------------------------------------------------
rule np_event_align:
# align the events to the reference genome
    input:
        sortedbam             = os.path.join( DIR_SORTED, config['RUN_ID']+'_{sample}.bwaligned.sorted.bam'),
        NOTCALLED_indexedbam  = os.path.join( DIR_SORTED, config['RUN_ID']+'_{sample}.bwaligned.sorted.bam.bai'),
	fastq_file            = os.path.join( config['PATHIN'], 'fastq', 'pass', config['RUN_ID']+'_{sample}.fastq'), 
        NOTCALLED_fastq_npi   = os.path.join( config['PATHIN'], 'fastq', 'pass', config['RUN_ID']+'_{sample}.fastq.index'), 
        refgenome_fasta  = os.path.join(config['ref']['Genome_DIR'] , config['ref']['Genome_version']+ ".fa" ),
        NOTCALLED_bwt    = os.path.join(config['ref']['Genome_DIR'] , config['ref']['Genome_version']+ ".fa.bwt"),
        NOTCALLED_pac    = os.path.join(config['ref']['Genome_DIR'] , config['ref']['Genome_version']+ ".fa.pac")
    output:
        Ealigned         = os.path.join( DIR_EVENTALIGN, 'Ealign_{sample}.cvs' )
    log:
        logfile  = os.path.join( DIR_EVENTALIGN, 'Ealign_{sample}.log')
    message: """---- align events from sample {wildcards.sample} to the genome ----"""
    shell:
        " {nanopolish} eventalign --reads {input.fastq_file} --bam {input.sortedbam} --genome {input.refgenome_fasta} --scale-events  > {output}  2> {log.logfile} "




#------------------------------------------------------
# rule quickcheck:
# TODO:


#------------------------------------------------------
rule index_sortedbam:
# index the sorted bam file
    input:
        sortedbam  = os.path.join( DIR_SORTED, config['RUN_ID']+'_{sample}.bwaligned.sorted.bam')
    output:
        indexedbam  = os.path.join( DIR_SORTED, config['RUN_ID']+'_{sample}.bwaligned.sorted.bam.bai')
    log:
        logfile  = os.path.join( DIR_SORTED, 'index_{sample}_bwaMemOnt2d.log')
    message: """---- index the bam files for sample {wildcards.sample} ----"""
    shell:
        " {SAMTOOLS} index  {input.sortedbam}  2> {log.logfile} "

#------------------------------------------------------

rule bwa_mem_ont2d:
# Align the reads to the reference
    input:
        refg_fasta = os.path.join(config['ref']['Genome_DIR'] , config['ref']['Genome_version']+ ".fa" ),
        refg_bwt   = os.path.join(config['ref']['Genome_DIR'] , config['ref']['Genome_version']+ ".fa.bwt"),
        reads      = os.path.join( config['PATHIN'], 'fastq', 'pass', config['RUN_ID']+'_{sample}.fastq'),
        npi        = os.path.join( config['PATHIN'], 'fastq', 'pass', config['RUN_ID']+'_{sample}.fastq.index')
    output:
        sortedbam  = os.path.join( DIR_SORTED, config['RUN_ID']+'_{sample}.bwaligned.sorted.bam')
    params:
        options    = " mem -x ont2d ",
        tempfile   = os.path.join( DIR_SORTED, config['RUN_ID']+'_{sample}.bwaligniment.log')
    log:
        logfile  = os.path.join( DIR_SORTED, 'alignment_{sample}_bwaMemOnt2d.log')
    message: """---- Align the reads from sample {wildcards.sample} to the reference ----"""
    shell:
        " {BWA} {params.options} {input.refg_fasta} {input.reads} | samtools sort -o {output.sortedbam} -T {params.tempfile}  2> {log.logfile} "

#------------------------------------------------------

rule np_index:
# Index the reads and the fast5 files themselves
    input:
        fast5_folder = os.path.join( config['PATHIN'], 'fast5', 'pass', '{sample}' ),
	fastq_file   = os.path.join( config['PATHIN'], 'fastq', 'pass', config['RUN_ID']+'_{sample}.fastq') 
    output:
        npi    = os.path.join( config['PATHIN'], 'fastq', 'pass', config['RUN_ID']+'_{sample}.fastq.index'), 
	fai    = os.path.join( config['PATHIN'], 'fastq', 'pass', config['RUN_ID']+'_{sample}.fastq.index.fai'),
	gzi    = os.path.join( config['PATHIN'], 'fastq', 'pass', config['RUN_ID']+'_{sample}.fastq.index.gzi'),
	readdb = os.path.join( config['PATHIN'], 'fastq', 'pass', config['RUN_ID']+'_{sample}.fastq.index.readdb')
    params:
        options    = " index -d "
    log:
        logfile  = os.path.join( config['PATHIN'], 'fastq', 'pass', config['RUN_ID']+'_{sample}_npi.log' )
    message: """---- index the reads from sample {wildcards.sample} against the fast5 files from the same. ----"""
    shell:
        " nice -19 {nanopolish} {params.options} {input.fast5_folder} {input.fastq_file} 2> {log.logfile} "
        
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
