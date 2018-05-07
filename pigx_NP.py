#!/usr/bin/env python3.5
import os

# set config file
configfile: "./pigx_NP_config.json"

MM2=config["progs"]["minimap"]
RefTranscriptome = config["ref"]["Transcriptome"]

DIR_ALIGNED  = config["PATHOUT"]+"01_aligned/"
DIR_FILTERED = config["PATHOUT"]+"02_filtered/"
DIR_SORTED   = config["PATHOUT"]+"03_sortedbam/"
DIR_REPORT   = config["PATHOUT"]+"04_report/"

GENOME_VERSION = config["ref"]["Genome_version"]

OUTPUT_FILES =  [ expand ( DIR_REPORT + config["SAMPLES"][sample]+"_report.html",) for sample in config["SAMPLES"] ]

# print("OUTPUT_FILES=")
# for x in OUTPUT_FILES: 
#    print( x )

#=========================================================================================================
#
#                                         BEGIN RULES    
#
#=========================================================================================================
rule all:
    input:
        [ OUTPUT_FILES ]
#------------------------------------------------------
rule make_report:
    input:
        aligned_reads = os.path.join( DIR_SORTED, "{sample}.sorted.bam"),
        transcriptome = RefTranscriptome
    output:
        os.path.join( DIR_REPORT, "{sample}_report.html")
    params:
        " readcov_THRESH = 10;   ",
        " yplotmax = 10000; "
    log:
        logfile = os.path.join( DIR_REPORT, "finale_report_{sample}.log")
    message: """--- producing final report."""

    shell:  ' Rscript -e  \'{params}  fin_Transcript    = "{input.transcriptome}"; fin_readalignment = "{input.aligned_reads}"; Genome_version="{GENOME_VERSION}" ; rmarkdown::render("Nanopore_report.Rmd", output_file = "{output}" ) \'  '

#------------------------------------------------------
rule convert_sort:
    input:
        aligned     = os.path.join( DIR_FILTERED, "{sample}.0filtered.sam")
    output:
        sortedbam   = os.path.join( DIR_SORTED, "{sample}.sorted.bam")
    params:
        options = "-ax splice "
    log:
        logfile = os.path.join( DIR_SORTED, "{sample}.sortbam.log")
    message: """ --- converting, sorting, and indexing bam file. --- """
    shell:
        'samtools view  -Sb  {input} | samtools sort > {output} && samtools index {output} 2> {log.logfile}'
#------------------------------------------------------
rule filter_nonaligned:
    input:
        aligned  = os.path.join( DIR_ALIGNED, "{sample}.sam" ) 
    output:
        aligned  = os.path.join( DIR_FILTERED, "{sample}.0filtered.sam" )
    log:
        log      = os.path.join( DIR_FILTERED, "{sample}_0filtering.log" )
    message: """--- filtering unaligned reads from alignment data ---"""
    shell:
        " cat {input} | perl -lane 'print if $F[1] ne 4'  >  {output}   2> {log}"
 
 
#------------------------------------------------------
rule align:
    input:
        mmiref   = os.path.join( config['ref']['Genome_DIR'] , config['ref']['Genome_version']+ ".mmi" ),
        sample   = os.path.join( config['PATHIN'], "{sample}.fastq" )
    output:
        aligned  = os.path.join( DIR_ALIGNED, "{sample}.sam" )
    params:
        options  = " -ax splice "
    log:
        log      = os.path.join( DIR_ALIGNED, "{sample}_alignment.log")
    message: """--- aligning fastq reads to indexed reference"""
    shell:
        "{MM2} {params} {input.mmiref} {input.sample} > {output}  2> {log}"
 
 
#------------------------------------------------------
# Create version of reference genome for fast alignment later:

rule minimizer:
    input:
        refgenome_fasta  = os.path.join(config['ref']['Genome_DIR'] , config['ref']['Genome_version']+ ".fa" )
    output:
        refgenome_mmiref = os.path.join(config['ref']['Genome_DIR'] , config['ref']['Genome_version']+ ".mmi")
    params:
        options = " -d  "
    log:
        os.path.join( config['ref']['Genome_DIR'], config['ref']['Genome_version'], "_mmi2_minimizer_creation.log")
    message: """--- creating minimizer index of reference genome for minimap2."""
    shell:
        "{MM2} {params.options}  {output} {input} 2> {log}"


