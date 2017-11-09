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


OUTPUT_FILES = [ expand ( config["PATHOUT"]+sample+"_report.html",) for sample in config["SAMPLES"] ] 

# print("OUTPUT_FILES=")
# for x in OUTPUT_FILES: 
#    print( x)
# 
# exit()
# ==============================================================================================================
#
#                                         BEGIN RULES    
#
# ==============================================================================================================

rule all:
    input:
        [ OUTPUT_FILES ]

#------------------------------------------------------
rule make_report:
    input:
        aligned_reads = DIR_SORTED+"{sample}.sorted.bam",
        transcriptome = RefTranscriptome
    output:
        DIR_REPORT+"{sample}_report.html"
    params:
        " readcov_THRESH = 100;   ",
        " yplotmax = 10000; "
    log:
        DIR_REPORT+"finale_report_{sample}.log"
    message: """--- producing final report."""

    shell:  ' Rscript -e  \'{params}  fin_Transcript    = "{input.transcriptome}"; fin_readalignment = "{input.aligned_reads}"; Genome_version=config["ref"]["Genome_version"] ; rmarkdown::render("Nanopore_report.Rmd", output_file = "{output}" ) \'  '

#------------------------------------------------------
rule convert_sort:
    input:
        aligned     = DIR_FILTERED+"{sample}.0filtered.sam"
    output:
        sortedbam   = DIR_SORTED+"{sample}.sorted.bam"
    params:
        options = "-ax splice "
    log:
        DIR_SORTED+"{sample}.sortbam.log"
    message: """ --- converting, sorting, and indexing bam file. --- """
    shell:
        'samtools view  -Sb  {input} | samtools sort > {output} && samtools index {output} '
 
#------------------------------------------------------
rule filter_nonaligned:
    input:
        aligned  = DIR_ALIGNED+"{sample}.sam"
    output:
        aligned  = DIR_FILTERED+"{sample}.0filtered.sam"
    log:
        log      = DIR_FILTERED+"{sample}_0filtering.log"
    message: """--- filtering unaligned reads from alignment data ---"""
    shell:
        " cat {input} | perl -lane 'print if $F[1] ne 4'  >  {output}   2> {log}"
 
 
#------------------------------------------------------
rule align:
    input:
        mmiref   = config['ref']['Genome_DIR']+config['ref']['Genome_version']+".mmi",
        sample   = config['PATHIN']+"{sample}.fastq.gz"
    output:
        aligned  = DIR_ALIGNED+"{sample}.sam"
    params:
        options  = " -ax splice "
    log:
        log      = DIR_ALIGNED+"{sample}_alignment.log"
    message: """--- aligning fastq reads to indexed reference"""
    shell:
        "{MM2} {params} {input.mmiref} {input.sample} > {output}  2> {log}"
 
 
# #------------------------------------------------------
rule minimizer:
    input:
        config['ref']['Genome_DIR']+config['ref']['Genome_version']+".fa"
    output:
        mmiref = config['ref']['Genome_DIR']+config['ref']['Genome_version']+".mmi" 
    params:
        options = " -d  "
    log:
        config['ref']['Genome_DIR']+config['ref']['Genome_version']+"_mmi2_minimizer_creation.log"
    message: """--- creating minimizer index of reference genome for minimap2."""
    shell:
        "{MM2} {params.options}  {output} {input} 2> {log}"


