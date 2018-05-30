#!/usr/bin/env python3.5
import os
# import IPython; 

# set config file
configfile: "./config.json"

BWA        = config["progs"]["BWA"]
MM2        = config["progs"]["minimap"]
SAMTOOLS   = config["progs"]["SAMTOOLS"]
nanopolish = config["progs"]["nanopolish"]

RefTranscriptome = config["ref"]["Transcriptome"]
GENOME_VERSION = config["ref"]["Genome_version"]

DIR_ALIGNED_MINIMAP    = config["PATHOUT"]+"01_MM_aligned/"
DIR_FILTERED_MINIMAP   = config["PATHOUT"]+"02_MM_filtered/"
DIR_SORTED_MINIMAPPED  = config["PATHOUT"]+"03_MM_sortedbam/"
DIR_SORTED_ALIGNED_BWA = config["PATHOUT"]+"04_BWA_sortedbam/"
DIR_EVENTALIGN         = config["PATHOUT"]+"05_BWA_eventalign/"
DIR_REPORT             = config["PATHOUT"]+"06_report/"

Sample_indices_int     = range(config["MAXSAMPLEi"] +1) 
Sample_indices_str     = [ str(item) for item in Sample_indices_int  ]


Ealign_FILES_list = list( chain( *[ expand ( os.path.join( DIR_EVENTALIGN, 'Ealign_'+sample+'.cvs' ), ) for sample in Sample_indices_str ] ) )
 
#                   [
#                    expand ( os.path.join( DIR_REPORT,  config["SAMPLES"][sample]+"_report.html"), ) for sample in config["SAMPLES"],
#                   ]
# print( Ealign_FILES_list )
# Ealign_FILES_ss  = ' '.join(Ealign_FILES_list) 


OUTPUT_FILES=  [
               os.path.join( DIR_EVENTALIGN, 'E_aligned_all.cvs'),
               # expand ( os.path.join( DIR_REPORT,  config["SAMPLES"][sample]+"_report.html"), ) for sample in config["SAMPLES"]
 
               ]
# os.path.join( config['PATHIN'], 'fastq', 'pass', config['RUN_ID']+'_'+sample+'.fastq.index'), )
# print("OUTPUT_FILES=")
# for x in OUTPUT_FILES: 
# IPython.embed()
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

rule consolidate_alignments:
# Pool the alignment files into a single table (with just one header) to gather all the statistics from
    input:
        Ealign_FILES_list
    output: 
        os.path.join( DIR_EVENTALIGN, 'E_aligned_all.cvs') 
    shell:
        " head -1 {Ealign_FILES_list[0]} > '{output}' && tail -q -n +2 {Ealign_FILES_list} >> '{output}' "

#------------------------------------------------------

rule np_event_align:
# Align the events to the reference genome
    input:
        sortedbam             = os.path.join( DIR_SORTED_ALIGNED_BWA, config['RUN_ID']+'_{sample}.bwaligned.sorted.bam'),
        NOTCALLED_indexedbam  = os.path.join( DIR_SORTED_ALIGNED_BWA, config['RUN_ID']+'_{sample}.bwaligned.sorted.bam.bai'),
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
# Index the sorted bam file
    input:
        sortedbam  = os.path.join( DIR_SORTED_ALIGNED_BWA, config['RUN_ID']+'_{sample}.bwaligned.sorted.bam')
    output:
        indexedbam  = os.path.join( DIR_SORTED_ALIGNED_BWA, config['RUN_ID']+'_{sample}.bwaligned.sorted.bam.bai')
    log:
        logfile  = os.path.join( DIR_SORTED_ALIGNED_BWA, 'index_{sample}_bwaMemOnt2d.log')
    message: """---- index the bam files for sample {wildcards.sample} ----"""
    shell:
        " {SAMTOOLS} index  {input.sortedbam}  2> {log.logfile} "

#------------------------------------------------------

rule align_bwa_mem_ont2d:
# Align the reads to the reference
    input:
        refg_fasta = os.path.join(config['ref']['Genome_DIR'] , config['ref']['Genome_version']+ ".fa" ),
        refg_bwt   = os.path.join(config['ref']['Genome_DIR'] , config['ref']['Genome_version']+ ".fa.bwt"),
        reads      = os.path.join( config['PATHIN'], 'fastq', 'pass', config['RUN_ID']+'_{sample}.fastq'),
        npi        = os.path.join( config['PATHIN'], 'fastq', 'pass', config['RUN_ID']+'_{sample}.fastq.index')
    output:
        sortedbam  = os.path.join( DIR_SORTED_ALIGNED_BWA, config['RUN_ID']+'_{sample}.bwaligned.sorted.bam')
    params:
        options    = " mem -x ont2d ",
        tempfile   = os.path.join( DIR_SORTED_ALIGNED_BWA, config['RUN_ID']+'_{sample}.bwaligniment.log')
    log:
        logfile  = os.path.join( DIR_SORTED_ALIGNED_BWA, 'alignment_{sample}_bwaMemOnt2d.log')
    message: """---- Align the reads from sample {wildcards.sample} to the reference ----"""
    shell:
        " {BWA} {params.options} {input.refg_fasta} {input.reads} | samtools sort -o {output.sortedbam} -T {params.tempfile}  > {log.logfile} 2>&1 "

#------------------------------------------------------

rule convert_sort_minimap:
# convert from sam to bam format and sort by position
    input:
        aligned     = os.path.join( DIR_FILTERED_MINIMAP, "{sample}.0filtered.sam")
    output:
        sortedbam   = os.path.join( DIR_SORTED_MINIMAPPED, "{sample}.sorted.bam")
    params:
        options = "-ax splice "
    log:
        logfile = os.path.join( DIR_SORTED_MINIMAPPED, "{sample}.sortbam.log")
    message: """ --- converting, sorting, and indexing bam file. --- """
    shell:
        'samtools view  -Sb  {input} | samtools sort > {output} && samtools index {output} 2> {log.logfile}'

#------------------------------------------------------

rule filter_nonaligned_minimap:
# Check for alignment filter in sam file: if != 4 then remove this read
    input:
        aligned  = os.path.join( DIR_ALIGNED_MINIMAP, "{sample}.sam" )
    output:
        aligned  = os.path.join( DIR_FILTERED_MINIMAP, "{sample}.0filtered.sam" )
    log:
        log      = os.path.join( DIR_FILTERED_MINIMAP, "{sample}_0filtering.log" )
    message: """--- filtering unaligned reads from alignment data ---"""
    shell:
        " cat {input} | perl -lane 'print if $F[1] ne 4'  >  {output}   2> {log}"

#------------------------------------------------------

rule align_minimap:
# use minimap2 to align the fastq reads to the reference genome
    input:
        mmiref   = os.path.join( config['ref']['Genome_DIR'] , config['ref']['Genome_version']+ ".mmi" ),
        sample   = os.path.join( config['PATHIN'], "{sample}.fastq" )
    output:
        aligned  = os.path.join( DIR_ALIGNED_MINIMAP, "{sample}.sam" )
    params:
        options  = " -ax splice "
    log:
        log      = os.path.join( DIR_ALIGNED_MINIMAP, "{sample}_alignment.log")
    message: """--- aligning fastq reads to indexed reference"""
    shell:
        "{MM2} {params} {input.mmiref} {input.sample} > {output}  2> {log}"

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
