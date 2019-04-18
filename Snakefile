#!/home/bosberg/projects/nanopore/dev/guix/.guix-profile/bin/python3
import os, sys, json, csv, yaml
import argparse

# ------------------------------------------------------
# --- Define Dependencies:

MM2        = config["progs"]["minimap"]
SAMTOOLS   = config["progs"]["SAMTOOLS"]
nanopolish = config["progs"]["nanopolish"]

RefTranscriptome   = config["ref"]["Transcriptome"]
GENOME_VERSION     = config["ref"]["Genome_version"]
RmdReportScript    = os.path.join(config["scripts"]["script_folder"],"final_report","Nanopore_report.Rmd")
input_data_type    = config["input_data_type"]

R_tables2GR_main     = os.path.join( config["scripts"]["script_folder"], config["scripts"]["Rmain_tsv2GRconv"] )
R_tables2GR_funcs    = os.path.join( config["scripts"]["script_folder"], config[ "scripts"]["Rfuncs_tsv2GRconv"] )

R_flattenreads_main  = os.path.join( config["scripts"]["script_folder"], config["scripts"]["Rmain_flattenreads"] )
R_flattenreads_funcs = os.path.join( config["scripts"]["script_folder"], config["scripts"]["Rfuncs_flattenreads"] )

R_build_histlist_main   = os.path.join( config["scripts"]["script_folder"], config["scripts"]["Rmain_build_histlist"] )

#------------------------------------------------------
# --- Define output directories

SUBDIR_SYMLINKS           = "01_symlinks_fastqindex_chunks/"
SUBDIR_ALIGNED_MINIMAP    = "02_MM_aligned_chunks/"
SUBDIR_FILTERED_MINIMAP   = "03_MM_filtered_chunks/"
SUBDIR_SORTED_MINIMAPPED  = "04_MM_sortedbam/"
SUBDIR_EVENTALIGN         = "05_eventalign/"
SUBDIR_GR                 = "06_GRobjects"
SUBDIR_REPORT             = "Final_report/"

DIR_REFGENOME             = config['ref']['Genome_DIR']

#------------------------------------------------------
# --- Include function definitions and rules

include: os.path.join( config["scripts"]["script_folder"], config["scripts"]["pyfunc_defs"] )
include: os.path.join( config["scripts"]["script_folder"], config["scripts"]["rules_chunks"] )

#------------------------------------------------------
# --- Define output (target) files:

# Initialize empty list
OUTPUT_FILES = []

# increment by sample:
for sampleLoopi_targets in config["samplelist"]:

   # @@@ TODO: implement sample-dependent targets with defaults.
   if ( config["target_out"] == "report" ):
      OUTPUT_FILES.extend(
                          [ os.path.join( config["PATHOUT"], config["samplelist"][sampleLoopi_targets]["sampledir"], SUBDIR_REPORT, "" + sampleLoopi_targets + "_report.html") ]
                          )
   elif ( config["target_out"] == "histlist" ):
      OUTPUT_FILES.extend(
                          [ os.path.join( config["PATHOUT"], config["samplelist"][sampleLoopi_targets]["sampledir"], SUBDIR_GR, sampleLoopi_targets + "_kmer_histlist.rds") ]
                          )
   elif ( config["target_out"] == "flatreads_GRL" ):
      OUTPUT_FILES.extend(
                          [ os.path.join( config["PATHOUT"], config["samplelist"][sampleLoopi_targets]["sampledir"], SUBDIR_GR,  sampleLoopi_targets + "_reads_flat_GRL.rds") ]
                         )
   elif ( config["target_out"] == "reads_GRL" ):
      OUTPUT_FILES.extend(
                          [ os.path.join( config["PATHOUT"], config["samplelist"][sampleLoopi_targets]["sampledir"], SUBDIR_GR, sampleLoopi_targets + "_reads_GRL.rds") ]
                          )
   elif ( config["target_out"] == "mergedbam" ):
      OUTPUT_FILES.extend(
                          [ os.path.join( config["PATHOUT"], config["samplelist"][sampleLoopi_targets]["sampledir"], SUBDIR_SORTED_MINIMAPPED, sampleLoopi_targets + ".sorted.bam") ]
                         )
   elif ( config["target_out"] == "aligned_chunks"):
      OUTPUT_FILES.extend(
                          get_chunkfiles( sampleLoopi_targets, os.path.join( config["PATHOUT"], config["samplelist"][sampleLoopi_targets]["sampledir"], SUBDIR_ALIGNED_MINIMAP) ,  "", ".sam", 0 )
                         )
   else:
      print("Unrecognized target output file format: ", config["target_out"], " ... Terminating.")
      exit(1)

#---  DEBUGGING:
#------------------------------------------------------
# print("input_data_type = " + config["input_data_type"])
# print("target out = " + config["target_out"])
# IPython.embed()
# print("---- last check before rules: ------ ")
# print ( "len(OUTPUT_FILES)=")
# print ( len( OUTPUT_FILES) )
# print(" OUTPUT_FILES=")
# for x in OUTPUT_FILES:
#   print(x)
# print("\n finished outputting output files \n\n ")
#
#
# ========================================================================
#
#   BEGIN RULES
#
# ========================================================================

rule all:
    input:
        OUTPUT_FILES
    message:
        fmt("Target output files:\n" + "\n".join(OUTPUT_FILES) + "\n" )

#------------------------------------------------------
rule make_report:
    # build the final output report in html format
    input:
        aligned_reads_bam = os.path.join( config["PATHOUT"], "{wcreport_sampleDir}", SUBDIR_SORTED_MINIMAPPED, "{wcreport_samplename}.sorted.bam"),
        transcriptome     = RefTranscriptome,
        GRLreads          = os.path.join( config["PATHOUT"], "{wcreport_sampleDir}", SUBDIR_GR, "{wcreport_samplename}_reads_GRL.rds")
    output:
        os.path.join( config["PATHOUT"], "{wcreport_sampleDir}", SUBDIR_REPORT, "{wcreport_samplename}_report.html")
    params:
        " readcov_THRESH = 10;   ",
        " yplotmax = 10000; "
    log:
        logfile = os.path.join( config["PATHOUT"], "{wcreport_sampleDir}", SUBDIR_REPORT, "final_report_{wcreport_samplename}.log")
    message:
        fmt("Compiling final report")
    shell:
        " Rscript -e  '{params} "
        " fin_readalignment_bam = \"{input.aligned_reads_bam}\"; "
        " fin_Transcript        = \"{input.transcriptome}\";"
        " fin_GRLreads          = \"{input.GRLreads}\";"
        " Genome_version        =\"{GENOME_VERSION}\"; "
        " rmarkdown::render(\"{RmdReportScript}\", output_file = \"{output}\" ) ' "

#------------------------------------------------------
rule np_event_align:
    # Align the events to the reference genome.
    # The wildcard "chunk" can simply be "full", in cases
    # where there are no chunks
    input:
        sortedbam             = os.path.join( config["PATHOUT"], "{wcEvalign_sampleDir}", SUBDIR_SORTED_MINIMAPPED, "bam_chunks", "{wcEvalign_samplename}_{wcEvalign_chunk}.sorted.bam"),
        NOTCALLED_indexedbam  = os.path.join( config["PATHOUT"], "{wcEvalign_sampleDir}", SUBDIR_SORTED_MINIMAPPED, "bam_chunks", "{wcEvalign_samplename}_{wcEvalign_chunk}.sorted.bam.bai"),
        fastq_file            = os.path.join( config["PATHOUT"], "{wcEvalign_sampleDir}", SUBDIR_SYMLINKS,  "{wcEvalign_samplename}_{wcEvalign_chunk}" +config["fastq_suffix"]),
        fastq_npi             = os.path.join( config["PATHOUT"], "{wcEvalign_sampleDir}", SUBDIR_SYMLINKS,  "{wcEvalign_samplename}_{wcEvalign_chunk}" + config["fastq_suffix"] + ".index"),
        refgenome_fasta       = os.path.join( DIR_REFGENOME, config['ref']['Genome_version']+ ".fa" ),
        NOTCALLED_bwt         = os.path.join( DIR_REFGENOME, config['ref']['Genome_version']+ ".fa.bwt"),
        NOTCALLED_pac         = os.path.join( DIR_REFGENOME, config['ref']['Genome_version']+ ".fa.pac")
    output:
        Evaligned         = os.path.join( config["PATHOUT"], "{wcEvalign_sampleDir}", SUBDIR_EVENTALIGN, "tsv_chunks", 'Ealign_{wcEvalign_samplename}_{wcEvalign_chunk}.tsv' )
    log:
        logfile  = os.path.join( config["PATHOUT"], "{wcEvalign_sampleDir}", SUBDIR_EVENTALIGN, "tsv_chunks", 'Ealign_{wcEvalign_samplename}_{wcEvalign_chunk}.log')
    message: """---- Align events from sample {wildcards.wcEvalign_samplename}, chunk {wildcards.wcEvalign_chunk} to the genome ----"""
    shell:
        " {nanopolish} eventalign --reads {input.fastq_file} --bam {input.sortedbam} --genome {input.refgenome_fasta} --scale-events  > {output}  2> {log.logfile} "

#------------------------------------------------------
# rule quickcheck: (TODO)
#------------------------------------------------------

rule index_sortedbam:
    # Index the sorted bam file with samtools
    input:
        sortedbam  = os.path.join( config["PATHOUT"], "{wcindexbam_sampleDir}", SUBDIR_SORTED_MINIMAPPED, "bam_chunks", "{wcindexbam_samplename}_{wcindexbam_chunk}.sorted.bam")
    output:
        indexedbam = os.path.join( config["PATHOUT"], "{wcindexbam_sampleDir}", SUBDIR_SORTED_MINIMAPPED, "bam_chunks", "{wcindexbam_samplename}_{wcindexbam_chunk}.sorted.bam.bai")
    log:
        logfile    = os.path.join( config["PATHOUT"], "{wcindexbam_sampleDir}", SUBDIR_SORTED_MINIMAPPED, "bam_chunks", '{wcindexbam_samplename}_{wcindexbam_chunk}_samtoolsindex.log')
    message: """---- index the bam files for {wildcards.wcindexbam_samplename} chunk {wildcards.wcindexbam_chunk} ----"""
    shell:
        " {SAMTOOLS} index  {input.sortedbam}  2> {log.logfile} "

#------------------------------------------------------

