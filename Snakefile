#!/home/bosberg/projects/nanopore/dev/guix/.guix-profile/bin/python3
import os, sys, json, csv, yaml
import argparse

# ------------------------------------------------------
# --- Define Dependencies:

MM2        = config["progs"]["minimap"]
SAMTOOLS   = config["progs"]["SAMTOOLS"]
nanopolish = config["progs"]["nanopolish"]

GENOME_VERSION     = config["ref"]["Genome_version"]
RmdReportScript    = os.path.join(config["scripts"]["script_folder"],"final_report","Nanopore_report.Rmd")

# Rscripts that convert table tsv's into GRL objs and then combine them:
Rmain_tsv2GRL     = os.path.join( config["scripts"]["script_folder"], config["scripts"]["Rmain_tsv2GRconv"] )
Rmain_combine_readchunks = os.path.join( config["scripts"]["script_folder"], config["scripts"]["Rmain_combine_read_chunks"] )
# Rscript with func defn's for the above
Rfuncs_tsv2GRL    = os.path.join( config["scripts"]["script_folder"], config[ "scripts"]["Rfuncs_tsv2GRconv"] )

# Rscript to create histograms of currents according to kmer
R_build_histlist_main   = os.path.join( config["scripts"]["script_folder"], config["scripts"]["Rmain_build_histlist"] )

# Rscript to overlap reads with RoIs
Rmain_overlap_reads_RsoI = os.path.join( config["scripts"]["script_folder"], config["scripts"]["Rmain_overlap_reads_RsoI"] )

#------------------------------------------------------
# --- Define output directories

SUBDIR_SYMLINKS           = "01_symlinks_fastqindex_chunks/"
SUBDIR_ALIGNED_MINIMAP    = "02_MM_aligned_chunks/"
SUBDIR_FILTERED_MINIMAP   = "03_MM_filtered_chunks/"
SUBDIR_SORTED_MINIMAPPED  = "04_MM_sortedbam/"
SUBDIR_EVENTALIGN         = "05_eventalign/"
SUBDIR_GR                 = "06_GRobjects/"
SUBDIR_GRproc             = "07_GRprocessing/"
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

   # TODO: implement sample-dependent targets with defaults.
   if ( config["execution"]["target_out"] == "report" ):
      OUTPUT_FILES.extend(
                          [ os.path.join( config["PATHOUT"], config["samplelist"][sampleLoopi_targets]["sampledir"], SUBDIR_REPORT, "" + sampleLoopi_targets + "_report.html") ]
                          )
   elif ( config["execution"]["target_out"] == "histlist" ):
      OUTPUT_FILES.extend(
                          [ os.path.join( config["PATHOUT"], config["samplelist"][sampleLoopi_targets]["sampledir"], SUBDIR_GR, sampleLoopi_targets + "_kmer_histlist.rds") ]
                          )
   elif ( config["execution"]["target_out"] == "ROI_olap" ):
      OUTPUT_FILES.extend(
                          [ os.path.join( config["PATHOUT"], config["samplelist"][sampleLoopi_targets]["sampledir"], SUBDIR_GRproc, sampleLoopi_targets + "_read_ROIolap_"+ region +".rds") for region in config["ref"]["RsoI"]  ]
                          )
   elif ( config["execution"]["target_out"] == "reads_GRL" ):
      OUTPUT_FILES.extend(
                          [ os.path.join( config["PATHOUT"], config["samplelist"][sampleLoopi_targets]["sampledir"], SUBDIR_GR, sampleLoopi_targets + "_reads_GRL.rds") ]
                          )
   elif ( config["execution"]["target_out"] == "mergedbam" ):
      OUTPUT_FILES.extend(
                          [ os.path.join( config["PATHOUT"], config["samplelist"][sampleLoopi_targets]["sampledir"], SUBDIR_SORTED_MINIMAPPED, sampleLoopi_targets + ".sorted.bam") ]
                         )
   elif ( config["execution"]["target_out"] == "aligned_chunks"):
      OUTPUT_FILES.extend(
                          get_chunkfiles( sampleLoopi_targets, os.path.join( config["PATHOUT"], config["samplelist"][sampleLoopi_targets]["sampledir"], SUBDIR_ALIGNED_MINIMAP) ,  "", ".sam", 0 )
                         )
   else:
      print("Unrecognized target output file format: ", config["execution"]["target_out"], " ... Terminating.")
      exit(1)

#---  DEBUGGING:
#------------------------------------------------------
# print("target out = " + config["execution"]["target_out"])
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
rule compile_report:
    # build the final output report in html format
    input:
        aligned_reads_bam = os.path.join( config["PATHOUT"], "{wc_sampleDir}", SUBDIR_SORTED_MINIMAPPED, "{wc_sampleName}.sorted.bam"),
        GRLreads          = os.path.join( config["PATHOUT"], "{wc_sampleDir}", SUBDIR_GR, "{wc_sampleName}_reads_GRL.rds")
    output:
        os.path.join( config["PATHOUT"], "{wc_sampleDir}", SUBDIR_REPORT, "{wc_sampleName}_report.html")
    params:
        " readcov_THRESH = 10;   ",
        " yplotmax = 10000; "
    log:
        logfile = os.path.join( config["PATHOUT"], "{wc_sampleDir}", SUBDIR_REPORT, "final_report_{wc_sampleName}.log")
    message:
        fmt("Compiling final report")
    shell:
        " Rscript -e  '{params} "
        " fin_readalignment_bam = \"{input.aligned_reads_bam}\"; "
        " fin_RsoI              = \"{input.transcriptome}\";"
        " fin_GRLreads          = \"{input.GRLreads}\";"
        " Genome_version        = \"{GENOME_VERSION}\"; "
        " rmarkdown::render(\"{RmdReportScript}\", output_file = \"{output}\" ) ' "

#------------------------------------------------------
rule bin_kmer_histlist:
    # Take a read-separated list of events, and assemble a histogram of
    # current values for each unique kmer observed.
    input:
        RDS_GRLreads       = os.path.join( config["PATHOUT"], "{wc_sampleDir}", SUBDIR_GR, "{wckmerhist_samplename}_reads_GRL.rds")
    output:
        RDS_histlist       = os.path.join( config["PATHOUT"], "{wc_sampleDir}", SUBDIR_GR, "{wckmerhist_samplename}_kmer_histlist.rds")
    params:
        current_histmin=config["execution"]["currenthist_minrange"],
        current_histmax=config["execution"]["currenthist_maxrange"],
        current_histres=config["execution"]["currenthist_res"],
        RDS_GRLreads_in=os.path.join( config["PATHOUT"], "{wc_sampleDir}", SUBDIR_GR, "{wckmerhist_samplename}_reads_GRL.rds"),
        RDS_histlist_out=os.path.join( config["PATHOUT"], "{wc_sampleDir}", SUBDIR_GR, "{wckmerhist_samplename}_kmer_histlist.rds"),
        k=5
    log:
        logfile=os.path.join( config["PATHOUT"], "{wc_sampleDir}", SUBDIR_GR, "{wckmerhist_samplename}_histlist.log")
    message:
        fmt("Build list of histograms for unique kmers in dataset.")
    shell:
        nice('Rscript', [ R_build_histlist_main,
                          "--current_histmin={params.current_histmin}",
                          "--current_histmax={params.current_histmax}",
                          "--current_histres={params.current_histres}",
                          "--rds_fin_readdat={params.RDS_GRLreads_in}",
                          "--rds_fout_histlist={params.RDS_histlist_out}",
                          "--k={params.k}",
                          "--logfile={log.logfile}",] )

#------------------------------------------------------

rule overlap_reads_w_RsoI:
    # Process the aligned reads and filter for only those
    # that overlap with the regions of interest.
    input:
        reads_in          = os.path.join( config["PATHOUT"], "{wc_sampleDir}", SUBDIR_GR, "{wc_sampleName}_reads_GRL.rds"),
        Refregion_in      = lambda wc: os.path.join( config["ref"]["RsoI_abspath"], config["ref"]["RsoI"][wc.wcReadROI_regionName])
    output:
        readROI_olaps     = os.path.join( config["PATHOUT"],  "{wc_sampleDir}",  SUBDIR_GRproc, "{wc_sampleName}_read_ROIolap_{wcReadROI_regionName}.rds" )
    params:
        sampleName        = "{wc_sampleName}",
        regionName        = "{wcReadROI_regionName}"
    log:
        logFile           = os.path.join( config["PATHOUT"], "{wc_sampleDir}", SUBDIR_GRproc, "{wc_sampleName}_{wcReadROI_regionName}_read_ROI_olap.log")
    message:
        fmt("Overlap reads from {input.reads_in} with region of interest {input.Refregion_in}.")
    shell:
        nice('Rscript', [ Rmain_overlap_reads_RsoI,
                          "--pathin_reads={input.reads_in}",
                          "--pathin_RsoI={input.Refregion_in}",
                          "--pathout_alignedreads={output.readROI_olaps}",
                          "--sampleName={params.sampleName}",
                          "--regionName={params.regionName}",
                          "--logFile={log.logFile}"] )

#------------------------------------------------------

rule minimizer:
    # Create indexed version of reference genome for fast
    # alignment with minimap2 later:
    input:
        refgenome_fasta  = os.path.join(DIR_REFGENOME , config['ref']['Genome_version']+ ".fa" )
    output:
        refgenome_mmiref = os.path.join(DIR_REFGENOME , config['ref']['Genome_version']+ ".mmi")
    params:
        options = " -d  "
    log:
        os.path.join( DIR_REFGENOME, config['ref']['Genome_version'], "_mmi2_minimizer_creation.log")
    message:
        fmt("Creating minimizer index of reference genome for minimap2.")
    shell:
        "{MM2} {params.options}  {output} {input} 2> {log}"

