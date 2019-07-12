# Plot current from reads along RsoI:

## Collect arguments
args <- commandArgs(TRUE)

## Default setting when no arguments passed
if(length(args) < 1) {
  args <- c("--help")
}

## Help section
if("--help" %in% args) {
  cat("
      Render to report

      Arguments:
      pathin_reads         --path to Granges List structure with reads aligned to reference genome,
      pathin_RsoI          --path to Granges list of regions of interest in the current study,
      pathout_ROIdat       --where should we send the plots

      logFile              -- self-explanatory

      Example:
      $Rscript ./Rmain_olap_reads_w_RsoI.R --arg1=1 --arg2='output.txt' --arg3=TRUE \n\n")

  q(save="no")
}

## Parse arguments (we expect the form --arg=value)
parseArgs <- function(x) strsplit(sub("^--", "", x), "=")
argsDF    <- as.data.frame(do.call("rbind", parseArgs(args)))
argsL     <- as.list(as.character(argsDF$V2))

names(argsL) <- argsDF$V1
# BiocManager::install("Biostrings")
# BiocManager::install("BSgenome.Hsapiens.UCSC.hg19")
suppressPackageStartupMessages( library(GenomicRanges) )
suppressPackageStartupMessages( library(Biostrings) )


# ======= DEBUGGING: DELETE THIS: ========
# PATHOUT="/scratch/AG_Akalin/bosberg/nanopore/pipeline_output/"
# sampleDir = "20180417_1233_HEK293_polyA_RNA/"
# # sampleDir = "20190702_1647_IVdBm6A_FAK57816_38717116/"
# sampleName="HEK_untreated"
# #sampleName = "IVdBm6A"
# ROIname="MihaM_2pO"
# # "barcodes", "homoP"
# pathin_RsoI = paste0( "/scratch/AG_Akalin/bosberg/nanopore/ref/Regions_of_interest/", "MihaM_2po_Ill", ".rds" )
# #pathin_RsoI = paste0( "scratch/AG_Akalin/bosberg/nanopore/ref/Regions_of_interest/IV_dB/", choose one: [IVdB_ROI_barcodes.rds, IV_dB/IVdB_ROI_homoP.rds] )
#
# argsL=list(
#
# "pathin_reads"   = paste0( PATHOUT, sampleDir, "/07_GRprocessing/", sampleName, "_read_ROIolap_", ROIname, ".rds"),
# "pathin_RsoI"    = pathin_RsoI,
#
# "path_funcdefs"  = "/clusterhome/bosberg/projects/nanopiper/scripts/08_GRdata_vis/Rfuncs_plot_currents_ROI.R",
# "pathin_refGen"  = "/fast/AG_Akalin/refGenomes/hg19_canon/hg19_canon.fa",                         # path to Granges list of regions of interest in the current study
# # "pathin_refGen"  = "/scratch/AG_Akalin/refGenomes/IV_DeBRuijn_barcoded/IV_dB_Bcoded.fa",
#
# "pathout_ROIplotdat" = paste0( PATHOUT, sampleDir, "08_plotdat_vis/", sampleName, "_olap_", ROIname, "_plotdat.rds"),
#
# "sampleName"     = sampleName,
# "poremodel_ref"  = "/clusterhome/bosberg/projects/nanopiper/dev/ref/poremodel_RNA.csv",
#
# "mincov_in"      = 10,
# "plotrange_in"   = 10,
#
# "logFile"        = paste0( PATHOUT, "20190702_1647_IVdBm6A_FAK57816_38717116/08_plotdat_vis/homoP.log")
# #"logFile"        = paste0( PATHOUT, "20180417_1233_HEK293_polyA_RNA/08_testing/test.log"                             # self-explanatory
# )
#
# ======= DEBUGGING: DOWN TO HERE: ========

source( argsL$path_funcdefs )

# overlaps CONVENTION: findoverlaps ( region of interest, reads        )
#                                      < queryHits >    | < subjectHits >
#                                     __________________|______________
#                                       locus_i         |   read_i
#                                       ...             |   ...

# ========================================
# Read in inputs:

writeLines( "Reading in reads, refGenome and poremodel data.",
            argsL$logFile )
ref_Genome    <- readDNAStringSet( argsL$pathin_refGen )
readdat       <- readRDS( argsL$pathin_reads ) # these are just all the reads with at least one overlap on a ROI:
                                               # no structural organization in place (as some reads might overlap multiple RsOI)
putloci       <- readRDS( argsL$pathin_RsoI  )
poremodel     <- read.table( file = argsL$poremodel_ref,
                             sep = "\t",
                             header = TRUE
                             )
# ========================================
# Read in inputs:

OLAP_skip_TOL = 3

writeLines( "Filtering loci for coverage.",
            argsL$logFile )
#
# reduce list of RsoI to those with sufficient coverage.
loci_filtered_for_coverage <- filter_loci_for_coverage (  loci   = putloci$loci,
                                                          reads  = readdat$aligned_reads,
                                                          mincov = argsL$mincov_in,
                                                          OLAP_skip_TOL = OLAP_skip_TOL
                                                        )

# expand loci slightly to ensure overlap is registered,
# even if the exact base is skipped.
loci_filtered_for_coverage_expanded <- loci_filtered_for_coverage

if ( identical( all( width ( loci_filtered_for_coverage_expanded )  < OLAP_skip_TOL ) , TRUE ) )
  {
  start( loci_filtered_for_coverage_expanded ) <-  ( start( loci_filtered_for_coverage ) - OLAP_skip_TOL )
  end(   loci_filtered_for_coverage_expanded ) <-  ( end(   loci_filtered_for_coverage ) + OLAP_skip_TOL )
  }


# NOW check overlaps of reads with only the _covered_ loci.
#TODO: We run this overlaps more than once; if speedup is required later, eliminate this repetition.
writeLines( "Finding ROI overlaps.",
            argsL$logFile )
overlaps = findOverlaps(  loci_filtered_for_coverage_expanded,
                          readdat$aligned_reads )

# Don't need this anymore. The loci should refer directly to the exact ROI
rm(loci_filtered_for_coverage_expanded)

# overlaps_by_group queryHits now references the indices of COVERED loci.

# TODO: add row/col names for the output data structures.
writeLines( "collecting sampleROI_dat.",
            argsL$logFile )
sampleROI_dat<- collect_sample_dat_over_ROIs(       SampleName           = argsL$sampleName,
                                                    ref_Genome           = ref_Genome,
                                                    aligned_reads        = readdat$aligned_reads,
                                                    ROI_overlaps         = overlaps,
                                                    group_loci_covered   = loci_filtered_for_coverage,
                                                    poremodel_in         = poremodel,
                                                    plotrange_in         = 10 )

writeLines( "Saving output.",
            argsL$logFile )

saveRDS( object = sampleROI_dat,
         file   = argsL$pathout_ROIplotdat )

writeLines( "Finished exporting RDS file for plot data.",
            argsL$logFile )
