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

cat( "Reading in reads, refGenome and poremodel data.",
     file = argsL$logFile,
     append = TRUE,
     sep = "\n" )

ref_Genome    <- readDNAStringSet( argsL$pathin_refGen )
readdat       <- readRDS( argsL$pathin_reads ) # these are just all the reads with at least one overlap on a ROI:
                                               # no structural organization in place (as some reads might overlap multiple RsOI)
putloci       <- readRDS( argsL$pathin_RsoI  )
poremodel     <- read.table( file = argsL$poremodel_ref,
                             sep = "\t",
                             header = TRUE
                             )

# ========================================

if( typeof( width( putloci$loci ) ) == "S4" )
  {
  # in this case, loci are full genes and we should not
  # attempt to produce "plotdat" over them.
  sampleROI_dat <- NULL
  cat( "width(loci) is in S4 format, indicating gene list, rather than specific loci. sampleROI_dat will be exported as NULL. Exiting.",
     file = argsL$logFile,
     append = TRUE,
     sep = "\n" )
  q(save="no")

  }

OLAP_skip_TOL = 3

cat( "Filtering loci for coverage.",
     file = argsL$logFile,
     append = TRUE,
     sep = "\n" )

#
# reduce list of RsoI to those with sufficient coverage.
loci_filtered_for_coverage <- filter_loci_for_coverage (  loci   = putloci$loci,
                                                          reads  = readdat$aligned_reads,
                                                          mincov = argsL$mincov_in,
                                                          OLAP_skip_TOL = OLAP_skip_TOL
                                                       )

if( length( loci_filtered_for_coverage ) == 0 )
  {
  # in this case, none of the loci had sufficient coverage, so the plotdat is empty.
  sampleROI_dat <- NULL
  cat( "length of sufficiently-covered loci = 0. With no coverage on any ROI, sampleROI_dat will be exported as NULL. Exiting.",
     file = argsL$logFile,
     append = TRUE,
     sep = "\n" )
  q(save="no")
  }


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
cat( "Finding ROI overlaps.",
     file = argsL$logFile,
     append = TRUE,
    sep = "\n" )

overlaps = findOverlaps(  loci_filtered_for_coverage_expanded,
                      readdat$aligned_reads )

# Don't need this anymore. The loci should refer directly to the exact ROI
rm(loci_filtered_for_coverage_expanded)

# overlaps_by_group queryHits now references the indices of COVERED loci.

# TODO: add row/col names for the output data structures.
cat( "collecting sampleROI_dat.",
     file = argsL$logFile,
     append = TRUE,
     sep = "\n" )

sampleROI_dat <- collect_sample_dat_over_ROIs(   ref_Genome      = ref_Genome,
                                                 aligned_reads   = readdat$aligned_reads,
                                                 ROI_overlaps    = overlaps,
                                                 loci_covered    = loci_filtered_for_coverage,
                                                 poremodel_in    = poremodel,
                                                 plotrange_in    = 10
                                                )

coverage   <- unlist( lapply( c(1:length(sampleROI_dat)),
                            function(ROI_index) dim( sampleROI_dat[[ROI_index]]$read_normdiff )[1] ))
coverage_order_list <- order( coverage,
                              decreasing = TRUE )

sampleROI_dat_sorted_by_coverage = sampleROI_dat[ coverage_order_list ]


plotdat_out = list( "sampleName" = argsL$sampleName,
                    "regionName" = argsL$regionName,
                    "ROIdat"     = sampleROI_dat_sorted_by_coverage
                    )

# ========================================

# -------------------------------------------------------------------------
cat( "Saving output.",
     file = argsL$logFile,
     append = TRUE,
     sep = "\n" )

saveRDS( object = plotdat_out,
         file   = argsL$pathout_ROIplotdat )

cat( paste( "Finished exporting RDS to output file:",
             argsL$pathout_ROIplotdat,
             " ... Exiting." ),
     file = argsL$logFile,
     append = TRUE,
     sep = "\n" )

