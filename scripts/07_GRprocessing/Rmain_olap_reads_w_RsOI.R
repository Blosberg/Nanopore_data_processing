# Compare observed reads with known reference locations of interest.
# produce an object containing only the reads that overlap with these loci

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
      pathin_reads --path to Granges List structure with reads aligned to reference genome,
      pathin_RsoI  --path to Granges list of regions of interest in the current study,
      pathout_alignedreads  -- path to which the output of this script should be directed
      sampleName            -- self-explanatory,
      regionName            -- name describing the reference region under consideration,
      logFile               -- self-explanatory

      Example:
      $Rscript ./Rmain_olap_reads_w_RsoI.R --arg1=1 --arg2='output.txt' --arg3=TRUE \n\n")

  q(save="no")
}

## Parse arguments (we expect the form --arg=value)
parseArgs <- function(x) strsplit(sub("^--", "", x), "=")
argsDF    <- as.data.frame(do.call("rbind", parseArgs(args)))
argsL     <- as.list(as.character(argsDF$V2))

names(argsL) <- argsDF$V1

suppressPackageStartupMessages( library(GenomicRanges) )
suppressPackageStartupMessages( library(dplyr)         )

# import the nanopore reads as GRL
Readstruct_all_in      <- readRDS ( argsL$pathin_reads )
reads                  <- Readstruct_all_in$Events_GRL_splitbyread

# Define the regions of interest from the reference genome:
# They will come in sets defined by Region_groups
if( ! file.exists( argsL$pathin_RsoI) )
  {
  cat( "NO loci information found. An empty olap file will be produced and this section will be ignored in the final report.",
       file = argsL$logFile,
       append = TRUE,
       sep = "\n" )

  saveRDS( NULL, argsL$pathout_alignedreads )
  q(save="no")
}else{
  RsoI_in           <- readRDS( argsL$pathin_RsoI )
}

cat( "Obtained ROI file; processing overlaps now.",
     file   = argsL$logFile,
     append = TRUE,
     sep    = "\n" )

output               = list()
output$sampleName    = argsL$sampleName
output$RefRegionName = argsL$regionName

OLAP_skip_TOL = 3
# ========================================================
# count how many groupings of loci we are considering.
# (groupings typically lump types of modifications in different regions)

# expand the ROIs slightly to ensure overlap is recorded, even if
# the read skips a base at the exact position of interest.

# This oddly-constructed boolean is in case one of the RsoI is actually a whole gene,
# and the list is a set exons. the width will then return a set of names.
# (in this case, we don't want to bother with expanding)
if ( identical( all( width ( RsoI_in$loci )  < OLAP_skip_TOL ) , TRUE ) )
  {
  start( RsoI_in$loci ) <-  ( start( RsoI_in$loci ) - OLAP_skip_TOL )
  end(   RsoI_in$loci ) <-  ( end(   RsoI_in$loci ) + OLAP_skip_TOL )

  cat( paste( "Enlarged ROI region for group", argsL$regionName ),
       file   = argsL$logFile,
       append = TRUE,
       sep    = "\n" )
  }


# Build aligned_reads list for alignment of reads to each subset of loci:
# and collect indices of reads that hit at least one ROI
read_indices_on_ROI   = queryHits ( findOverlaps( reads,
                                                  RsoI_in$loci
                                                 )
                                    )
# filter for uniqueness (we don't care how many ROIs the read overlaps with):
read_indices_on_ROI = unique( read_indices_on_ROI )

# Collect just those reads, and put them together as an output object.
output$aligned_reads = reads[ read_indices_on_ROI ]

output$sampleName <- argsL$sampleName
output$regionName <- argsL$regionName

# ======  SAVE OUTPUT ===========

cat( paste( "Storing output for sample:",
                   argsL$sampleName,
                   ", region:", argsL$regionName ),
     file   = argsL$logFile,
     append = TRUE,
     sep    = "\n" )


saveRDS( output, file = argsL$pathout_alignedreads )

cat( "Saved RDS file. Program complete.",
      file   = argsL$logFile,
      append = TRUE,
      sep    = "\n" )
