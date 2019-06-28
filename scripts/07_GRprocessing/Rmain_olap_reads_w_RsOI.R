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
  writeLines("NO loci information found. An empty olap file will be produced and this section will be ignored in the final report.", argsL$logFile )
  saveRDS( NULL, argsL$pathout_alignedreads )
  q(save="no")
}

  writeLines( "Obtained ROI file; processing overlaps now.",
              argsL$logFile )


RsoI_in           <- readRDS(  argsL$pathin_RsoI )

output               = list()
output$aligned_reads = list()
output$N_g_filtered  = list()
output$sampleName    = argsL$sampleName
output$RefRegionName = argsL$regionName


# ========================================================
# count how many groupings of loci we are considering.
# (groupings typically lump types of modifications in different regions)
N_locus_groupings = length( length( RsoI_in$Region_groups ) )

# list the indices of reads that cover at least one ROI:
read_indices_on_ROI = list()

# Build aligned_reads list for alignment of reads to each subset of loci:
for ( group in  names( RsoI_in$Region_groups )  )
  {
  read_indices_on_ROI[[group]]   = queryHits ( findOverlaps( reads,
                                                                    RsoI_in$Region_groups[[group]]
                                                                   )
                                                     )
  # collect indices of reads that hit at least one ROI:
  read_indices_on_ROI[[ group ]] = unique( read_indices_on_ROI[[ group ]] )

  # filter for uniqueness:
  output$aligned_reads[[group]] = reads[ read_indices_on_ROI[[ group ]] ]

  # store the length.
  output$N_g_filtered[[group]]   = length( RsoI_in$Region_groups[[group]] )
}

names( output$aligned_reads ) <- names( RsoI_in$Region_groups )
names( output$N_g_filtered  ) <- names( RsoI_in$Region_groups )


# ======  SAVE OUTPUT ===========

saveRDS( output, file = argsL$pathout_alignedreads )
