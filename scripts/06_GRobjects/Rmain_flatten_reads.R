#================================================================================
# Rsript to flatten reads from an rds object with read-separated current events
# B. Osberg, MDC Berlin 2019
#================================================================================

suppressPackageStartupMessages( library(GenomicRanges) )

args <- commandArgs(trailingOnly=TRUE)

   ## Default setting when no arguments passed
   if(length(args) < 1) {
     args <- c("--help")
   }

   ## Help section
   if("--help" %in% args) {
     cat("
         Render to report

         Arguments:

	 --Rfuncs_flattenreads = script of necessary functions

	 --readsGRL_in         = the input file with an rds object containing the read-partitioned event list.

	 --flatreadsGRL_out = output file for the histogram list to be stored.

	 --logFile             = file to print the logs to

         --help               ==>  print this text

         Example:
         ./test.R --arg1=1 --arg2='output.txt' --arg3=TRUE \n\n")

     q(save="no")
   }


   ## Parse arguments (we expect the form --arg=value)
   parseArgs <- function(x) strsplit(sub("^--", "", x), "=")

   argsDF <- as.data.frame(do.call("rbind", parseArgs(args)))
   argsL  <- as.list(as.character(argsDF$V2))

   names(argsL) <- argsDF$V1

 # argsL <- list (
# "Rfuncs_flattenreads" = "/home/bosberg/projects/nanopore/scripts/06_GRobjects/Rfuncs_flatten_reads.R",
# "readsGRL_in"         = "/scratch/AG_Akalin/bosberg/nanopore/pipeline_output/testset/06_GRobjects/TESTSET0_reads_GRL.rds",
# "flatreadsGRL_out"    = "/scratch/AG_Akalin/bosberg/nanopore/pipeline_output/testset/06_GRobjects/TESTSET0_reads_flat_GRL.rds",
# "logFile"             = "/scratch/AG_Akalin/bosberg/nanopore/pipeline_output/testset/06_GRobjects/TESTSET0_flattenreads_GRL.log",
# "samplename"          = "TESTSET0" )

# ==================================================================

source( argsL$Rfuncs_flattenreads )

# read in the rds object of all experimental events split by read.
  readdat            <- readRDS( argsL$readsGRL_in )
  GRLin_splitby_read <- readdat$Events_GRL_splitbyread

  reads_flattened_GRL <- lapply( GRLin_splitby_read, function(r) flatten_read( read_GR_in = r) )

 # Now export this output:
 saveRDS( reads_flattened_GRL,
          file = argsL$flatreadsGRL_out )

 write( x = paste0( " flatten_reads_main complete. ",
                   as.character(length(reads_flattened_GRL)),
                   " histograms written to file ",
                   argsL$rds_fout_flattened ),
        file   = argsL$logFile,
        append = TRUE
        )
