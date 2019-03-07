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
         --rds_fin_readdat     = the input file with an rds object containing the read-partitioned event list.
         --rds_fout_flattened  = output file for the histogram list to be stored.
         --logFile             = file to print the logs to
	 
	 --flattenscript_funcs = script of necessary functions

         --help               ==>  print this text

         Example:
         ./test.R --arg1=1 --arg2='output.txt' --arg3=TRUE \n\n")

     q(save="no")
   }

#  argsL=list( "rds_fin_readdat"     = "/clusterhome/bosberg/projects/nanopore/scripts/06_GRobjects/testreads_50_GRL.rds",
#              "rds_fout_flattened"  = "/clusterhome/bosberg/projects/nanopore/scripts/06_GRobjects/testreads_50_flattened_GRL.rds",
#              "flattenscript_funcs" = "/clusterhome/bosberg/projects/nanopore/scripts/06_GRobjects/flatten_reads_funcs.R",
#              "logfile"             = "/clusterhome/bosberg/projects/nanopore/scripts/06_GRobjects/testreads_50_flattening.log"     )

   ## Parse arguments (we expect the form --arg=value)
   parseArgs <- function(x) strsplit(sub("^--", "", x), "=")

   argsDF <- as.data.frame(do.call("rbind", parseArgs(args)))
   argsL  <- as.list(as.character(argsDF$V2))

   names(argsL) <- argsDF$V1

source( argsL$flattenscript_funcs )
   # ==================================================================

  # read in the rds object of all experimental events split by read.
  readdat            <- readRDS( argsL$rds_fin_readdat)
  GRLin_splitby_read <- readdat$Events_GRL_splitbyread

  start_time <- Sys.time()
  flatten_read( read_GR_in = testread_lean )
  end_time <- Sys.time()

  
  start_time <- Sys.time()
  flatten_read( read_GR_in = testread_fat )
  end_time <- Sys.time()

  
  start_time <- Sys.time()
  reads_flattened_GRL <- lapply( GRLin_splitby_read, function(r) flatten_read( read_GR_in = r) )
  end_time <- Sys.time()


 # Now export this output:
 saveRDS( reads_flattened_GRL,
          file = argsL$rds_fout_flattened )


 write( x = paste0( " flatten_reads_main complete. ",
                   as.character(length(reads_flattened_GRL)),
                   " histograms written to file ",
                   argsL$rds_fout_flattened ),
        file = argsL$logfile,
        append = TRUE
        )
