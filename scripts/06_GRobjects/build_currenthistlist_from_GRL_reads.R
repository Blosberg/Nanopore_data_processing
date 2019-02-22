#================================================================================
# Rsript to take an rds object with read-separated current events, and generate a list of histograms
# of current observations for each unique kmer observed.
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
         --rds_fin_readdat = the input file with an rds object containing the read-partitioned event list.
         --rds_fout_histlist         = output file for the histogram list to be stored.
         --logFile               = file to print the logs to
         --k                     = number of bases in sequences.
         --current_histmin       = min value considered in the current histogram,
         --current_histmax       = max   .... ^^ ,
         --current_histres       = resolution of ... ^^

         --help              - print this text

         Example:
         ./test.R --arg1=1 --arg2='output.txt' --arg3=TRUE \n\n")

     q(save="no")
   }

   ## Parse arguments (we expect the form --arg=value)
   parseArgs <- function(x) strsplit(sub("^--", "", x), "=")

   argsDF <- as.data.frame(do.call("rbind", parseArgs(args)))
   argsL <- as.list(as.character(argsDF$V2))

   names(argsL) <- argsDF$V1

   # ==================================================================
   # --- return normalized histogram for a given set of breaks and GR
   get_normalized_current_hist <- function( GR_kmer_in      = stop("kmer-split event GR obj must be provided"),
                                            current_histmin = 50,
                                            current_histmax = 150,
                                            current_histres = 0.5
                                            )
   {
     breakset = seq( current_histmin, current_histmax, current_histres)

     result <- hist( GR_kmer_in$event_mean,
                     breaks  = breakset,
                     plot    = FALSE
                     )
     return(result)
   }
   # ==================================================================

  # read in the rds object of all experimental events split by read.
  readdat            <- readRDS( argsL$rds_fin_readdat)
  GRLin_splitby_read <- readdat$Events_GRL_splitbyread

  # undo read-partitioning to have all events together
  GR_all_events  <- unlist( GRLin_splitby_read )

  # now split by model_kmer:
  GRL_splitby_kmer     = split( GR_all_events,    GR_all_events$model_kmer     )

  # filter out "NNNNN":
  GRL_splitby_kmer     = GRL_splitby_kmer[ which ( names(GRL_splitby_kmer) != "NNNNN" )   ]

  # Now take this kmer-separated list, and generated a histogram of current observations for each kmer:
  allkmer_histlist     <- lapply( GRL_splitby_kmer, function(x)   get_normalized_current_hist ( GR_kmer_in      = x,
                                                                                                current_histmin = as.numeric( argsL$current_histmin),
                                                                                                current_histmax = as.numeric( argsL$current_histmax),
                                                                                                current_histres = as.numeric( argsL$current_histres)
                                                                                                )
                                  )

 # Now export this list of histograms to an rds output file for easy manipulation later:
 saveRDS( allkmer_histlist,
          file = argsL$rds_fout_histlist )


 write( x = paste0( " build_currenthistlist complete. ",
                   as.character(length(allkmer_histlist)),
                   " histograms written to file ",
                   argsL$rds_fout_histlist),
        file = argsL$logfile,
        append = TRUE
        )
