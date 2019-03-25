# Environment from this execution is stored in : "/home/bosberg/projects/nanopore/signal_processing_Rworkspace.RData"
# import nanopore data and run some statistics on them:
# rm(list=ls()) # CLEAN UP EVERYTHING


# @@@ TODO: generalize this as an input boolean (see issue #9 on GH):
data_is_RNA=TRUE

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
      Rfuncs_table2GRconv --script with function definitions used here.
      output_reads_GRL    -- rds filename read-separated GRL output.
      output_poremodel    -- rds filename with model date (for each model kmer) saved.
      samplename          -- name of sample for documentation
      Flatten_reads       -- Boolean: should we flatten the read events that overlap.
      logFile             -- filename to pipe output to
      Ealign_files        -- array of files to use as input

      Example:
      ./test.R --arg1=1 --arg2='output.txt' --arg3=TRUE \n\n")

  q(save="no")
}

## Parse arguments (we expect the form --arg=value)
parseArgs <- function(x) strsplit(sub("^--", "", x), "=")
argsDF    <- as.data.frame(do.call("rbind", parseArgs(args)))
argsL     <- as.list(as.character(argsDF$V2))

names(argsL) <- argsDF$V1

# catch output and messages into log file
# out <- file(argsL$logFile, open = "wt")
# sink(out, type = "output")
# sink(out, type = "message")

# Run Functions -----------------------------------------------------------

# e.g. (replace this list with actual arguments)
Rfuncs_table2GRconv  <- argsL$Rfuncs_table2GRconv
output_poremodel     <- argsL$output_poremodel
output_reads_GRL     <- argsL$output_reads_GRL
samplename           <- argsL$samplename
logFile              <- argsL$logFile
Ealign_files         <- unlist( strsplit(argsL$Ealign_files,",")  )

#===============================================================
suppressPackageStartupMessages( library(GenomicRanges) )
suppressPackageStartupMessages( library(dplyr)         )
source(Rfuncs_table2GRconv)

dat_all         = get_event_dat( Event_file_list = Ealign_files, 
                                 logFile         = logFile )

read_list_final = unique( dat_all$read_index  )
Nreads          = length( read_list_final     )

if( !( identical( as.numeric(c(1:Nreads)) , as.numeric(read_list_final) )  ))
  {  writeLines("ERROR: final read list is not step-wise increasing", logFile );
     exit(1)

  }

# ========================================================================================
# output the "poremodel" data - so we have the model data saved
if ( ! is.null( output_poremodel ))
  { writeout_pore_model( dat_all, output_poremodel ) }

# ========================================================================================
# remove non-finite entries:

dat_finite = dat_all [ which (! is.na (dat_all$event_level_mean) ), ]

# =======================
# Add strand information:

dat_finite_stranded = assign_strand( dat_finite, perform_sanity_checks = TRUE )

# ================================================
# Split by read

ReadList_finite_stranded  <- split( dat_finite_stranded, 
                                    dat_finite_stranded$read_index )

# ================================================
# Flatten overlapping read-events:

if(Flatten_reads)
{
  # tic <- Sys.time()
  Flatreadlist <-
    lapply(ReadList_finite_stranded, function(x)
      flatten_read_tbl (read_tbl_in   = x,
                        perform_sanity_checks = FALSE))
  # toc <- Sys.time()
  
  GRL_out <- GRangesList(lapply(Flatreadlist, function(x)
    GRanges(
      seqnames = x$contig,
      strand   = x$strand,
      range    = IRanges(start = x$position + 1,
                         end   = x$position + 5),
      read_index   = x$read_index,
      event_index  = x$event_index,
      event_mean   = x$event_level_mean,
      event_stdv   = x$event_stdv,
      event_length = x$event_length,
      model_kmer   = x$model_kmer
    )))
} else {

  GRL_out <- GRangesList( lapply( ReadList_finite_stranded, function(x)
    GRanges(
      seqnames = x$contig,
      strand   = x$strand,
      range    = IRanges(start = x$position + 1,
                         end   = x$position + 5),
      read_index   = x$read_index,
      event_index  = x$event_index,
      event_mean   = x$event_level_mean,
      event_stdv   = x$event_stdv,
      event_length = x$event_length,
      model_kmer   = x$model_kmer
    )))
  
}

# ================================================
# BUILD GRanges OBJECT TO Process

# output the GRL object to store reads:
saveRDS( list( "samplename"             = samplename,
               "Events_GRL_splitbyread" = GRL_out ),
         file = output_reads_GRL  )

writeLines("script tables2GR complete: read list and model catalogue printed to output files. Now exiting normally.", logFile )
