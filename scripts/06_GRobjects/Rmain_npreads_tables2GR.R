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
      output_reads_GRL  -- rds filename read-separated GRL output.
      output_poremodel  -- rds filename with model date (for each model kmer) saved.
      samplename   -- name of sample for documentation
      Ealign_files -- array of files to use as input
      logFile      -- filename to pipe output to

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

dat_all         = get_event_dat( Ealign_files )

read_list_final = unique( dat_all$read_index  )
Nreads          = length( read_list_final     )

if( !( identical( as.numeric(c(1:Nreads)) , as.numeric(read_list_final) )  ))
  {  writeLines("ERROR: final read list is not step-wise increasing", logFile );
     exit(1)

  }


# ========================================================================================
# remove non-finite entries:

dat_finite = dat_all [ which (! is.na (dat_all$event_level_mean) ), ]

# =======================
# Add strand information:
dat_finite_stranded = assign_strand( dat_finite )

# Check that they are all attributed correctly
Watson = dat_finite_stranded[ dat_finite_stranded$strand =="+", ];
Crick  = dat_finite_stranded[ dat_finite_stranded$strand =="-", ];
Unk    = dat_finite_stranded[ dat_finite_stranded$strand =="*", ];

if( dim(Unk)[1] > 0 )
  { writeLines("ERROR: reads being assigned unknown strand", logFile );
    stop(paste("ERROR: reads being assigned unknown strand")) }

# ================================================
# BUILD GRanges OBJECT TO Process

allevents_GR = GRanges( seqnames     = dat_finite_stranded$contig,
                        strand       = dat_finite_stranded$strand,
                        range        = IRanges( start = dat_finite_stranded$position+1,
                                                end   = dat_finite_stranded$position+5 ), # --- nopolish provides position before beginning of 5bp window.

#                        reference_km = dat_finite_stranded$reference_kmer,
                        read_index   = dat_finite_stranded$read_index,
                        event_index  = dat_finite_stranded$event_index,
                        event_mean   = dat_finite_stranded$event_level_mean,
                        event_stdv   = dat_finite_stranded$event_stdv,
                        event_length = dat_finite_stranded$event_length,
                        model_kmer   = dat_finite_stranded$model_kmer,
                        model_mean   = dat_finite_stranded$model_mean,
                        model_stdv   = dat_finite_stranded$model_stdv
)

# output the "poremodel" data - so we have the model data saved
model_GRL = split( allevents_GR,
                   allevents_GR$model_kmer )

model_mean_list <- lapply( model_GRL, function(x) unique(x$model_mean) )
model_stdv_list <- lapply( model_GRL, function(x) unique(x$model_stdv) )

test_length     <- unlist( lapply( 1:length(model_mean_list), function(x) length(model_mean_list[[x]]) + length(model_stdv_list[[x]]) ) )

if( min(test_length != 2 ) || max (test_length !=2 ))
  { stop("ERROR: irregular model values found in table data. Terminating.") }

model_dat <- cbind( unlist(model_mean_list), unlist(model_stdv_list))
colnames(model_dat) <- c("mean", "stdv")

model_dat <- model_dat[ rownames(model_dat) != "NNNNN", ]

if( data_is_RNA )
  { rownames(model_dat) <- gsub("T", "U", rownames(model_dat) ) }

write.table( model_dat,
	     file = output_poremodel,
	     sep="\t")

# ------------------------------------------
# split events up into list separated by read.

allevents_GR$model_mean = NULL;
allevents_GR$model_stdv = NULL;


reads_GRL = split( allevents_GR,
                   allevents_GR$read_index )


# output the GRL object to store reads:
saveRDS( list( "samplename"             = samplename,
               "Events_GRL_splitbyread" = reads_GRL),
         file = output_reads_GRL  )


writeLines("script tables2GR complete: read list and model catalogue printed to output files. Now exiting normally.", logFile )
