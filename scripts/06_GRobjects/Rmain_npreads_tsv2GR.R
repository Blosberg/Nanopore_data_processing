# Rmain_npreads_tsv2GR
# import nanopore data in .tsv format and output read-separated GRanges objects:
# For optimization, model_mean/model_stdv columns are omitted to conserve space.

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
      Rfuncs_tsv2GRconv --script with function definitions used here.
      output_reads_GRL    -- rds filename read-separated GRL output.
      output_poremodel    -- rds filename with model date (for each model kmer) saved.
      sampleName          -- name of sample for documentation
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

# Run Functions -----------------------------------------------------------

# e.g. (replace this list with actual arguments)
Rfuncs_tsv2GRconv  <- argsL$Rfuncs_tsv2GRconv
output_reads_GRL   <- argsL$output_reads_GRL
output_poremodel   <- argsL$output_poremodel
sampleName         <- argsL$sampleName
Flatten_reads      <- argsL$Flatten_reads
logFile            <- argsL$logFile
Ealign_files       <- unlist( strsplit(argsL$Ealign_files,",")  )

#===============================================================
suppressPackageStartupMessages( library(GenomicRanges) )
suppressPackageStartupMessages( library(dplyr)         )
source(Rfuncs_tsv2GRconv)

table_dat       = get_event_dat( Event_file_list = Ealign_files,
                                 logFile         = logFile )

# catalogue reads by index, not by name, since
# chimeric alignments are treated as distinct reads
# readname preserves this info
read_list_final = unique( table_dat$read_index  )
Nreads          = length( read_list_final     )

if( !( identical( as.numeric(c(1:Nreads)) , as.numeric(read_list_final) )  ))
  {  writeLines("ERROR: final read list is not step-wise increasing", logFile );
     exit(1)
  }

# ========================================================================================
# output the "poremodel" data - so we have the model data saved
if ( ! is.null( output_poremodel ))
  { writeout_pore_model( table_dat, output_poremodel, strand_type="RNA" ) }

# ========================================================================================
# remove non-finite entries:

# table_dat = table_dat [ which (! is.na (table_dat$event_level_mean) ), ]

# =======================
# Add strand information:

table_dat = assign_strand( Datin                 = table_dat,
                           perform_sanity_checks = TRUE )

# ================================================
# Split by read

ReadList_finite_stranded  <- split( table_dat,
                                    table_dat$read_index )

rm( table_dat )
MemLog = gc( verbose = FALSE )

# ================================================
# Convert table to GRanges list; flatten overlapping read-events iff Flatten_reads == TRUE

GRL_out <- convert_tbl_readlist_to_GRL( Readlist_in   = ReadList_finite_stranded,
                                        Flatten_reads = Flatten_reads )
rm( ReadList_finite_stranded )
MemLog = gc( verbose = FALSE )

# ================================================
# BUILD GRanges OBJECT TO Process

# output the GRL object to store reads:
saveRDS( list( "sampleName"             = sampleName,
               "Flattened_reads"        = Flatten_reads,
               "Events_GRL_splitbyread" = GRL_out ),
         file = output_reads_GRL  )

writeLines("script tables2GR complete: read list and model catalogue printed to output files. Now exiting normally.", logFile )
