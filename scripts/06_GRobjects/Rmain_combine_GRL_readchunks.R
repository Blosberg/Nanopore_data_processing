# Rmain_combine_GRL_readchunks.R
# Take many lists of reads and combine them into one.

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
                         "--Rfuncs_tsv2GRconv="+Rfuncs_tsv2GRL,
                         "--GRL_reads_combined_out={output.GRL_reads_combined}",
                         "--output_poremodel={output.poremodel}",
                         "--sampleName={params.sampleName}",
                         "--logFile={log}",
                         "--poremodel_chunks={input.poremodel_chunks}",
                         "--GRL_chunk_files={input.GRL_chunk_files}"] )


      Arguments:
      Rfuncs_tsv2GRconv --script with function definitions used here.
      GRL_reads_combined_out  -- combined (across all chunks) rds filename read-separated GRL output.
      output_poremodel    -- rds filename with model date (for each model kmer) saved.
      samplename          -- name of sample for documentation
      logFile             -- filename to pipe output to
      poremodel_chunks    -- Array of unique kmer models from each chunk.
      GRL_chunk_files     -- array of files to use as input

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
Rfuncs_tsv2GRconv       <- argsL$Rfuncs_tsv2GRconv
GRL_reads_combined_out  <- argsL$GRL_reads_combined_out
# output_poremodel        <- argsL$output_poremodel
sampleName              <- argsL$sampleName
logFile                 <- argsL$logFile
GRL_chunk_files         <- unlist( strsplit(argsL$GRL_chunk_files,",")  )

#===============================================================
suppressPackageStartupMessages( library(GenomicRanges) )
source(Rfuncs_tsv2GRconv)

Nchunks         = length( readchunk_files )

Combined_reads  = readRDS( readchunk_files[1] )

if ( Nchunks >= 2 )
  {
  for ( i  in c(2:Nchunks) )
     {
     # count up the total current number of reads:
     Nreads  = length( Combined_reads$Events_GRL_splitbyread )

     # Read in the next chunk of data:
     tempdat = readRDS( readchunk_files[i] )

     # sanity-check that the names are in order:
     # TODO: comment this line out for efficiency once we're convinced this works.
     if( ! identical( as.character( unlist( lapply( tempdat_copy$Events_GRL_splitbyread, function(x) unique( x$read_index ) ) ) ), as.character( c(1:length( tempdat_copy$Events_GRL_splitbyread ) )  )  ) ){ stop(paste("Disordered names in reads of file: ", readchunk_files[i] ) ) }

     # offset the read index for each one of these:
     tempdat$Events_GRL_splitbyread  = lapply( tempdat$Events_GRL_splitbyread, function(x) offset_read_indices ( GRLchunk_in  = x, Nread_offset = Nreads ) )
    
     # Likewise offset the corresponding names: 
     names( tempdat$Events_GRL_splitbyread )  <-  as.character( as.numeric( names(tempdat$Events_GRL_splitbyread) ) + Nreads )
    
     # Now, absorb these reads into the "combined" data set
     Combined_reads$Events_GRL_splitbyread <- c( Combined_reads$Events_GRL_splitbyread, tempdat$Events_GRL_splitbyread ) 
     
     # cleanup memory
     MemLog = gc( verbose = FALSE )
     # and move on to the next one in the loop...
     }
  }

# ========================================================================================
# output the "poremodel" data - so we have the model data saved
if ( ! is.null( output_poremodel ))
  { writeout_pore_model( GRLdat_readsplit  = Combined_reads$Events_GRL_splitbyread,
                         fout              = output_poremodel, 
                         strand_type       = "RNA" ) }

# ========================================================================================
# remove non-finite entries:

table_dat = table_dat [ which (! is.na (table_dat$event_level_mean) ), ]

# =======================
# Add strand information:

table_dat = assign_strand( table_dat, perform_sanity_checks = TRUE )

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
saveRDS( list( "samplename"             = samplename,
               "Flattened_reads"        = Flatten_reads,
               "Events_GRL_splitbyread" = GRL_out ),
         file = output_reads_GRL  )

writeLines("script tables2GR complete: read list and model catalogue printed to output files. Now exiting normally.", logFile )
