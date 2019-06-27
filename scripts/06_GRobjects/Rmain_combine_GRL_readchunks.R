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

      Arguments:
      Rfuncs_tsv2GRconv       --script with function definitions used here.
      GRL_reads_combined_out  -- combined (across all chunks) rds filename read-separated GRL output.
      output_poremodel        -- rds filename with model date (for each model kmer) saved.
      samplename              -- name of sample for documentation
      logFile                 -- filename to pipe output to
      poremodel_chunks        -- Array of unique kmer models from each chunk.
      GRL_chunk_files         -- array of files to use as input

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

Rfuncs_tsv2GRconv       <- argsL$Rfuncs_tsv2GRconv
GRL_reads_combined_out  <- argsL$GRL_reads_combined_out
# output_poremodel        <- argsL$output_poremodel
sampleName              <- argsL$sampleName
logFile                 <- argsL$logFile
GRL_chunk_files         <- unlist( strsplit(argsL$GRL_chunk_files,",")  )

#===============================================================
suppressPackageStartupMessages( library(GenomicRanges) )
source(Rfuncs_tsv2GRconv)

# number of individual readGRL objects to be merged:
Nchunks         = length( GRL_chunk_files  )

# what will ultimately be the merger of all of them:
Combined_reads  = readRDS( GRL_chunk_files[1] )

if ( Nchunks >= 2 )
  {
  for ( i  in c(2:Nchunks) )
     {
     # count up the total current number of reads:
     Nreads  = length( Combined_reads$Events_GRL_splitbyread )

     # Read in the next chunk of data:
     tempdat = readRDS( GRL_chunk_files[i] )

     # sanity-check that the names are in order:
     # TODO: comment this line out for efficiency once we're convinced this works.
     if( ! identical( as.character( unlist( lapply( tempdat$Events_GRL_splitbyread, function(x) unique( x$read_index ) ) ) ), as.character( c(1:length( tempdat$Events_GRL_splitbyread ) )  )  ) )

       { stop(paste("Disordered names in reads of file: ", GRL_chunk_files[i] ) ) }

     # offset the read index for each one of these:
     tempdat$Events_GRL_splitbyread  = lapply( tempdat$Events_GRL_splitbyread, function(x) offset_read_indices ( GRLchunk_in  = x, Nread_offset = Nreads ) )

     # likewise, offset the corresponding read "names"
     names(tempdat$Events_GRL_splitbyread) <- as.character( lapply( tempdat$Events_GRL_splitbyread, function(x) unique( x$read_index ) ) )

     # Do not offset name. We allow for chimeric alignments
     # As such, there may be multiple possible alignments per read:
     # names( tempdat$Events_GRL_splitbyread )  <-  as.character( as.numeric( names(tempdat$Events_GRL_splitbyread) ) + Nreads )

     # Now, absorb these reads into the "combined" data set
     Combined_reads$Events_GRL_splitbyread <- c( Combined_reads$Events_GRL_splitbyread, tempdat$Events_GRL_splitbyread )

     # cleanup memory
     MemLog = gc( verbose = FALSE )
     # and move on to the next one in the loop...
     }
  }

# TODO: need to consolidate the individual poremodel files.
# ========================================================================================
# output the "poremodel" data - so we have the model data saved
# if ( ! is.null( output_poremodel ))
#   { writeout_pore_model( GRLdat_readsplit  = Combined_reads$Events_GRL_splitbyread,
#                          fout              = output_poremodel,
#                          strand_type       = "RNA" ) }


# output the GRL object to store reads:
saveRDS( Combined_reads,
         file = GRL_reads_combined_out  )

writeLines("script tables2GR complete: read list and model catalogue printed to output files. Now exiting normally.", logFile )
