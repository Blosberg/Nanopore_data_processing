# ===============================================================
# take a read (as a GR object)  for a single read and eliminate redundance (average overlapping positions)
extract_sequence_from_npread  <- function( read_GR_in     = stop("read_GR_in must be provided")
                                        )
{
  read_input_GR_called <- read_GR_in[ read_GR_in$reference_kmer != "NNNNN" ]
  # -------------------------------------------------------
  read_input_GRL_splitby_startpos <- split( read_input_GR_called, start( read_input_GR_called) )
  redundancy_removed              <- unlist( GRangesList( lapply( read_input_GRL_splitby_startpos, function(x) x[1])))

  kmer_set =  redundancy_removed[ order(redundancy_removed$event_index) ]$model_kmer

  result <- paste( substr( kmer_set, 1,1), collapse = '')
  return(result)
}
# ===============================================================
# collect current signal over a region from _all_ reads that overlap with it.
get_flat_signal_over_all_reads  <- function( reads_GRL_in     = stop("reads_GRL_in must be provided"),
                                             target_range_GR  = stop("target_range must be provided")
                                        )
{
  olaps             <-  findOverlaps( reads_GRL_in, target_range_GR )
  reads_of_interest <-  reads_GRL_in[ queryHits( olaps) ]

  # to "flatten" is to average over multiple events that align to the same position; do this for each read independently
  flattened_reads_of_interest_GRL <- lapply( reads_of_interest, function(x) flatten_read( read_GR_in = x ) )

  # collapse everything back into a single GRanges object
  all_events_all_overlapping_reads_GR <- unlist( GRangesList(  flattened_reads_of_interest_GRL ) )

  # and now split them up by position, ignoring read_index
  all_events_all_overlapping_reads_splitby_startpos_GRL <- split( all_events_all_overlapping_reads_GR,
                                                                  start(all_events_all_overlapping_reads_GR) )
  # -------------------------------------------------------

  chr         = unique( seqnames( all_events_all_overlapping_reads_GR ) )
  str         = unique( strand(   all_events_all_overlapping_reads_GR ) )
  if( length( chr) != 1 || length( str) != 1 )
    { stop("irregular seqnames and strands in read") }

  # Consdier adding another metadata column for # of distinct read_indices hits

  all_reads_flattened  <- unlist( GRangesList( lapply( all_events_all_overlapping_reads_splitby_startpos_GRL,
                                                       function(x) flatten_GR_at_spec_position( readposn_GR_in = x,
                                                                                                chr_in         = chr,
                                                                                                strand_in      = str
                                                                                                )
                                                      )))

  return(all_reads_flattened)
}


# ===============================================================
# take a GR object for a single read and eliminate redundance (average overlapping positions)
flatten_read  <- function( read_GR_in            = stop("read_GR_in must be provided"),
                           perform_sanity_checks = FALSE
                          )
{
  # --- takes a specific read and combines all events associated with the same position
  # --- producing a single value for each position.
  # --- skipped positions are left as "NA"s

  if ( perform_sanity_checks )
    {
    chr         = unique( seqnames( read_GR_in ) )
    str         = unique( strand(   read_GR_in ) )
    if( length( chr) != 1 || length( str) != 1 )
      { stop("irregular seqnames and strands in read") }
    read_index      = unique( read_GR_in$read_index )
    if( length( read_index ) != 1 )
      { stop("in flatten_read: read_GR_in is not associated with a single read") }
  } else {
    chr         = as.character( seqnames( read_GR_in[1] ) )
    str         = as.character( strand(   read_GR_in[1] ) )
    read_index  = read_GR_in[1]$read_index
    }

  read_GR_in_Ns_removed <- read_GR_in[ read_GR_in$model_kmer != "NNNNN" ]
  # ----------

  # create a list, with each element containing groups of GR objects with the same start position
  read_GRL_in_splitby_startpos <- split(   read_GR_in_Ns_removed,
                                           start( read_GR_in_Ns_removed ) )

    if( perform_sanity_checks )
    {
    # If all is ok, it should return 0 everywhere.
    sanity  <- unlist( lapply( read_GRL_in_splitby_startpos,
                       function(x) sanity_check_spec_position ( readposn_GR_in = x)
                              )
                       ) # length( which( sanity != 0 ) )
    }


    # for each starting position, compact the events into a single observation.
    return( unlist (GRangesList(  lapply( read_GRL_in_splitby_startpos,
                                                       function(x) flatten_GR_at_spec_position( readposn_GR_in = x,
                                                                                                chr_in         = chr,
                                                                                                strand_in      = str
                                                                                                )
                                )
                      ) ) )

    # Above condenses individual commands provided here:
    # read_GRL_in_splitby_startpos_flattened  <- lapply( ... )
    # result <- unlist ( GRangesList( read_GRL_in_splitby_startpos_flattened ) )
    # return(  result )

}

# ===============================================================
# make sure the specific position parameters make sense.
sanity_check_spec_position  <- function( readposn_GR_in        = stop("read_GR_in must be provided")
                                       )
{
    S               = unique( start(    readposn_GR_in ) )
    E               = unique( end(      readposn_GR_in ) )
    reference_kmer  = unique( readposn_GR_in$reference_kmer )
    model_kmer      = unique( readposn_GR_in$model_kmer     )
    model_mean      = unique( readposn_GR_in$model_mean     )
    model_stdv      = unique( readposn_GR_in$model_stdv     )

      if ( max ( length(S),
                 length(E),
                 length(reference_kmer),
                 length(model_kmer),
                 length(model_mean),
                 length(model_stdv)   ) > 1
          )
    {
      stop("non-unique values for position-constant parameters")
      }
    else{
      return( 0 )
      }
}
# ===============================================================
# Define, precisely, the one-element GRanges object for each position in the above function.
flatten_GR_at_spec_position  <- function( readposn_GR_in        = stop("read_GR_in must be provided"),
                                          chr_in                = stop("chr_in must be provided"),
                                          strand_in             = stop("strand_in must be provided"),
                                          average_currentdat    = TRUE,
                                          perform_sanity_checks = FALSE
                                          )
{
if( length( readposn_GR_in ) == 1 )
  { return( readposn_GR_in ) }

    S               = start(  readposn_GR_in[1] )
    E               = end (   readposn_GR_in[1] )

    read_index       = readposn_GR_in[1]$read_index
    event_index      = paste0( as.character( min ( readposn_GR_in$event_index ) ),
                               "-",
                               as.character(max ( readposn_GR_in$event_index ) ) )

    model_kmer       = readposn_GR_in[1]$model_kmer

#    reference_kmer  = readposn_GR_in[1]$reference_kmer
#    model_mean      = readposn_GR_in[1]$model_mean
#    model_stdv      = readposn_GR_in[1]$model_stdv
   # ^ redundant metadata columns omitted for efficiency.

  # duration of overlapping events:
  event_length = sum(  readposn_GR_in$event_length)

  event_mean   = sum(  readposn_GR_in$event_mean * readposn_GR_in$event_length ) /event_length
  event_stdv   = sqrt( sum( readposn_GR_in$event_stdv*readposn_GR_in$event_stdv * readposn_GR_in$event_length ) / event_length )


  result <- GRanges( seqnames        = chr_in,
                     ranges          = IRanges ( start  = S,
                                                 end    = E ),
                     strand          = strand_in,

		     read_index      = read_index,
                     event_index     = event_index,

 #                    reference_kmer  = reference_kmer,
 #                    model_mean      = model_mean,
 #                    model_stdv      = model_stdv,
                     # ^ As noted above: omitted for efficiency

                     event_mean      = event_mean,
                     event_stdv      = event_stdv,
                     event_length    = event_length,
                     model_kmer      = model_kmer
                    )

  return(result)
}
