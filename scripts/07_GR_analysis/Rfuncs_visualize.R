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