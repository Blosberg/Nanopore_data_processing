# ===============================================================
# take a read (as a GR object)  for a single read and eliminate redundance (average overlapping positions)
extract_sequence_from_npread  <- function( read_GR_in     = stop("read_GR_in must be provided")
                                        )
{
  read_input_GR_called <- read_GR_in[ read_GR_in$model_kmer != "NNNNN" ]
  # -------------------------------------------------------
  # read_input_GRL_splitby_startpos <- split( read_input_GR_called, start( read_input_GR_called) )

  if( unique( as.character( strand( read_input_GR_called ) ) ) == "+"  ) 
    { 
    
    read_ordered <- sort( read_input_GR_called, decreasing = FALSE )
    
    } else if( unique( as.character( strand( read_input_GR_called ) ) ) == "-"  ){

    read_ordered <- sort( read_input_GR_called, decreasing = TRUE )
    
  }  else{
    stop("ERROR extracting sequence; irregularity in strand information.")
  }

  kmer_set =  read_ordered$model_kmer

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
# --- return normalized histogram for the log-dwell time
get_normalized_lt_hist <- function( GR_kmer_in      = stop("kmer-split event GR obj must be provided"),
                                    lt_histmin = -7,
                                    lt_histmax = 2,
                                    lt_histres = 0.1,
                                    eps = 0.0000001
                                    )
{
 breakset = seq( lt_histmin, lt_histmax, lt_histres)

 dat = log( GR_kmer_in$event_length)
 
 if( min(dat) < lt_histmin || max(dat) > lt_histmax )
   {
   print("WARNING: data lies partially outside histogram range")
   dat[ dat < lt_histmin ] = lt_histmin+eps;
   dat[ dat > lt_histmax ] = lt_histmax-eps;
   }
 
 result <- hist( log( GR_kmer_in$event_length),
                 breaks  = breakset,
                 plot    = FALSE
                )
 return(result)
}
   
