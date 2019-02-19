
# readsdat_splitby_kmer = split( readsdat_in, readsdat_in$model_kmer )

get_synthdat <- function( N                     = stop("Num reads must be provided"),
                          seqs                  = stop("sequence must be provided"),
                          events_splitby_kmer   = stop("unsorted events must be provided")
 )
{
 
list_of_mean_arrays     <- lapply( seqs, function(x)  get_Nvector_for_given_seq( N       = N, 
                                                                                 kmer_GR = events_splitby_kmer[[x]] ) )   

synthdat_result    <- t( simplify2array( list_of_mean_arrays ) )

row.names( synthdat_result ) <- seqs

return( synthdat_result  )

}

#=============
get_Nvector_for_given_seq <-function(  N       = stop("Need size N"),
                                       kmer_GR = stop("Need seqin"))
{
  GR_splitby_kmer_and_read = split( kmer_GR, kmer_GR$read_index )[1:N]
  
  mean_array  <- lapply( GR_splitby_kmer_and_read, function(x)  mean( x$event_mean )  ) 
  # sdv_array   <- lapply( GR_splitby_kmer_and_read, function(x)  mean( x$event_stdv )  ) 
  
  return( unlist(mean_array) )
  
  }