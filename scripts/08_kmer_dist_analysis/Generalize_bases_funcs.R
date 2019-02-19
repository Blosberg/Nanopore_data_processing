# These are functions for handling the conversion of bases to more generalized constructs:
# e.g., if you want to input ACGNW, and you want the "N" to trace over all of ACGT, and the "W" to trace over the "Weak" bases: AT, etc. See IUPAC naming convention for a comprehensive list.

# ==================================================================
# --- produce list of sequences that trace over a given unknown base
produce_list_unknown_trace <- function( sequence_in  = stop("sequence must be provided"),
                                        base_set     = c("A","C","G","T"),
                                        unknown_mark = "N"
)
{
  if ( nchar(unknown_mark) != 1 )
    stop("Unknown_mark must be a single character")

  N_unknown = str_count(sequence_in, unknown_mark) # how many bases need to be traced over.

  # k            = nchar(sequence_in) # number of bases under consideration.
  # locs_unknown = str_locate_all(teststring, "unknown_mark")
  # result       = list()

  if( N_unknown < 0)
  { stop("ERROR: N_unknown <0")}else if(N_unknown == 0){
    return(sequence_in) } else{

      seqset_1 = list()
      seqset_2 = list()

      seqset_2[[1]] = sequence_in

      for (loop_iteration in c(1:N_unknown))
      {
        Nseqs  = 1
        seqset_1 = list()

        for ( seq_i in seqset_2 )  # Loop through existing sequences (with (N_unknown-loop_iteration) N's)
        {
          for ( base in base_set ) # loop through base being substituted in: A,C,G,T
          {
            temp_sub          <- sub( unknown_mark,
                                      base,
                                      seq_i) # Changes only the 1st pattern match per string
            seqset_1[[Nseqs]] <- temp_sub
            Nseqs             <- Nseqs +1;
          }
        }

        seqset_2  <- seqset_1

      }

      remaining_unknowns = FALSE

      for ( i in c(1:length(seqset_2)))
      {

        if( str_count( seqset_2[[i]], unknown_mark) > 0  )
        {
          remaining_unknowns = TRUE
          stop("unknowns remaining at termination of produce_list_unknown_trace ")
        }
      }

      if( length( unique(unlist(seqset_2 ))) !=  length( unlist(seqset_2 ))   || length( unlist(seqset_2 )) != length(base_set)^N_unknown  )
      { stop("Output of produce_list_unknown_trace either has non-unique values, or unexpected length. Debug!") }

      result = unlist(seqset_2)
      return(result)

    }

}

# ==================================================================
# --- take average of bases in a given position
trace_over_unknown_base <- function( histlist_in  = stop("reads must be provided"),
                                     sequence_in  = stop("sequence must be provided"),
                                     breaks_in    = seq(50, 150, 0.5),
                                     base_set     = c("A","C","G","T"),
                                     unknown_mark = "N"
)
{
  k             = nchar(sequence_in) # number of bases under consideration.

  sequence_list = produce_list_unknown_trace( sequence_in = sequence_in,
                                              base_set    = base_set  )
  # No longer generate the histograms here --- use the ones provided in histlist_in
  # histlist <- lapply(sequence_list,  function(x)   hist( GRL_input[[x]]$event_mean,
  #                                                        breaks=breaks_in,
  #                                                        plot = FALSE) )

  cumulative_hist=matrix( 0, (length(breaks_in)-1), 1)
  for( seq in sequence_list )
    {
    cumulative_hist = cumulative_hist + histlist_in[[seq]]
    }

  return( (1/sum(cumulative_hist))*cumulative_hist )
}

# ==================================================================
# --- plot_traces_over_bases
plot_background_hist <- function( xdat     = stop("xdat must be provided"),
                                  hist_in  = stop("histogram must be provided"))
{
  plot( c(min(xdat), xdat, max(xdat)),
        c(0, hist_in, 0),
        col   = 'grey43',
        lty   = "blank",
        type  = "l",
        xlab  = "Current [pA]",
        ylab  = "freq",
        main  = paste("Current distribution for Homopolymers")
  )
  polygon(c(min(xdat), xdat, max(xdat)),
          c(0, hist_in, 0),
          lty  ="blank",
          col  = adjustcolor('grey43',alpha.f = 0.3) )
}


