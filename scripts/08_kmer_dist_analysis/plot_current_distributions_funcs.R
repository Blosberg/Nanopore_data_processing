# Various functions for plotting current distributions (with background and transparency features, etc.)

# ==================================================================
# --- add colored areas for specific bases with various unknowns
plot_specific_sequence_currenthist <- function( xdat_in     = stop("xdat must be provided"),
                                                histlist_in = stop("histogram list must be provided"),
                                                seq_in      = stop("sequence must be provided"),
                                                scale_HP    = 0.2,
                                                col_in      = rgb( 0,      0,     1,   alpha=0.5 ),
                                                add         = TRUE
                                                )
{
  if( add )
  {
    polygon( c( min(xdat_in),
                xdat_in,
                max(xdat_in)
    ),
    scale_HP * c( 0,
                  trace_over_unknown_base ( histlist_in  = histlist_in,
                                            sequence_in  = seq_in
                  ), # N.B. this need not _necessarily_ contain unknown bases (i.e. "N"'s -- but it CAN handle them.)
                  0),
    col  = col_in,
    lty  = "blank"
    )
  } else {
    plot( c( min(xdat_in),
             xdat_in,
             max(xdat_in)
    ),
    scale_HP * c( 0,
                  trace_over_unknown_base ( histlist_in  = histlist_in,
                                            sequence_in  = seq_in ), # N.B. this need not _necessarily_ contain unknown bases (i.e. "N"'s -- but it CAN handle them.)
                  0),
    col   = col_in,
    lty   = "blank",
    type  = "l",
    xlab  = "Current [pA]",
    ylab  = "freq",
    main  = ""
    )
    polygon( c( min(xdat_in),
                xdat_in,
                max(xdat_in)
               ),
              scale_HP * c( 0,
                            trace_over_unknown_base ( histlist_in  = histlist_in,
                                                      sequence_in  = seq_in ), # N.B. this need not _necessarily_ contain unknown bases (i.e. "N"'s -- but it CAN handle them.)
                              0),
              lty  ="blank",
              col  = col_in
            )
    }

  legend( "topright", 
          legend = seq_in,
          fill   = col_in, 
          lty    = "blank", 
          cex    = 0.8
  )
  
  
}

# ==================================================================
# --- add colored areas for specific bases with various unknowns
plot_homopolymers <- function( xdat_in     = stop("xdat must be provided"),
                               histlist_in = stop("histogram list must be provided"),
                               pos_N       = NULL,
                               scale_HP    = 0.2,
                               col_A = rgb( 0,      0,     1,   alpha=0.5 ),
                               col_C = rgb( 0.3548, 0.4484, 0.1967742,   alpha=0.5 ),
                               col_G = rgb( 1,         0.5,         0,   alpha=0.75 ),
                               col_T = rgb( 1,           0,         0,   alpha=0.5 ),
                               k     = 5,
                               isRNA = TRUE)
{

  seq_A = paste(replicate(k, "A"), collapse = "");
  seq_C = paste(replicate(k, "C"), collapse = "");
  seq_G = paste(replicate(k, "G"), collapse = "");
  seq_T = paste(replicate(k, "T"), collapse = "");

  if( length(pos_N) > 0)
  {
  for(i in c(1:length(pos_N) ) )
    {
    if( i < 1 || i > k )
      { stop( paste("ERROR: trying to substitute base ", as.character(i)," into string of length ", as.character(k))) }

    substr(seq_A, pos_N[i], pos_N[i]) <- unknown_mark;
    substr(seq_C, pos_N[i], pos_N[i]) <- unknown_mark;
    substr(seq_G, pos_N[i], pos_N[i]) <- unknown_mark;
    substr(seq_T, pos_N[i], pos_N[i]) <- unknown_mark;
    }
  }

  seqset <- c( seq_A, seq_C, seq_G, seq_T)
  col_array = c(col_A,
                col_C,
                col_G,
                col_T )
  names(col_array ) <-  seqset

  for ( seq in seqset )
    {
    polygon( c( min(xdat_in),
                xdat_in,
                max(xdat_in)
               ),
              scale_HP * c(0,
                        trace_over_unknown_base ( histlist_in  = current_histlist_by_kmer,
                                                  sequence_in  = seq
                                                  ),
                           0),
            col  = col_array[seq],
            lty  = "blank"
            )
  }
  if( isRNA )
    {
    displayseq = gsub("T","U", seqset)
    } else{
    displayseq = seqset
    }

  legend( "topright",
          legend = c( displayseq, "NNNNN"),
          fill   = c( col_array,
                      adjustcolor('grey43',alpha.f = 0.3) ),
          lty    = "blank",
          cex    = 0.8)
}




