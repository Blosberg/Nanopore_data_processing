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

# ==================================================================
# --- very general (probably obsolete compared to subsequent functions)
# --- returns normalized histogram for a given set of breaks and GR
get_normalized_current_hist <- function( GR_input   = stop("reads must be provided"),
                                         breaks_in  = seq(50,150,0.5)
                                         )
{
  temp <- hist( GR_input$event_mean,
                breaks=breaks_in,
                plot = FALSE)
  result <- (1/sum( temp$counts )) * temp$counts
  return(result)
}

# ==================================================================
# --- select specific sequence from two histlists
seq_spec_compare <-  function( seq        = stop("seq  must be provided"),
                               criterion  = "undefined",
                               SOI_GR     = stop("sites of interest must be provided in GRanges format"),
                               control_GR = stop("control sites must be provided in GRanges format")
                               )
{

SOI_seq_GR     = SOI_GR[     SOI_GR$reference_kmer     == seq ]  # --- sites of interest with corresponding sequence.
control_seq_GR = control_GR[ control_GR$reference_kmer == seq ]  # --- corresponding sequence in control locations (i.e. everywhere else.)

result = list( "SOI_seq_GR" =  SOI_seq_GR,
               "control_seq_GR" = control_seq_GR,
               "seq"=seq,
               "part_criterion" = criterion )

return( result )

}

# ==================================================================
# --- plot the data in the list obtained from the previous function ^
plot_seq_spec_comparison  <- function( seq_spec_list = stop("seq_list  must be provided"),
                                       pore_model    = NULL,
                                       mincurrent    = 50,
                                       maxcurrent    = 150,
                                       res           = 1,
                                       scale         = FALSE )
{

  temp1 = seq_spec_list$SOI_seq_GR[  seq_spec_list$SOI_seq_GR$event_mean  > mincurrent ]
  temp1 = temp1[ temp1$event_mean < maxcurrent ]

  temp2 = seq_spec_list$control_seq_GR[  seq_spec_list$control_seq_GR$event_mean  > mincurrent ]
  temp2 = temp2[ temp2$event_mean < maxcurrent ]

  # current_window = c( mincurrent, maxcurrent )
  # par( mfrow = c(2,1) ) # === putting things into the same window now.

  breakset = seq( from = mincurrent,
                to   = maxcurrent,
                by   = res )

  if( scale ) # scale by (difference from model_mean)/model_stddev
  {
    stat_params = unlist( pore_model[ seq_spec_list$seq, ] )
    if (length(stat_params) != 2)
      {
      stop("Failed to obtain statistical parameters from pore_model reference")
      }
    breakset          <- ((breakset-stat_params[1])/stat_params[2] );
    temp1$event_mean  <- ((temp1$event_mean-stat_params[1])/stat_params[2] );
    temp2$event_mean  <- ((temp2$event_mean-stat_params[1])/stat_params[2] );
    }

  # plot control first in blue, so it appears behind the SOI
  hist( temp2$event_mean,
      freq = FALSE,
      lty="blank",
      col=rgb( 0.0, 0.0, 1, 0.5 ),
      main=paste("sequence=",
                 seq_spec_list$seq,
                 "\nputative vs. bulk alignments (",
                 seq_spec_list$part_criterion,
                 ")"),
      breaks = breakset,
      xlab = "current [pA]",
      ylab = "prob" )

  # plot SOI on top in red
  hist( temp1$event_mean,
      freq = FALSE,                   # we want probabilities
      lty="blank",
      col=rgb( 1, 0.0, 0.0, 0.5 ),
      breaks = breakset,
      add=T)
  # --------
 # legend( 120, 0.08,
 #         inset  = 0.2,
 #         legend = c("putative", "control"),
 #         cex=0.8,
 #         col=c( "red", "blue" ),
 #         box.lwd=1
 #        )

 legend( maxcurrent - (0.4*(maxcurrent-mincurrent)), -0.01, # places a legend at the appropriate place c("Health","Defense"), # puts text in the legend
         lty=c(1,1), # gives the legend appropriate symbols (lines)
         legend =c("putative", "control"),
         lwd=c(2.5,2.5),col=c("red","blue")) # g

return(0);
}

# ==================================================================
# --- calculate overlap of two histograms.
calculate_histogram_overlap <- function(  hist1    = stop("hist1 must be provided"),
                                          hist2    = stop("hist2 must be provided")
                                          )
{
  if( !identical( hist1$mids, hist2$mids) )
    stop("histograms do not share common x-axis binning positions.")

  # two profiles alongside each-other now from each histogram.
  profiles = rbind( diff( hist1$breaks) * hist1$density,  diff( hist2$breaks) * hist2$density )
  
  # take the min of each one, and sum the result.
  result = sum( apply(profiles, 2, FUN=min) )

  return(result);
}

# ==================================================================
# --- Plot overlap of two histograms.
plot_histogram_overlap <- function(   hist1    = stop("hist1 must be provided"),
                                      hist2    = stop("hist2 must be provided"),
                                      seq_in      = stop("seq_in must be provided"),
                                      col_in1  = rgb( 1,      0,     0,   alpha=0.5 ),
                                      col_in2  = rgb( 0,      0,     1,   alpha=0.5 )
)
{
  if( !identical( hist1$mids, hist2$mids) )
    stop("histograms do not share common x-axis binning positions.")
  if( is.null(hist1) || is.null(hist2) )
    stop("Null passed as histogram. Stopping. ")
  
  xdat = hist1$mids

  plot( c( min(xdat),
           xdat,
           max(xdat) ),
        c( 0,
           hist1$density, 
           0),
        col = "black",
        type="l",
        xlab = "Current [pA]",
        ylab = "prob.",
        main = paste( seq_in, " : ", round( calculate_histogram_overlap( hist1 = hist1,  
                                                                                  hist2 = hist2), 
                                                     digits = 4)  )
  )
  
  polygon( c( min(xdat),
              xdat,
              max(xdat)
  ),
  c( 0,
     hist1$density,
     0),
  col  = col_in1,
  border ="black",
  lty  = "blank"
  )
  
  lines( c( min(xdat),
            xdat,
            max(xdat) ),
         c( 0,
            hist2$density,
            0),
         col = "black"
  )
  
  polygon( c( min(xdat),
              xdat,
              max(xdat)
  ),
  c( 0,
     hist2$density,
     0),
  col  = col_in2,
  border ="black",
  lty  = "blank"
  )

}


# ==================================================================
# --- Plot just one histogram.
plot_single_histogram <- function(   hist     = stop("hist1 must be provided"),
                                     seq      = NULL,
                                     col_in   = rgb( 1,      0,     0,   alpha=0.5 ),
                                     add      = FALSE
)
{
  
  if ( is.null(seq) )
    {   
    y = hist } else{
    y = hist[[seq]] 
    }
  
  xdat = y$mids
  
  if( add ) 
  {  lines( c( min(xdat),
           xdat,
           max(xdat) ),
        c( 0,
           y$density, 
           0),
        col = "black",
        type="l",
        xlab = "Current [pA]",
        ylab = "prob.",
        main = seq
  )  } else {
  plot( c( min(xdat),
           xdat,
           max(xdat) ),
        c( 0,
           y$density, 
           0),
        col = "black",
        type="l",
        xlab = "Current [pA]",
        ylab = "prob.",
        main = seq
  )
  }
    # ----------------------
  
  polygon( c( min(xdat),
              xdat,
              max(xdat)
  ),
  c( 0,
     y$density,
     0),
  col  = col_in,
  border ="black",
  lty  = "blank"
  )
}
  
# ==================================================================
# --- get mean
hist_mean <- function(   hist_in    = stop("hist must be provided"),
                         seq_in     = stop("seq must be provided")
                         )
{
  return( sum( hist_in$density* hist_in$mids * diff(hist_in$breaks) ) )
  #                P(x)             x              dx 
  }


######@@@@@@@ HAVENT LOOKED BELOW HERE: @@#$######


# Various functions for plotting current distributions (with background and transparency features, etc.)

# ==================================================================
# --- add colored areas for specific bases with various unknowns
plot_specific_sequence_currenthist <- function( 
                                                histlist_in = stop("histogram list must be provided"),
                                                seq_in      = stop("sequence must be provided"),
                                                scale_HP    = 1,
                                                col_in      = rgb( 0,      0,     1,   alpha=0.5 ),
                                                add         = TRUE
                                                )
{
  xdat = histlist_in[[ 1 ]]$mids
  
  if( add )
  {
    polygon( c( min(xdat),
                xdat,
                max(xdat)
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
    plot( c( min(xdat),
             xdat,
             max(xdat)
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
    polygon( c( min(xdat),
                xdat,
                max(xdat)
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




