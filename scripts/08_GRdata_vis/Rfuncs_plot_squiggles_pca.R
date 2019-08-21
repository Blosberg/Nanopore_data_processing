# ===============================================================
# --- move GRanges ob downstream:
shift_GR_streamwise <- function( GRange_ob_in = stop("query must be specified"),
                                 n = 1,
                                 direction = "down"
                                )
{
result =  GRanges( seqnames = seqnames( GRange_ob_in),
                   ranges   = IRanges( start = start( GRange_ob_in ),
                                       end   = end(   GRange_ob_in )),
                   strand   = strand( GRange_ob_in )
                   )

 strand_p_list = which( as.character(strand(GRange_ob_in) ) == "+") 
 strand_n_list = which( as.character(strand(GRange_ob_in) ) == "-") 
 
 if( direction == "up")
  {
   start( result[ strand_p_list ] ) =  start( result[ strand_p_list ] ) - n
   end(   result[ strand_p_list ] ) =  end(   result[ strand_p_list ] ) - n

   end(   result[ strand_n_list ] ) =  end( result[ strand_n_list ] ) + n
   start( result[ strand_n_list ] ) =  start( result[ strand_n_list ] ) + n
   } else if( direction == "down") {
     
   end(   result[ strand_p_list ] ) =  end(   result[ strand_p_list ] ) + n
   start( result[ strand_p_list ] ) =  start( result[ strand_p_list ] ) + n

   start( result[ strand_n_list ] ) =  start( result[ strand_n_list ] ) - n
   end(   result[ strand_n_list ] ) =  end(   result[ strand_n_list ] ) - n
   } else{ stop("invalid direction") }
 
return( result)
}
# ===============================================================
# --- get the reference sequence for a given GRob:
get_sequence_at_reference <- function( GRange_ob_in = stop("query must be specified"),
                                       refGen_in    = stop("reference must be specified"),
                                       upstream_incl      = 2,
                                       dnstream_incl      = 2
                                       )
{
if( as.character( strand( GRange_ob_in ) ) == "+" )
  {
  leftseg  = upstream_incl;
  rightseg = dnstream_incl;
  } else if(as.character( strand( GRange_ob_in ) ) == "-") {
  rightseg = upstream_incl;
  leftseg  = dnstream_incl;
  } else{
    stop("invalid strand type.")
    }
  
  
Biostrings_seq = eval( parse( text= paste0( "refGen_in$",
                                            seqnames(GRange_ob_in), 
                                            "[", as.character( start(GRange_ob_in) - leftseg  ),
                                            ":", as.character( end(GRange_ob_in)   + rightseg ), "]"
                                           )
                            )
                      )
if( as.character( strand( GRange_ob_in ) ) == "-" )
  { Biostrings_seq <- reverseComplement( Biostrings_seq ) }

return( as.character( Biostrings_seq ))
}
# ===============================================================
# --- plot the current lines from all the reads of a sample over a given region:
plot_hist_with_edge_trimming <- function ( data      = stop("data must be provided"),
                                           breaks_in = seq( 50, 150, 1),
                                           col       = rgb( 0, 0, 1, 0.5 ),
                                           add       = F,
                                           eps       = 0.0001,
                                           main      = "Untitled"
                                           )
{
  
  data [ data < min( breaks_in ) ] <- min( breaks_in ) + eps
  data [ data > max( breaks_in ) ] <- max( breaks_in ) - eps
  breakset = breaks_in 
  
  hist( data, 
        freq = FALSE,
        breaks = breakset,
        lty="blank",
        col=col, 
        add = add,
        main = main
  )
}
# ===============================================================
# --- plot the current lines from all the reads of a sample over a given region:
plot_samplesignal_over_ROI  <- function( SampleName       = "unnamed",
                                         sampleROI_dat    = stop("aligned sampledat must be provided"),
                                         refgen           = stop("refgen must be provided"),
                                         mincurrent       = 50,
                                         maxcurrent       = 150,
                                         normed           = FALSE,
                                         squiggle_type    = "none",
                                         shading_darkness = 0.3,
                                         mean_darkness    = 0.5,
                                         line_darkness    = 0.1
                                         )
{
  Nxpos  = length( sampleROI_dat$xposn_vals)
  Nreads = dim( sampleROI_dat$posn_means )[1]
  chr    = as.character( unique( seqnames( sampleROI_dat$ROI )))

  normdev_range = 4

  # =====================================================================
  if( normed ){ # plotting deviations from normal (i.e. (signal - mean)/stddev)
    # Plot a dashed line for the mean (==0). Everything after that should be "lines"

      if( squiggle_type =="full" ){ stop("TODO: Haven't handled this case yet. Code this.") }

    plot(  sampleROI_dat$xposn_vals,
           rep( 0, length( sampleROI_dat$xposn_vals ) ),
           col=rgb(0, 0, 0, mean_darkness ),
           type="l",
           lty = "dashed",
           lwd=2,
           xlim = c( min(sampleROI_dat$xposn_vals),
                     max(sampleROI_dat$xposn_vals) ),

           main= paste(  SampleName,
                         "; normalized \n", chr, ":",
                         as.character( min(sampleROI_dat$xposn_vals) ), "-",
                         as.character( max(sampleROI_dat$xposn_vals)  ) ),
           xlab = "reference sequence",
           xaxt="n",
           ylab = "current [pA]",
           ylim = c( (-1*normdev_range), normdev_range )
      )

     # plot a transparent polygon around the +/-1 std.dev region of the standard (expected) current values.
     polygon( c( sampleROI_dat$xposn_vals, rev(sampleROI_dat$xposn_vals) ),
              c( rep(-1, Nxpos), rep(1,Nxpos)),
              col  = rgb(0, 0, 0, shading_darkness ),
              lty  = "blank"
            )
  #------  NOW PLOT THE READS THEMSELVES ------
  # positions in each read without current data are left "NA"

  for (i in c(1:Nreads) )
    {
    lines( sampleROI_dat$xposn_vals,
           sampleROI_dat$read_normdiff[i, ],
           col    = rgb(0, 0, 1, line_darkness),
           type   = "l",
           lwd    = 1
          )
  }

  # add reference sequence, complement, and directionality.
  add_sequence_tick_marks( refgen   = refgen,
                           Locus    = sampleROI_dat$ROI,
                           min_ypos = (-1*normdev_range)  )

  # =====================================================================
  }else { # i.e. not the deviation from normal.

  poremodel_Icurves = matrix( 0,
                              Nxpos,
                              3 )
  # now fill in the min/mean/max values into three columns of the output respectively.
  poremodel_Icurves[,1] = sampleROI_dat$poremodel_metrics[,1] - sampleROI_dat$poremodel_metrics[,2]
  poremodel_Icurves[,2] = sampleROI_dat$poremodel_metrics[,1]
  poremodel_Icurves[,3] = sampleROI_dat$poremodel_metrics[,1] + sampleROI_dat$poremodel_metrics[,2]


  if( squiggle_type =="full" || squiggle_type == "event" )
  { xmean_dashed = c( rep(sampleROI_dat$xposn_vals,each=2)[2:(2*Nxpos) ],  sampleROI_dat$xposn_vals[Nxpos]+1 )
    ymean_dashed = rep(poremodel_Icurves[ , 2],each=2)

    shaded_boundary_x = c( xmean_dashed, rev(xmean_dashed) )
    shaded_boundary_y = c( rep( (poremodel_Icurves[ , 1]),each=2),
                             rev( rep( poremodel_Icurves[ , 3],each=2) ) )

  } else {
     xmean_dashed = sampleROI_dat$xposn_vals;
     ymean_dashed = poremodel_Icurves[ , 2];

     shaded_boundary_x = c( sampleROI_dat$xposn_vals, rev(sampleROI_dat$xposn_vals) )
     shaded_boundary_y = c( (poremodel_Icurves[ , 1]),
                             rev(poremodel_Icurves[ , 3]) )

    }
  # Plot a dashed line for the mean. Everything after that should be "lines"
  plot(  xmean_dashed,
         ymean_dashed,
         col=rgb(0, 0, 0, mean_darkness),
         type="l",
         lty = "dashed",
         lwd=2,
         xlim = c( min(sampleROI_dat$xposn_vals),
                   max(sampleROI_dat$xposn_vals) ),

         main= paste(  SampleName, "\n",
                       chr,
                       ":",
                       as.character( min(sampleROI_dat$xposn_vals) ),
                       "-",
                       as.character( max(sampleROI_dat$xposn_vals)  ) ),
         xlab = "reference sequence",
         xaxt="n",
         ylab = "current [pA]",
         ylim = c( mincurrent,
                   maxcurrent )
      )

  # plot a transparent polygon around the +/-1 std.dev region of the standard (expected) current values.
  polygon( shaded_boundary_x,
           shaded_boundary_y,
           col  = rgb(0, 0, 0, 0.5*shading_darkness),
           lty  = "blank"
          )

  #------  NOW PLOT THE READS THEMSELVES ------
  # positions in each read without current data are left "NA"
  if( squiggle_type == "full" || squiggle_type=="event" ){
      for (read_i in c(1:Nreads) )
      {
      squiggle_dat = get_single_squiggle( sampleROI_dat = sampleROI_dat,
                                          squiggle_type = squiggle_type,
                                          read_i        = read_i )
      color_set = rgb(0, 0, 1, line_darkness) 
      
      lines( squiggle_dat$x,
             squiggle_dat$y,
             col    = color_set,
             type   = "l",
             lwd    = 1
            )
      }

  } else{
    for (read_i in c(1:Nreads) )
      {
      
      # if( *Set condition as needed* )
      # { color_set = rgb( 1,0,0, line_darkness)} else{ 
      color_set = rgb(0, 0, 1, line_darkness) 
      # }
      
      lines( sampleROI_dat$xposn_vals,
             sampleROI_dat$posn_means[read_i, ],
             col    = color_set,
             type   = "l",
             lwd    = 1
            )
    }
      polygon( shaded_boundary_x,
           shaded_boundary_y,
           col  = rgb(0, 0, 0, 0.5*shading_darkness),
           lty  = "blank"
          )
       lines( xmean_dashed,
              ymean_dashed,
              col=rgb(0, 0, 0, mean_darkness),
              type="l",
              lty = "dashed",
              lwd=2 )

  }
  # add reference sequence, complement, and directionality.
  add_sequence_tick_marks( refgen   = refgen,
                           Locus    = sampleROI_dat$ROI,
                           min_ypos = mincurrent )
  }# done the "if" as to whether the plot should be normed diff or not.

}
# ===============================================================
# --- plot the passage time of a sample over a given region:
plot_dwelltime_over_ROI <-  function( SampleName     = "Unnamed",
                                      sampleROI_dat  = stop("aligned sampledat must be provided"),
                                      refgen         = stop("Reference genome must be provided"),
                                      plot_logdwell  = FALSE )
{
  Nxpos  = length( sampleROI_dat$xposn_vals)
  Nreads = dim( sampleROI_dat$Event_duration_chars )[1]
  chr    = as.character( unique( seqnames( sampleROI_dat$ROI )))

  dwell_times = matrix( NA, Nreads, Nxpos )
  for ( i in c(1:Nreads))
    {
    for ( j in c(1:Nxpos ) )
      {
      dwell_times[i,j] = sum( extract_numericArray_from_charlist( sampleROI_dat$Event_duration_chars[i,j] ) )
      }
    }

  if( plot_logdwell )
    { dwell_times[!is.na(dwell_times)] <- log( dwell_times[!is.na(dwell_times)]  ) }

  min_T          = min( dwell_times, na.rm = TRUE )
  max_T          = max( dwell_times, na.rm = TRUE )


  plot(  sampleROI_dat$xposn_vals,
         dwell_times[1, ],
         col=rgb(0, 0, 0, 0.5),
         type="l",
         lwd=1,
         xlim = c( start(sampleROI_dat$ROI),
                   end(sampleROI_dat$ROI) ),

         main= paste(  SampleName, "\n",
                       chr, ":",
                       as.character( start(sampleROI_dat$ROI) ), "-",
                       as.character( end(sampleROI_dat$ROI)  ) ),
         xlab = "reference sequence",
         xaxt = "n",
         ylab = "dwell time [s]",
         ylim = c( min_T, max_T)
      )

  # add reference sequence, complement, and directionality.
  add_sequence_tick_marks( refgen   = refgen,
                           Locus    = sampleROI_dat$ROI,
                           min_ypos = min_T  )

    # and plot the mean current vals over this range:
    # positions in each read without current data are left "NA"
  # plot the rows of this matrix for each read
  for (i in c(2:Nreads) )
    {
    lines( sampleROI_dat$xposn_vals,
           dwell_times[i, ],
           col    = rgb(0,0,1,0.5),
           type   = "l",
           lwd    = 1
         )
  }

}
# ===============================================================
# --- Add sequence markers to the horizontal axis.
add_sequence_tick_marks   <- function(  refgen        = stop("Reference genome must be defined"),
                                        Locus         = stop("Locus must be defined"),
                                        Col_sense     = rgb( 0, 0, 0, 1  ),
                                        Col_antisense = rgb( 0, 0, 0, 0.5),
                                        min_ypos      = 50,
                                        max_ypos      = 150,
                                        strandtype    = "RNA"
                                        )
{
  # chromosome coordinates, and sequence values:
  xvals             = c( start( Locus):end(Locus) )
  chrom_seq         = eval( parse( text=paste0( "refgen$", as.character( unique( seqnames(Locus))) ) ) )
  seqchars_plus     = chrom_seq[ xvals ]
  seqchars_neg      = chartr("ACGT", "TGCA", seqchars_plus )

  # Now label give plus/neg strand appropriate emphasis:
    if ( as.character( strand( Locus )) == "+" )
  { x0 = start(  Locus )
    x1 = end(Locus)
    Col_plusstrand = Col_sense
    Col_negstrand  = Col_antisense  } else if( as.character( strand(Locus)) =="-" ){

    x0 = end(   Locus )
    x1 = start( Locus )
    Col_plusstrand = Col_antisense
    Col_negstrand  = Col_sense

    } else{ stop("ill-defined strand.")}

  # Now label the sequence in place:
  for ( xi  in c(1:length(xvals)) )
  { axis(1,
         at = xvals[xi],
         seqchars_plus[ xi ],
         padj = 1,
         col.axis = Col_plusstrand
         )
     axis(1,
         at = xvals[xi],
         seqchars_neg[ xi ],
         padj = -0.5,
         col.axis = Col_negstrand
         )
  }

  # Finally, draw an arrow indicating direction of transcription.
  arrows( x0, min_ypos, x1, min_ypos,
          length = 0.25, angle = 20,
          code = 2, col = par("fg"), lty = par("lty"),
          lwd = 3 )

}

# ===============================================================
# --- get a single-read current line trace:
get_single_squiggle <- function ( sampleROI_dat = stop("sampleROI_dat required"),
                                  squiggle_type = "full", # "full"  = all time sample observatiosn plotted
                                                          # "event" = horiontal lines for each event plotted
                                                          # "none"  = single value for all positions (events flattened).
                                  read_i        = stop("read_index required.")
                                )
{
  xpos_in = sampleROI_dat$xposn_vals

  x_out = c()
  y_out = c()

  # =============================================================
  if( squiggle_type == "full")
  {
  ydat_char_in = sampleROI_dat$Isamples_chars[read_i, ]

    for (xi in c(1:length(xpos_in) ) )
    {
    if (! is.na( ydat_char_in[xi] )) {

      temp_y = extract_numericArray_from_charlist( ydat_char_in[xi]  )
      temp_x = seq( from = xpos_in[xi],
                    to   = xpos_in[xi] +1,
                    by   = 1/length(temp_y) )
      temp_x = temp_x[1:(length(temp_x)-1)]
    }else{
      temp_x = NA
      temp_y = NA
    }

    x_out = c( x_out, temp_x )
    y_out = c( y_out, temp_y )
    }
  # =============================================================
} else if ( squiggle_type == "event" ) {

  ydat_char_in = sampleROI_dat$Event_means_chars[read_i, ]
  ydat_vals = sampleROI_dat$Event_means_chars[read_i, ]
  ydat_durs = sampleROI_dat$Event_duration_chars[read_i, ]

    for (xi in c(1:length(xpos_in) ) )
    {
    if (! is.na( ydat_char_in[xi] )) {

      # convert character strings into floats:
      temp_y   = extract_numericArray_from_charlist( ydat_vals[xi]  )
      temp_y   = rep( temp_y,
                      each=2)

      # cumulative times of each events gives us coordinates:
      temp_t   = cumsum( extract_numericArray_from_charlist( ydat_durs[xi]  ) )

      # divide the duration of events into fractions.
      position_total_duration = max( temp_t )
      temp_x = c( 0, rep( temp_t,
                          each=2) )/position_total_duration

      # create an array from this position to the next representing time.
      temp_x = xpos_in[xi] + temp_x[1: (length(temp_x)-1) ]

    }else{
      temp_x = NA
      temp_y = NA
    }

    x_out = c( x_out, temp_x )
    y_out = c( y_out, temp_y )
    }

}# end "if" for which type of squiggle we want.

  return( list( "x" = x_out,
                "y" = y_out) )
}


# ===============================================================
# --- take numerical array data from a string of ,-separated numbers in character format.
extract_numericArray_from_charlist <- function ( sample_char_in  = stop("sample_char required")
                                     )
{ 
  if ( nchar ( sample_char_in ) == 0 )
  { return( NA ) } else{
  return( as.numeric( unlist(strsplit( sample_char_in, "," ) ) ) )
  }
}
# ===============================================================
# --- filter a GRanges object for overlap with a putative location of intererest
iterate_pca_clusters  <- function( sampleROI_dat_in  = stop("GRanges must be provided"),
                                   k_in              = 5,
                                   method            = "pca", # TODO: add a KL-div based method.
                                   shouldplot        = TRUE,
                                   should_zeropadd   = FALSE,
                                   should_prune_outliers = TRUE,
                                   should_prune_cuts     = TRUE,
                                   maxiter           = 1000,
                                   min_outli_dist    = 3,
                                   min_outli_NN      = 1,
                                   mincov            = 10,
                                   retrial_num       = 10
                                   )
{ Flag = ""

  # length of the locus at hand
  m = width( sampleROI_dat_in$ROI )

  Nreads  = dim( sampleROI_dat_in$read_normdiff )[1]
  Nxpos   = dim( sampleROI_dat_in$read_normdiff )[2]
  
  # --- Sanity checks: ------
  if( Nreads < mincov || Nxpos < 10 )
    { stop("Inadequate input data to iterate") }
  
  margin = ((Nxpos-m)/2)
  if( margin %% 1 > 0  || (Nxpos-2*margin) < m )
  { stop("irregular margin size near ROI") }
  
  if( k_in < 3 )
    { stop("k<3 doesn't make any sense. Check your k value.") }
  test_range =  c( (margin - k_in + 1) : ( margin + m - 1) )

  # ---
  # gather chromosome locations that refer to data that includes the ROI
  measurement_locations = as.character( c( ( start(sampleROI_dat_in$ROI) - k_in ) : (end(sampleROI_dat_in$ROI)-1) )  )
  
  # Then convert the above into indices in the plotrange matrix.
  # this is the range of indices from our matrices that will be included in the calculations
  # because they refer to measurements whose range of effects overlap with the ROI
  kcomp_range = unlist( lapply( c(1:length(measurement_locations)),
                                 function( i ) 
                                   which ( colnames( sampleROI_dat_in$read_normdiff ) == measurement_locations[i] )
                              ) 
                      )
  # Sanity check: ensure the location "names" are correct:
  
  
  if (  ! identical( test_range, kcomp_range ) )
    { 
    stop("invalid kcomp range set obtained.")
    }
  
  
  # get the subset of data points near the centre of the ROI:
  ROI_dat_lin = sampleROI_dat_in$read_normdiff[ , kcomp_range, drop = F ]
  ROI_dat_lin[ is.na( ROI_dat_lin )  ] <- 0
  Nplotpoints = dim(ROI_dat_lin)[2]

  # ----------- plot raw pca -----------------
  if( shouldplot )
    {
     clusterdat  <- get_clusters ( 
                                   ROI_dat_lin     = ROI_dat_lin,
                                   should_zeropadd = should_zeropadd
                                  )
     pca_raw <- get_pca( ROI_dat_lin      = ROI_dat_lin,
                         shouldplot       = shouldplot,
                         cluster          = clusterdat$cluster )
  }
  # ----------- Remove reads that cut within the plot region  i.e. "partial overlaps" -----------
  if( should_prune_cuts )
    {
    olap_cuts <- identify_olap_cuts( plotdat_locus_in   = sampleROI_dat_in,
                                     kcomp_range        = kcomp_range )
    
    if( length( olap_cuts ) > 0.5* Nreads ) {
      # IF more than half of the reads cut out near the ROI, then add a warning to the flag
      Flag = paste( Flag, "SPLICE", sep = "+") }
    
    pruned_dat = ROI_dat_lin[ -c(olap_cuts), , drop = F ]
    
  } else{ pruned_dat = ROI_dat_lin[ -c(olap_cuts), , drop = F]   }
  
  # V -- precation in case nrows = 1. Shouldn't be necessary, but is in R.  
  if(  is.null( dim( pruned_dat ) ) ||  ( dim( pruned_dat )[1] < mincov) )
    {
    Flag = paste( Flag, "BELOWMINCOV", sep = "+")
    return( list( "Flag" = Flag ) )
    } 
  # ----------- Remove outliers ----------------------------------------
  outliers   <- identify_outliers( dat_in         = pruned_dat, 
                                   min_outli_dist = min_outli_dist,
                                   min_outli_NN   = min_outli_NN
                                  )
  if( should_prune_outliers && length( outliers ) >= 1 )
    {
    pruned_dat = pruned_dat[ -c(outliers), , drop = F]
    
    Flag = paste( Flag, "OUTLIERS", sep = "+")
    
  } # else do nothing.
  
  Nreads_pruned = dim( pruned_dat )[1]
  if( Nreads_pruned < mincov )
    {
    Flag = paste( Flag, "BELOWMINCOV", sep = "+")
    return( list( "Flag" = Flag ) )
  }
  
  dist_dat = dist( pruned_dat )
  
    # DEBUGGING: CHECKING THE RESULTS WITHOUT PRUNING OUTLIERS:
    # clusterdat  <- get_clusters ( ROI_dat_lin      = pruned_dat,
    #                             should_zeropadd  = should_zeropadd
    #                            )
    # 
    # pca <- get_pca(   ROI_dat_lin    = pruned_dat,
    #                 shouldplot       = TRUE,
    #                 cluster          = clusterdat$cluster
    #             )
  # ----------- gather pca components clusters -----------------
  pca <- get_pca(   ROI_dat_lin      = pruned_dat,
                    shouldplot       = shouldplot ) # clusterdat is NULL by default.
  
  # ----------- gather initial clusters -----------------
  clusterdat  <- get_clusters ( ROI_dat_lin      = pruned_dat,
                                should_zeropadd  = should_zeropadd
                               )
  Silh_score <- mean( silhouette( x    = clusterdat$cluster,
                                  dist = dist_dat )[, 3 ])
  
  # ----------- see if clustering can be improved: -----------------
  
  retrial_counter = 0
  while( retrial_counter < retrial_num  )
    {

    clusterdat_maybe_improved <- get_clusters ( ROI_dat_lin      = pruned_dat,
                                               should_zeropadd  = should_zeropadd
                                              )
    Silh_maybe_improved <- mean(  silhouette( x    = clusterdat_maybe_improved$cluster,
                                              dist = dist_dat )[,3] )
    if( Silh_maybe_improved > Silh_score )
      {
      # print( paste( "silhouette updated from ", 
      #               as.character( Silh_score ), 
      #               " to ", 
      #               as.character(Silh_maybe_improved) ) 
      #        )
      clusterdat <- clusterdat_maybe_improved
      Silh_score <- Silh_maybe_improved
      retrial_counter <- 0
    } else { retrial_counter = retrial_counter + 1 }
    
  #  plot_pca_clusters( pca       = pca, 
  #                     cluster   = clusterdat_maybe_improved$cluster,
  #                     PlotTitle = "Testing")
  #  Sys.sleep(1)

  }    
  
  return( 
         list( "clusterdat" = clusterdat,
               "Silh_score" = Silh_score,
               "pca" = pca,
               "ROI" = sampleROI_dat_in$ROI,
               "xposn_vals" = sampleROI_dat_in$xposn_vals,
               "poremodel_metrics" = sampleROI_dat_in$poremodel_metrics,
               "Flag" = Flag
               )
          )
}

# ===============================================================
# ---  remove reads that slip off the pore at exactly this position.
identify_olap_cuts   <- function( plotdat_locus_in = stop("plotdat_locus must be provided"),
                                  kcomp_range      = stop("plotrange must be specified."), 
                                  min_frac         = 0.5 )
{
  Nreads     =  dim ( plotdat_locus_in$read_normdiff )[1]
  filter_cut = rep( FALSE, Nreads)
  
  full_normdat <- plotdat_locus_in$read_normdiff
  
  # get the subset of data points near the centre of the ROI:
  ROI_dat_lin = plotdat_locus_in$read_normdiff[ , kcomp_range ]
  Nplotpoints = dim( ROI_dat_lin )[2]
  
  # count how many of the components are measured (i.e. not-NaN ) in each read.
  krange_count <- rowSums( as.matrix( !is.nan( ROI_dat_lin ) ) )
  
  # those that contain fewer components than the min frac should be "cut"
  filter_cut[ which( (krange_count / Nplotpoints) < min_frac ) ] <-  TRUE
  
  # now check that the read continues beyond on either side:
  minpos_finite_reading <- unlist( lapply( c(1:Nreads), 
                                     function(read_i) 
                                        min( which( !is.na( full_normdat[ read_i , ] ) ) ) 
                                          ) # get the lowest non-NaN position value 
                                  ) # i.e. the left-most current observation
    
  maxpos_finite_reading <- unlist( lapply( c(1:Nreads), 
                                      function(read_i) 
                                        max( which( !is.na( full_normdat[ read_i , ] ) ) ) 
                                          ) # get the highest non-NaN position value 
                                  ) # i.e. the right-most current observation
  
  # If the read was cut-off (spliced/terminated/etc) within a nt. of the ROI, then we expect weird shit to happen there,
  # and should ignore deviations from the canonical model (that are unlikely to have nothing to do with modifications.) 
  Read_was_cutoff_mid_krange <- ( minpos_finite_reading >= min(kcomp_range) | maxpos_finite_reading < max(kcomp_range) )
  
  # So remove them.
  filter_cut[ which( Read_was_cutoff_mid_krange ) ] <- TRUE
  
  # return the indices of reads that should be "cut" from clustering consideration 
  return( which(filter_cut ) )
}
# ===============================================================
# --- remove reads that slip off the pore at exactly this position.
get_range  <-function()
{
    # lapply get_range()
    # get max/min of where the current value is not nan.
  

  
}
# ===============================================================
# --- remove outliers from data
identify_outliers  <- function( dat_in         = stop("datin must be provided"),
                                min_outli_dist = 2,
                                min_outli_NN   = 2
                               )
{ # N.B. Dat_in is in units of std. dev's from the standard model. 
  # We define an outlier as any point that is further from ALL other points by at least
  # min_outli_dist number of standard deviations.
  
  dist_mat <- as.matrix( dist( dat_in )) # N.B. this will be a symmetric matrix with 0s on the main diagonal.
  diag( dist_mat ) <- min_outli_dist + 1   # we set this to a dummy value to make calculation below easier.
  Nreads <- dim(dist_mat)[1]

  outlier_status <- unlist ( lapply( c(1:Nreads), 
                                function( read_i ) 
                                     length( which( dist_mat[ read_i, ]  < min_outli_dist ) ) <= min_outli_NN   
                                    #   ^^^ the number of other reads that are within the min_outli_dist of this read
                                    )
                            )
    
  return ( which( outlier_status ) )  
  
  # We don't bother with changing around these values anymore.
  # result = dat_in 
  # result$read_normdiff        <-        dat_in$read_normdiff[ -c(read_indexes_to_delete), ] 
  # result$posn_means           <-           dat_in$posn_means[ -c(read_indexes_to_delete), ]
  # result$Isamples_chars       <-       dat_in$Isamples_chars[ -c(read_indexes_to_delete), ]
  # result$Event_means_chars    <-    dat_in$Event_means_chars[ -c(read_indexes_to_delete), ]
  # result$Event_duration_chars <- dat_in$Event_duration_chars[ -c(read_indexes_to_delete), ]
  # 
  return( result )
}

# ===============================================================
# --- Collect cluster data:
get_clusters  <- function( ROI_dat_lin       = stop("GRanges must be provided"),
                           method            = "pca", # TODO: add a KL-div based method.
                           should_zeropadd   = FALSE,
                           maxiter           = 1000
                          )
{
  Nreads      = dim( ROI_dat_lin )[1]
  Nplotpoints = dim( ROI_dat_lin )[2]
  N0pad  = 10 * Nreads
  
if( should_zeropadd )
    {
    zero_array     = matrix(  0, 
                              N0pad, 
                              dim( ROI_dat_lin )[2] )
    
    padded_normdat = rbind( zero_array, 
                            ROI_dat_lin )
    row.names( padded_normdat )[1: Nreads] <- rep( "dummy_0pad", Nreads )
    
    clust_kmeans_padded <- kmeans( x        = padded_normdat,
                                   centers  = 2,
                                   iter.max = maxiter )
    # clust_kmeans       <- kmeans( x        = ROI_dat_lin,
    #                               centers  = clust_kmeans_padded$centers,
    #                               iter.max = 100 )
    
    padded_cluster <- unique( clust_kmeans_padded$cluster[ c(1:N0pad) ] )
    if( length( padded_cluster )!= 1 )
      { stop("padded cluster points are being split into separate clusters. This should never happen. Something is wrong with get_pca_clusters") }
    
    clust_kmeans=clust_kmeans_padded
    clust_kmeans$cluster <- clust_kmeans$cluster[ (N0pad+1) : (N0pad+Nreads) ] 
    
    # N.B. we don't yet concern ourselves with which one is at the origin.
    clust_kmeans$centers[1, ] <- get_cluster_mean( data_in     = ROI_dat_lin [ which( clust_kmeans$cluster ==1 ) , ],
                                                   Nplotpoints = Nplotpoints )
    clust_kmeans$centers[2, ] <- get_cluster_mean( data_in     = ROI_dat_lin [ which( clust_kmeans$cluster ==2 ) , ],
                                                   Nplotpoints = Nplotpoints)
      
    clust_kmeans$totss <- collect_sum_squared_diff ( data        = ROI_dat_lin,
                                                     Nplotpoints = Nplotpoints )
    clust_kmeans$withinss[1]  <- collect_sum_squared_diff ( data        = ROI_dat_lin[ which( clust_kmeans$cluster ==1 ), ],
                                                            Nplotpoints = Nplotpoints )
    clust_kmeans$withinss[2]  <- collect_sum_squared_diff ( data        = ROI_dat_lin[ which( clust_kmeans$cluster ==2 ), ], 
                                                            Nplotpoints = Nplotpoints )
    
    clust_kmeans$tot.withinss <- sum( clust_kmeans$withinss)
    clust_kmeans$betweenss    <- clust_kmeans$totss - clust_kmeans$tot.withinss
    
    clust_kmeans$size[1] <- length( which ( clust_kmeans$cluster == 1 ))
    clust_kmeans$size[2] <- length( which ( clust_kmeans$cluster == 2 ))
    
  } else{ 
    clust_kmeans <- kmeans( x        = ROI_dat_lin,
                            centers  = 2,
                            iter.max = maxiter  )    
    }
  
  if( min( clust_kmeans$size) > 0 )
  {
    dev = sqrt( rowSums( clust_kmeans$centers * clust_kmeans$centers  )  )
    # the cluster labels are arbitrary. We impose that cluster 1 is the one closer to the origin
    # and cluster 2 is the other one. If this is not satisfied, we reverse the order.

    if ( dev[2] < dev[1] )
    { # reverse the order labelling of the clusters s.t. center[1] < center[2]
    clust_kmeans_copy = clust_kmeans

    clust_kmeans_copy$cluster[ clust_kmeans$cluster == 1 ] =2
    clust_kmeans_copy$cluster[ clust_kmeans$cluster == 2 ] =1

    clust_kmeans_copy$centers  <- clust_kmeans$centers[  c(2,1),]
    clust_kmeans_copy$withinss <- clust_kmeans$withinss[ c(2,1) ]
    clust_kmeans_copy$size     <- clust_kmeans$size[     c(2,1) ]

    clust_kmeans <- clust_kmeans_copy
    rm( clust_kmeans_copy )
    } # else do nothing.
  }
  
return( clust_kmeans) 
  
}
# ===============================================================
# --- filter a GRanges object for overlap with a putative location of intererest
# --- long-term goal: make this take sets of locations
get_pca   <- function( ROI_dat_lin      = stop("GRanges must be provided"),
                       shouldplot       = TRUE,
                       PlotTitle        = "Current data",
                       cluster          = NULL
                      )
{
  Nreads      = dim( ROI_dat_lin )[1]
  Nplotpoints = dim( ROI_dat_lin )[2]

  # grab the principle components
  pca <- prcomp( ROI_dat_lin,
                 center = FALSE,
                 scale  = FALSE)

  # if you do pca <- prcomp ( A )
  # A_copy <- pca$rotation %*% t(pca$x)
  # ^ gives you back the original value "A"

  # to get the matrix that performs the forward operation (i.e. rotates from original space _into_ PCA space):
  # rot_in2_pca_mat = solve ( pca$rotation )
  # ("solve" in R means "invert", because R is terrible.)

  # pcax_copy <- t( rot_in2_pca_mat %*% t( ROI_dat_lin ) )
  # ^ This will be equal to pca$x (within machine tolerance)
  # synthdat_normed_rot2pca = rot_in2_pca_mat %*% synthdat_normed
  # these are now _row_ -based vectors

  if( shouldplot )
  {
    if( is.null( cluster ) )
      { stop("Attempting to plot with null cluster data. Exiting.") }
    plot_pca_clusters( pca       = pca,
                       cluster   = cluster,
                       PlotTitle = PlotTitle )
  }
 return( pca )  
}

# ===============================================================
# --- plot pca with cluster info:
plot_pca_clusters <- function( pca        = stop("pca dat must be provided"),
                               cluster    = NULL,
                               PlotTitle  = "Current data"
                              )
{
  Nreads <- dim( pca$x )[1]
  
    plot( # pca$x[,1] * (1/sqrt(Nxpos)),
          # pca$x[,2] * (1/sqrt(Nxpos)),
          pca$x[,1],
          pca$x[,2],
          col  = "red",
          xlab = "PC 1",
          ylab = "PC 2",
          main = PlotTitle )

    # draw a unit circle to represent a single std. dev.
    pi=3.14159265358979;
    angles = seq(0, 2*pi, 0.001);
    xcirc  = cos(angles);
    ycirc  = sin(angles);
    lines(   xcirc,    ycirc, col="black", lty = 1)
    lines( 2*xcirc,  2*ycirc, col="black", lty = 2)
    lines( 3*xcirc,  3*ycirc, col="black", lty = 3)

    if( ! is.null( cluster ))
      { 
      if( length( cluster ) != Nreads )
      { stop("cluster array supplied to get_pca incompatible with number of observed reads.") }
      # else:
      points( # pca$x[,1] * (1/sqrt(Nxpos)),
             # pca$x[,2] * (1/sqrt(Nxpos)),
             pca$x[ cluster == 1 ,1],
             pca$x[ cluster == 1 ,2],
             col  = "blue",
             lw   = 3 )
      }
  

  return( NULL )
}
# ===============================================================
# --- collect sum of squared diffs of a cluster:
collect_sum_squared_diff <- function( data        = stop("data set must be provided"),
                                      Nplotpoints = stop("Nplotpoints must be provided")
                                      )
{
  centre       <- get_cluster_mean ( data_in     = data,
                                     Nplotpoints = Nplotpoints)
  
  if ( is.null(dim( centre) ) )
    {
    vector_diffs <- rep(0,Nplotpoints)
    } else {
  
      vector_diffs <- sweep( data,
                             2,
                             centre )
    }

  result <- sum( vector_diffs ^2 )

  return( result )
}
# ===============================================================
# --- get the mean of a cluster, and handle edge-cases:
get_cluster_mean <- function( data_in     = stop("data_input required"),
                               Nplotpoints = stop("Nplotpoints must be provided") 
                               )
{
if( is.null( dim( data_in )))
   { # matrix is a single row. 
    result  = data_in } else if( dim( data_in )[1] == 0 ) {
    result = rep(0,Nplotpoints) # must ensure that dim() of this vector returns NULL
   } else {
     # matrix has at least two entries:
     result = colMeans( data_in  )
     }
   
return( result)
}

# ===============================================================
# --- filter for coverage
filter_for_coverage  <- function( plotdat_in    = stop("data_input required"),
                                  mincov        = stop("mincov must be provided"),
                                  remove_partial_coverings = FALSE
                               )
{

coverage_array <- unlist( lapply( c(1: length( plotdat_in )), 
                            function(locus) 
                              dim( plotdat_in[[locus]]$read_normdiff )[1]
                            )
                    )
result_covered <- plotdat_in[ which (coverage_array >= mincov ) ]

}

# ===============================================================
# --- collect centres and handles nulls
get_centres  <- function( x = stop("x must be provided"),
                          which_centre = stop("must specify which centre to take")
                               )
{
  if( is.null( x ))
  {return (NULL) } else{
    
    return( sqrt( 
                 rowSums( x$clusterdat$centers^2 )[which_centre] 
                 ) 
            )
    }
}
# ===============================================================
# --- collect centres and handles nulls
get_Silh <- function( clusterdat_in = stop("data_input required")   )
{
  if( is.null(clusterdat_in) ) {
    return (NULL) } else {
      return( clusterdat_in$Silh_score) 
      } 
}
# ===============================================================
# --- Collate the plot data for the higher-order:
gather_metaplot_clusterdat  <- function( clusterdat_in = stop("data_input required") 
                               )
  {
  #  null_spots <- which( unlist( lapply( clusterdat_in, function(x) is.null(x) )) )
  #  clusterdat_nullsremoved <- clusterdat_in[ -c(null_spots) ]
  #  clusterdat <- clusterdat_nullsremoved
  
  Flags        = unlist( lapply( clusterdat_in, 
                               function(x)  x$Flag  ) )
  
  filter_passed <- which( nchar( Flags ) == 0 )
  
  # take the set without flags
  clusterdat <- clusterdat_in[ filter_passed ]
 
  Silh <- lapply( clusterdat, 
                            function(x) get_Silh(x)
                        )
                
  c1_means     = lapply( clusterdat, 
                               function(x)  get_centres( x, 1) )
                       
  c2_means     = lapply( clusterdat, 
                               function(x)  get_centres( x, 2) )
  
  return( list( "Silh"      = Silh,
                "c1_means"  = c1_means,
                "c2_means"  = c2_means,
                "Flags"     = Flags
              )
          )
  }
