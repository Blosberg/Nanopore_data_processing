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
                                         line_darkness    = 0.1 )
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
      lines( squiggle_dat$x,
             squiggle_dat$y,
             col    = rgb(0, 0, 1, line_darkness),
             type   = "l",
             lwd    = 1
            )
      }

  } else{
    for (read_i in c(1:Nreads) )
      {
      lines( sampleROI_dat$xposn_vals,
             sampleROI_dat$posn_means[read_i, ],
             col    = rgb(0, 0, 1, line_darkness),
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
plot_dwelltime_over_ROI <-  function( sampleROI_dat  = stop("aligned sampledat must be provided"),
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
# --- long-term goal: make this take sets of locations
get_pca_clusters  <- function( sampleROI_dat_in  = stop("GRanges must be provided"),
                               k                 = 5,
                               method            = "pca", # TODO: add a KL-div based method.
                               shouldplot        = TRUE
                              )
{
  Nxpos = dim( sampleROI_dat_in$read_normdiff )[2]

  # get the subset of data points near the centre of the ROI:
  ROI_dat_lin = sampleROI_dat_in$read_normdiff[ , c( (ceiling(Nxpos/2)-k+1) : (ceiling(Nxpos/2)+k-1) ) ]


  ROI_dat_lin[ is.na( ROI_dat_lin ) ] <- 0

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

  clust_kmeans <- kmeans( x = ROI_dat_lin,
                          centers = 2,
                          iter.max = 100 )


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

  if( shouldplot )
  {
    plot( # pca$x[,1] * (1/sqrt(Nxpos)),
          # pca$x[,2] * (1/sqrt(Nxpos)),
          pca$x[,1],
          pca$x[,2],
          col  = "red",
          xlab = "PC 1",
          ylab = "PC 2",
          main = "Deviation from model" )

    # draw a unit circle to represent a single std. dev.
    pi=3.14159265358979;
    angles = seq(0, 2*pi, 0.001);
    xcirc  = cos(angles);
    ycirc  = sin(angles);
    lines(   xcirc,    ycirc, col="black", lty = 1)
    lines( 2*xcirc,  2*ycirc, col="black", lty = 2)
    lines( 3*xcirc,  3*ycirc, col="black", lty = 3)

   points( # pca$x[,1] * (1/sqrt(Nxpos)),
           # pca$x[,2] * (1/sqrt(Nxpos)),
           pca$x[ clust_kmeans$cluster == 1 ,1],
           pca$x[ clust_kmeans$cluster == 1 ,2],
           col  = "blue",
           lw   = 3 )
  }

  Silh <- silhouette( x = clust_kmeans$cluster,
                      dist = dist(pca$x) )

  # TODO: incorporate cluster robustness

  return( list("pca"          = pca,
               "clust_kmeans" = clust_kmeans,
               "Silh"         = Silh )
         )
}

