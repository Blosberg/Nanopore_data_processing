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

# ===============================================================
# Get the loci that are adequately covered, and the reads that cover them
filter_loci_for_coverage <- function( loci   = stop("Loci list must be provided"),
                                      reads  = stop("readlist must be provided"),
                                      mincov = 10)
{
overlaps_by_group = lapply ( names( loci ), function(x)
                               findOverlaps(  loci[[x]],
                               reads[[x]] )
                               )
names( overlaps_by_group ) <- names( loci )

# collect a list of how many hits each loci accumulated:
loci_coverage_by_group  = lapply( names( loci ),
                                  function(group) collect_ROI_coverage ( group_loci     = loci[[group]],
                                                                         group_overlaps = overlaps_by_group[[group]]  )
)
names( loci_coverage_by_group ) <- names(loci)


# collect a list of which loci are above threshold coverage:
covered_loci_list <- lapply( loci_coverage_by_group,
                                  function(group) which( group >= mincov ) )

# now filter those above threshold:
loci_filtered_for_coverage <- lapply( names( loci ),
                                      function(group) loci[[group]][ covered_loci_list[[group]] ] )
names( loci_filtered_for_coverage ) <- names(loci)

return( loci_filtered_for_coverage )
}
# ===============================================================
# --- collect coverage:
collect_ROI_coverage <-  function( group_loci     = stop("group loci must be provided"),
                                   group_overlaps = stop("group olaps must be provided"))
{
  result <- unlist( lapply( c(1:length( group_loci )),
                            function(x) length( which( queryHits( group_overlaps ) == x ) )
  )
  )
  return( result )
}

# ===============================================================
# for efficiency, implement this function to screen overlaps and reads, as well as loci.
# # --- Filter olaps and reads:
# filter_olaps_and_reads <- function()
# {
# # find out which overlaps correspond to these loci:
# covered_overlaps_by_group = lapply( names( loci ),
#                                     function(group) overlaps_by_group[[group]][
#
#                                       which( ! is.na( match(  queryHits(  overlaps_by_group[[group]] ),
#                                                               covered_loci_list[[group]] )
#                                       )
#                                       ) ]
# )
# names( covered_overlaps_by_group ) <- names(loci)
#
# }
# ===============================================================
# --- plot the current lines from all the reads of a sample over a given region:
get_sampledat_over_ROI <-  function( SampleName         = stop("Sample name must be provided"),
                                     overlapping_reads  = stop("reads must be provided"),
                                     refgen             = stop("Reference genome must be provided"),
                                     ROI_raw            = stop("Region of Interest Grange must be provided"),
                                     poremodel_ref      = stop("standard reference data must be provided."),
                                     k                  = 5,
                                     plotrange          = 10,
                                     mincurrent         = 50,
                                     maxcurrent         = 150 )
{
  # for these inputs, we know that the reads overlap at least once on the ROI, but we don't know where
  Nreads     = length( overlapping_reads )
  ROI        = ROI_raw
  start(ROI) <- start(ROI_raw) - plotrange;
  end(ROI)   <- end(ROI_raw)   + plotrange;

  # get number of x positions in the plotting domain:
  Nxpos      = end( ROI ) - start( ROI ) + 1;
  xposn_vals = c( start( ROI):end(ROI) )
  chr        = as.character( unique( seqnames( ROI ) ) )

  # Get standard plot
  # get_reference_strand # write a function here:
  # need to take reverse-compliment of ref genome for "-" strands
  poremodel_metrics <- extract_poremodel_metrics_from_sequence ( poremodel_ref = poremodel_ref,
                                                                 refgen        = refgen,
                                                                 ROI           = ROI,
                                                                 strandtype    = "RNA"
                                                                )

  # populate a matrix of current values at each position.
  read_currents <- matrix( NA, Nreads, Nxpos )

  # and grab dwell times while we're at it.
  dwell_times <- matrix( NA, Nreads, Nxpos )

  for (i in c(1:Nreads) )
    {
    # get just the segment of the read that overlaps
    readseg_olaps <- findOverlaps( ROI,
                                   overlapping_reads[[i]] );
    # and then take just the segment of the read that is within the ROI:
    readseg_ROI  <- overlapping_reads[[i]][ subjectHits( readseg_olaps ) ]

    # assign corresponding values to the matrix (with possible NA breaks).
    read_currents [i, start( readseg_ROI ) - xposn_vals[1] + 1 ]     = readseg_ROI$event_mean

    dwell_times[i, start( readseg_ROI ) - xposn_vals[1] + 1 ] = readseg_ROI$event_length
  }

  # return the read-data in matrix form.
  return( list( "SampleName"         = SampleName,
                "ROI"                = ROI,
                "xposn_vals"         = xposn_vals,
                "poremodel_metrics"  = poremodel_metrics,
                "read_currents"      = read_currents,
                "dwell_times"        = dwell_times) )

}

# ===============================================================
# --- plot the current lines from all the reads of a sample over a given region:
plot_samplesignal_over_ROI<-  function( sampleROI_dat  = stop("aligned sampledat must be provided"),
                                        refgen         = stop("refgen must be provided"),
                                        mincurrent     = 50,
                                        maxcurrent     = 150 )
{
  Nxpos = length( sampleROI_dat$xposn_vals)
  Nreads = dim( sampleROI_dat$read_currents )[1]

  poremodel_Icurves = matrix( 0,
                              Nxpos,
                              3 )
  # now fill in the min/mean/max values into three columns of the output respectively.
  poremodel_Icurves[,1] = sampleROI_dat$poremodel_metrics[,1] - sampleROI_dat$poremodel_metrics[,2]
  poremodel_Icurves[,2] = sampleROI_dat$poremodel_metrics[,1]
  poremodel_Icurves[,3] = sampleROI_dat$poremodel_metrics[,1] + sampleROI_dat$poremodel_metrics[,2]


  # Then plot a grey region between the +/- 1 std. dev. Everything after that should be "lines"
  plot(  sampleROI_dat$xposn_vals,
         poremodel_Icurves[ , 2],
         col=rgb(0, 0, 0, 0.5),
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

  # add reference sequence, complement, and directionality.
  add_sequence_tick_marks( refgen = refgen,
                           Locus  = sampleROI_dat$ROI )

  # plot a transparent polygon around the +/-1 std.dev region of the standard (expected) current values.
  polygon( c( sampleROI_dat$xposn_vals, rev(sampleROI_dat$xposn_vals) ),
           c(  (poremodel_Icurves[ , 1]),
               rev(poremodel_Icurves[ , 3]) ),
           col  = rgb(0, 0, 0, 0.4),
           lty  = "blank"
          )

  #========  NOW PLOT THE READS THEMSELVES =========
  # Get the particular position within each read where it overlaps with the ROI


    # and plot the mean current vals over this range:
    # positions in each read without current data are left "NA"
  # plot the rows of this matrix for each read
  for (i in c(1:Nreads) )
    {
    lines( sampleROI_dat$xposn_vals,
           sampleROI_dat$read_currents[i, ],
           col    = rgb(0,0,1,0.5),
           type   = "l",
           lwd    = 1
          )
    }

}
# ===============================================================
# --- plot the passage time of a sample over a given region:
plot_dwelltime_over_ROI <-  function( sampleROI_dat  = stop("aligned sampledat must be provided"),
                                      refgen         = stop("Reference genome must be provided"),
                                      log            = FALSE )
{
  Nxpos  = length( sampleROI_dat$xposn_vals)
  Nreads = dim( sampleROI_dat$dwell_times )[1]

  dwell_dat = sampleROI_dat$dwell_times
  if( log )
    { dwell_dat[!is.na(dwell_dat)] <- log( dwell_dat[!is.na(dwell_dat)]  ) }

  min_T          = min( dwell_dat, na.rm = TRUE )
  max_T          = max( dwell_dat, na.rm = TRUE )


  plot(  sampleROI_dat$xposn_vals,
         dwell_dat[1, ],
         col=rgb(0, 0, 0, 0.5),
         type="l",
         lwd=1,
         xlim = c( start(sampleROI_dat$ROI),
                   end(sampleROI_dat$ROI) ),

         main= paste(  SampleName, "\n",
                       chr,
                       ":",
                       as.character( start(ROI) ),
                       "-",
                       as.character( end(ROI)  ) ),
         xlab = "reference sequence",
         xaxt = "n",
         ylab = "dwell time [s]",
         ylim = c( min_T,
                   max_T)
      )

  # add reference sequence, complement, and directionality.
  add_sequence_tick_marks( refgen = refgen,
                           Locus  = sampleROI_dat$ROI )

    # and plot the mean current vals over this range:
    # positions in each read without current data are left "NA"
  # plot the rows of this matrix for each read
  for (i in c(2:Nreads) )
    {
    lines( sampleROI_dat$xposn_vals,
           dwell_dat[i, ],
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
                                        mincurrent    = 50,
                                        maxcurrent    = 150,
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
  arrows( x0, mincurrent, x1, mincurrent,
          length = 0.25, angle = 20,
          code = 2, col = par("fg"), lty = par("lty"),
          lwd = 3 )

}
# ===============================================================
# --- Get the "Standard" current curve from the reference genome sequence:
extract_poremodel_metrics_from_sequence   <- function ( poremodel_ref = poremodel_ref,
                                                        refgen        = refgen,
                                                        ROI           = ROI,
                                                        strandtype    = "RNA",
                                                        k = 5
                                                       )
{
  chr = as.character( unique(seqnames(ROI )) )
  if( as.character( strand(ROI)) == "+" )
    { sequence = get_sequence( refgen = refgen,
                             chr    = chr,
                             range  = c( start(ROI)+1, end(ROI)+k ) )
    } else if( as.character( strand(ROI)) == "-" ){
      sequence = reverseComplement( get_sequence( refgen = refgen,
                                                   chr    = chr,
                                                   range  = c( start(ROI)+1, end(ROI)+k ) )
                                   )

    } else{
     stop("invalid strand.")
    }

  if( length(sequence)-k < 1 )
    { stop("ERROR: invalid sequence")}

  # transliterate the "poremodel" data, as needed, to compare to DNA reference genome.
  strand_ref = poremodel_ref
  if( strandtype =="RNA" )
    { row.names(strand_ref) <- gsub("U", "T", row.names(strand_ref) )  }

  # now extract the right row from this poremodel reference for each position in the ROI
  std_current=matrix( 0,
                      length(sequence)-k+1,
                      2 )

  std_current[,1] <- unlist( lapply( c(1:(length(sequence)-k+1)), function(x)
                                     poremodel_ref[ as.character(substr( sequence, x, (x+k-1) )) , 1 ]
                                     ) )
  std_current[,2] <- unlist( lapply( c(1:(length(sequence)-k+1)), function(x)
                                     poremodel_ref[ as.character(substr( sequence, x, (x+k-1) )) , 2 ]
                                     ) )

  # reverse order of values, if we are on the reverse strand.
  if( as.character( strand(ROI)) == "-" ) {
    std_current[,1] = rev( std_current[,1] )
    std_current[,2] = rev( std_current[,2] )
    }

  colnames( std_current )   <- c( "mean", "stddev")
  row.names( std_current )  <- unlist( lapply( c(1:(length(sequence)-k+1)), function(x)
                                     as.character( substr( sequence, x, (x+k-1) ) )
                                     ) )

  return( std_current )
}
# ===============================================================
# --- Get the "Standard" current curve from the reference genome sequence:
get_sequence   <- function ( refgen  = stop("refgenome must be provided"),
                             chr     = stop("chr must be provided"),
                             range   = stop("range must be provided")
                            )
{
  if( length(range) != 2)
    stop("irregular range submitted to get_sequence")
  result = eval( parse( text=paste0( "refgen$",
                                      chr,"[",as.character(range[1]),
                                      ":",
                                      as.character(range[2]),"]" )) )
  return( result )
}
