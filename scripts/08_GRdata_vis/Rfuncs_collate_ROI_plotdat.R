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
                                      mincov = 10,
                                      OLAP_skip_TOL = OLAP_skip_TOL )
{
#  if( ! identical( names(reads), names( loci ) ) )
#  { stop ("ERROR: inconsistent names between reads and loci") }

# expand loci slightly to ensure overlap is registered,
# even if the exact base is skipped.
loci_expanded <- loci

if ( identical( all( width ( loci )  < OLAP_skip_TOL ) , TRUE ) )
        {
        start( loci_expanded ) <-  ( start( loci ) - OLAP_skip_TOL )
        end(   loci_expanded ) <-  ( end(   loci ) + OLAP_skip_TOL )
        }

overlaps = findOverlaps(  loci_expanded,
                          reads )
rm( loci_expanded )

# collect a list of how many hits each loci accumulated:
loci_coverage  =  collect_ROI_coverage ( loci         = loci,
                                         ROI_overlaps = overlaps  )

# collect a list of which loci are above threshold coverage:
covered_loci_list <-  which( loci_coverage >= mincov )

# now filter those above threshold:
loci_filtered_for_coverage <- loci[ covered_loci_list ]

return( loci_filtered_for_coverage )

}
# ===============================================================
# --- collect coverage:
collect_ROI_coverage <-  function( loci          = stop("loci must be provided"),
                                   ROI_overlaps  = stop("olaps must be provided"))
{
  result <- unlist( lapply( c(1:length( loci )),
                            function(x) length( which( queryHits( ROI_overlaps ) == x ) )
  )
  )
  return( result )
}

# ===============================================================
# --- plot the current lines from all the reads of a sample over a given region:
get_sampledat_over_ROI <-  function( overlapping_reads  = stop("reads must be provided"),
                                     refgen             = stop("Reference genome must be provided"),
                                     ROI_raw            = stop("Region of Interest Grange must be provided"),
                                     poremodel          = stop("standard reference data must be provided."),
                                     k_in               = 5,
                                     plotrange          = 10,
                                     mincurrent         = 50,
                                     maxcurrent         = 150,
                                     perform_sanity_checks = FALSE )
{
  # for these inputs, we know that the reads overlap at least once on the ROI, but we don't know where
  Nreads             <- length( overlapping_reads )
  plotRegion         <- ROI_raw
  start( plotRegion) <- start(plotRegion) - plotrange - k_in;
  end(   plotRegion) <- end(  plotRegion) + plotrange + k_in;

  # get number of x positions in the plotting domain:
  Nxpos      = end( plotRegion ) - start( plotRegion ) + 1;
  xposn_vals = c( start( plotRegion):end(plotRegion) )
  
  if( perform_sanity_checks )
    {
    chr        = as.character( unique( seqnames( plotRegion ) ) ) # makes sure they're the same.
    } else {
     chr       = as.character( seqnames( plotRegion[1] )  )
    }

  # Get standard plot
  # get_reference_strand # write a function here:
  # need to take reverse-compliment of ref genome for "-" strands
  poremodel_metrics <- extract_poremodel_metrics_from_sequence ( poremodel_in  = poremodel,
                                                                 refgen        = refgen,
                                                                 ROI           = plotRegion,
                                                                 chr_in        = chr,
                                                                 strandtype    = "RNA"
                                                                )
  # ----- declare and allocate output.----------
  # populate a matrix of current values at each position (averaged over all events and samples
  
  posn_means <- matrix( NA, Nreads, Nxpos )
  colnames( posn_means ) <- as.character( xposn_vals )
  if( length( overlapping_reads ) == 0 )
    { stop("ERROR: in function get_sampledat_over_ROI: overlapping_reads has length 0.") }
  row.names( posn_means ) <- as.character( paste0("read_", names( overlapping_reads ) ) )

  # Also tabulate the "squiggle" current readings (keep in character format):
  Isamples_chars              <- matrix( NA, Nreads, Nxpos )
  colnames( Isamples_chars )  <- as.character( xposn_vals )
  row.names( Isamples_chars ) <- as.character( paste0("read_", names( overlapping_reads ) ) )

  # Also tabulate an array of the "event-averaged" current values (multiple possible events per position)
  Event_means_chars  <- matrix( NA, Nreads, Nxpos )
  colnames( Event_means_chars ) <- as.character( xposn_vals )
  row.names( Event_means_chars ) <- as.character( paste0("read_", names( overlapping_reads ) ) )

  # And tabulate the duration of each of these events.
  Event_duration_chars <- matrix( NA, Nreads, Nxpos )
  colnames( Event_duration_chars ) <- as.character( xposn_vals )
  row.names( Event_duration_chars ) <- as.character( paste0("read_", names( overlapping_reads ) ) )

  # -----------------------------------------------
  # scan through each read overlapping this ROI
  
  for (i in c(1:Nreads) )
    {
    # get just the segment of the read that overlaps
    readseg_olaps <- findOverlaps( plotRegion,
                                   overlapping_reads[[i]] );
    # and then take just the segment of the read that overlaps with the plotRegion:
    readseg_plotRegion   <- overlapping_reads[[i]][ subjectHits( readseg_olaps ) ]

    # re-initialize:
    read_posn_dat <- list()
    read_posn_dat <- collect_singlebase_res_read_data ( readseg_plotRegion = readseg_plotRegion,
                                                        xposn_vals  = xposn_vals )

    # now pass these values over to the matrix by x-position
    # assign corresponding values to the matrix (with possible NA breaks).
    posn_means[ i,  ]            <- read_posn_dat$posn_means[1: Nxpos]
          # numeric; single value for each position.
    Isamples_chars[ i,  ]        <- read_posn_dat$Isamples_chars[1: Nxpos]
          # string; each position has an array of sample-values, separated by ","
    Event_means_chars [ i,  ]    <- read_posn_dat$Event_means_chars[1: Nxpos]
          # string; each position as an array of means for the corresponding event.
          # Multiple possible events per position --separated by ","
    Event_duration_chars [ i,  ] <- read_posn_dat$Event_duration_chars[1: Nxpos]
          # string; each position has an array of duration values.
          # should be same length as above --separated by ","
  }

  # ========  SANITY CHECK =======================
  if( perform_sanity_checks )
    {
    posn_means_xcheck = matrix( NA, Nreads, Nxpos )
    xcheck_compare    = matrix( NA, Nreads, Nxpos )
    for (i in c(1:Nreads) )
      {
      for (j in c(1:Nxpos ) )
        {
        if ( !is.na( Isamples_chars   [i,j ] ) )
          { posn_means_xcheck[i,j] <- mean( extract_numericArray_from_charlist(Isamples_chars   [i,j ])) }
        }
      }

    xcheck_compare = (posn_means_xcheck - posn_means)/posn_means
    xcheck_compare[ is.na(posn_means)  ] <- 0
    if( max ( xcheck_compare) > 0.01 )
      {stop("current mean differs from sample mean by more than 1%.")}
  }
  # ========  SANITY CHECK COMPLETE ==============

  # for each read, take the difference from the standard model mean.
  read_diff_from_model <- sweep( posn_means,
                                 2,
                                 as.array( poremodel_metrics[,1]) )
  colnames( read_diff_from_model )  <- colnames( posn_means )
  row.names( read_diff_from_model ) <- row.names( posn_means )

  # and then divide it by the std. dev.
  read_normdiff  <- sweep( read_diff_from_model,
                           2,
                           as.array( poremodel_metrics[,2] ),
                           FUN = "/" )

  # return the read-data in matrix form.
  return( list( "ROI"                  = ROI_raw,
                "plotRegion"           = plotRegion,
                "xposn_vals"           = xposn_vals,
                "poremodel_metrics"    = poremodel_metrics,
                "read_normdiff"        = read_normdiff,
                "posn_means"           = posn_means,
                "Isamples_chars"       = Isamples_chars,
                "Event_means_chars"    = Event_means_chars,
                "Event_duration_chars" = Event_duration_chars)
          )
}
# ===============================================================
# --- gather a list of matrices (of fixed size) that tabulate base-position specific data for each read.
collect_singlebase_res_read_data <-function( readseg_plotRegion = stop("readseg_plotRegion must be provided"),
                                             xposn_vals  = stop("xposn_vals must be provided"))
{
  Nxpos  <- length( xposn_vals )
  # Nreads == 1 (we perform this function on each read, just using the segment that overlaps with ROI)

  # populate a matrix of current values at each position (averaged over all events and samples
  posn_means <- matrix( NA, 1, Nxpos )
#  colnames( posn_means )          <- as.character( xposn_vals )

  # Also tabulate the "squiggle" current readings (keep in character format):
  Isamples_chars  <- matrix( NA, 1, Nxpos )
#  colnames( Isamples_chars )      <- as.character( xposn_vals )

  # Also tabulate an array of the "event-averaged" current values (multiple possible events per position)
  Event_means_chars  <- matrix( NA, 1, Nxpos )
#  colnames( Event_means_chars )   <- as.character( xposn_vals )

  # And tabulate the duration of each of these events.
  Event_duration_chars <- matrix( NA, 1, Nxpos )
#  colnames( Event_duration_chars ) <- as.character( xposn_vals )

  for ( x in c(1:Nxpos))
    {
    readseg_at_bploci = readseg_plotRegion[ start( readseg_plotRegion ) == xposn_vals[x] ]

    # Depending on convention as to sequence "location", this may need to be offset by one.
    # At latest check (July 2019), the convention below is consistent with observations:
    posn_means[x]           = sum( readseg_at_bploci$event_mean * readseg_at_bploci$event_length )/sum(readseg_at_bploci$event_length)
    Isamples_chars[x]       = paste( readseg_at_bploci$samples,
                                  collapse = "," )

    Event_means_chars[x]    = paste( readseg_at_bploci$event_mean, collapse = "," )
    Event_duration_chars[x] = paste( readseg_at_bploci$event_length, collapse = "," )
    }

return( list(
            "posn_means"           = posn_means,
            "Isamples_chars"       = Isamples_chars,
            "Event_means_chars"    = Event_means_chars,
            "Event_duration_chars" = Event_duration_chars
            )
      )

  }

# ===============================================================
# --- Get the "Standard" current curve from the reference genome sequence:
extract_poremodel_metrics_from_sequence   <- function ( poremodel_in  = stop("poremodel data required"),
                                                        refgen        = stop("refgen data required"),
                                                        ROI           = stop("ROI data required"),
                                                        chr_in        = stop("seqname must be provided"),
                                                        strandtype    = "RNA",
                                                        k = 5
                                                       )
{
  chr = chr_in
  if( as.character( strand(ROI)) == "+" )
    { RefSequence = get_sequence( refgen = refgen,
                             chr    = chr,
                             range  = c( start(ROI)+1, end(ROI)+k ) )
    } else if( as.character( strand(ROI)) == "-" ){
      RefSequence = reverseComplement( get_sequence( refgen = refgen,
                                                   chr    = chr,
                                                   range  = c( start(ROI)+1, end(ROI)+k ) )
                                      )
    } else{
     stop("invalid strand.")
    }

  if( length(RefSequence)-k < 1 )
    { stop("ERROR: invalid RefSequence")}

  # transliterate the "poremodel" data, as needed, to compare to DNA reference genome.
  poremodel_ref = poremodel_in
  if( strandtype =="RNA" )
     { row.names(poremodel_ref) <- gsub("U", "T", row.names(poremodel_ref) )  }
  # now poremodel_ref is always written with "T"s to compare to reference,
  # regardless of the initial input model poremodel_in

  # now extract the right row from this poremodel reference for each position in the ROI
  std_current=matrix( 0,
                      length(RefSequence)-k+1,
                      2 )

  # N.B. Refsequence is already reverse-complimented (as necessary) according to strand.
  seq_kmers       <-  unlist( lapply( c(1:(length(RefSequence)-k+1)), function(x)
                                     as.character( substr( RefSequence, x, (x+k-1) ) ) ) )

  std_current[,1] <- unlist( lapply( c(1:(length(RefSequence)-k+1)), function(x)
                                     poremodel_ref[ seq_kmers[x] , 1 ]
                                     ) )
  std_current[,2] <- unlist( lapply( c(1:(length(RefSequence)-k+1)), function(x)
                                     poremodel_ref[ seq_kmers[x], 2 ]
                                     ) )

  # reverse order of values, if we are on the reverse strand.
  if( as.character( strand(ROI)) == "-" ) {
    std_current[,1] = rev( std_current[,1] )
    std_current[,2] = rev( std_current[,2] )
    }

  colnames( std_current )   <- c( "mean", "stddev")
  row.names( std_current )  <- seq_kmers

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
# ===============================================================
# --- For a given "group" of reads/RsOI, collect read data in a given range nearby
collect_sample_dat_over_ROIs <- function ( ref_Genome     = stop("refGenome must be supplied."),
                                           aligned_reads  = stop("reads must be supplied"),
                                           ROI_overlaps   = stop("overlaps must be specified"),
                                           loci_covered   = stop("loci must be provided"),
                                           poremodel_in   = stop("poremodel must be specified."),
                                           plotrange_in   = 10,
                                           k_in           = 5
)
{
  ROIread_dat = lapply( c(1:length( loci_covered)),
                              function(locus_i)
                                get_sampledat_over_ROI ( overlapping_reads  = aligned_reads[
                                                              subjectHits( ROI_overlaps[
                                                                  queryHits( ROI_overlaps ) == locus_i ] ) ],
                               # To understand the above lines.
                               # for plotting, we want to collect:
                               #                         ^ The subset of the reads
                               #                              ^ That are referenced as the subject
                               #                                  ^ in an overlap-pair for which the query is the ROI we are looking at (i.e. ="locus_i").

                               refgen             = ref_Genome,
                               ROI_raw            = loci_covered[locus_i],
                               plotrange          = plotrange_in,
                               poremodel          = poremodel_in,
                               k_in               = k_in
                                                       )
                            )

  return(  ROIread_dat )
  }

