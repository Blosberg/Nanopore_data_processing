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
get_list_of_read_olaps_with_ROIs <- function()
{ ... }
   
# ===============================================================
# --- plot the current lines from all the reads of a sample over a given region:
plot_samplesignal_over_ROI <-  function( SampleName = stop("Sample name must be provided"),
                                         overlapping_reads  = stop("reads must be provided"),
                                         refgen = stop("Reference genome must be provided"),
                                         ROI_in = stop("Region of Interest Grange must be provided"),
                                         standard_dat = stop("standard reference data must be provided."),
                                         k=5,
                                         mincurrent = 50,
                                         maxcurrent = 150 )
{ 
  # for these inputs, we know that the reads overlap at least once on the ROI, but we don't know where
  Nreads = length( overlapping_reads )
  
  # get number of x positions in the plotting domain: 
  Nxpos     = end( ROI_in ) - start( ROI_in ) + 1;
  xvals     = c( start( ROI_in):end(ROI_in) )
  chr       = as.character( unique( seqnames( ROI_in ) ) )

  
  # Get standard plot
  # get_reference_strand # write a function here: 
  # need to take reverse-compliment of ref genome for "-" strands
  # Also need to transliterate "T" to "U" for RNA, 
  # use lapply to take ref values for mean and stand. dev.
  # ymean   = unlist ( lapply( Datin , function(x) c(x$event_level_mean, x$event_level_mean ) ) )
  # ystddev = unlist ( lapply( Datin , function(x) c(x$event_stdv,       x$event_stdv ) ) )

  # Then plot a grey region between the +/- 1 std. dev. Everything after that should be "lines"
  plot(  xvals, ymean,
         col="black",
         type="l",
         lty = "dashed",
         lwd=2,      
         xlim = c( start(ROI_in), 
                   end(ROI_in) ),
         
         main= paste(  SampleName, "\n",
                       chr, 
                       ":", 
                       as.character( start(ROI_in) ),
                       "-", 
                       as.character( end(ROI_in)  ) ),
         xlab = "reference position",
         xaxt="n",
         ylab = "current [pA]",
         ylim = c( mincurrent, 
                   maxcurrent )
      )
  
  # add reference sequence, complement, and directionality.
  add_sequence_tick_marks( refgen = refgen,
                           Locus  = ROI_in )
  
  # plot a transparent polygon around the +/-1 std.dev region of the standard (expected) current values.
  polygon( c( xvals, rev(xvals) ), 
           c(  (ymean-ystddev), 
               rev(ymean+ystddev) ),
           
           col  = rgb(0, 0, 0, 0.1),
           lty  = "blank"
          )

  #========  NOW PLOT THE READS THEMSELVES =========
  # Get the particular position within each read where it overlaps with the ROI

  # populate a matrix of current values at each position.  
  plotdat <- matrix( NA, Nreads, Nxpos )
  for (i in c(1:Nreads) )
    {
    xivals <- start( overlapping_reads[[i]] ) - xvals[1] + 1;
    plotdat[ i, xivals ] <-  overlapping_reads[[i]]$event_level_mean;
    # positions in each read without current data are left "NA"
  # plot the rows of this matrix for each read
  lines( xvals, 
         plotdat[i, ],
         col    =rgb(0,0,1,0.5),
         type   = "l",
         lwd    = 1
       )
  }

}
# -------------------------------------------------------------

add_sequence_tick_marks   <- function(  refgen        = stop("Reference genome must be defined"),
                                        Locus         = stop("Locus must be defined"),
                                        Col_sense     = rgb( 0, 0, 0, 1  ),
                                        Col_antisense = rgb( 0, 0, 0, 0.5),
                                        mincurrent    = 50,
                                        maxcurrent    = 150,
                                        RNA_strand    = TRUE
                                        )
{
  # chromosome coordinates, and sequence values:
  xvals             = c( start( Locus):end(Locus) )
  chrom_seq         = eval( parse( text=paste0( "refgen$", as.character( unique( seqnames(Locus))) ) ) )
  seqchars_plus     = chrom_seq[ xvals ]
  seqchars_neg      = chartr("ACGT", "TGCA", seqchars_plus ) 
  
  # Now label give plus/neg strand appropriate emphasis:
    if ( strand( Locus ) == "+" )
  { x0 = start(  Locus )
    x1 = end(Locus)
    Col_plusstrand = Col_sense 
    Col_negstrand  = Col_antisense  } else if( strand(Locus) =="-" ){
      
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
