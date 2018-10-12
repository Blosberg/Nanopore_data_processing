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
# --- get arbitrary sequences from an established reference

get_refgen_seqs <-  function(  refgen     = stop("seq  must be provided"),
                               ROI_GR     = stop("Region of interest must be provided in GRanges format"),
                               lead       = 0,
                               trail      = 0,
                               RNAstrand  = FALSE
)
{
  if ( as.character(strand(ROI_GR) ) == "+")
  {

    lo_end   <- start(ROI_GR) - lead;
    hi_end   <- end(ROI_GR)   + trail;
    flipseq = FALSE;

  }else if( as.character(strand(ROI_GR) )  == "-") {

    lo_end   <- start(ROI_GR) - trail;
    hi_end   <- end(ROI_GR)   + lead;    # it will be flipped later.
    flipseq = TRUE;

  } else {
    stop("undefined strand in get_refgen_seqs");
  }

  result = eval (
                parse(
                     text=paste0(
                                "refgen$",seqnames(ROI_GR),"[",as.character(lo_end),":",as.character(hi_end),"]"
                                )
                    )
                )

  # ------------------------------
  if ( flipseq )
  { # --- if the read is from the - strand, then reverse and compliment it
    result = chartr("ATGC","TACG", stringi::stri_reverse( result ) )
  }

  if( RNAstrand )
  {
    result = sub( "T", "U", result )
  }

  return( result )
}


# ==================================================================
# --- Go through a set of read sequences and reverse the order
# --- (and kmer) of only the ones aligned to the reverse strand

strand_align <- function( GRin = stop("GRobj must be provided"),
                          k = 6)
{ # k is the length of the kmer sequence (minION devices spit out k=5, but we can extend that further.)
  temp = GRin[strand(GRin)=="-"]

  GRin[ strand(GRin)=="-"]  <- strand_reversal( temp )

  return(GRin)
}
# ==================================================================
# --- reverse the strand of a given read

strand_reversal  <- function( range_in = stop("range to be flipped must be provided"),
                              k = 5 )
{
  temp                <- range_in
  temp$reference_kmer <- reverse( chartr("ATGC","TACG", range_in$reference_kmer ) )

  if( !is.null( temp$modposition) ) # i.e. check if the metadata column exists.
    {
    temp$modposition    <- (k+1) -range_in$modposition
    }

  return(temp)
}

# ==================================================================
# --- get list (separated by position of modification) of differences
# --- between current with modified base vs. control. Each list entry
# --- tabulates

get_moddiff_vs_modpos_seq <- function ( GRmod      = stop("GRmod must be provided"),
                                        GRcontrol  = stop("GRcontrol must be provided"),
                                        k=5,
                                        sigthresh = 5 )
{
  # --- temporary kludge fix for strand alignment. Figure out a more efficient way to get both strands.
  GRcontrol <- GRcontrol[strand(GRcontrol)=="+"]

  SEQHITS = sort( unique(GRmod$reference_kmer) )
  Nunique_sequences = length( SEQHITS )

  i=1
  result = matrix(0,Nunique_sequences, 3  )
  row.names(result) <- SEQHITS

  GRmodL_modpos_split <- split(GRmod, GRmod$modposition)
  GRcontrolL_seqsplit <- split(GRcontrol, GRcontrol$reference_kmer)

  difflist = list()

  for ( i in c (1:k))
    {
    temp = split( GRmodL_modpos_split[[i]], GRmodL_modpos_split[[i]]$reference_kmer)

    m=1;n=1;
    temp2=GRangesList()
    for ( m in c (1:length(temp)))
      {
      if ( length( temp[[m]] ) > sigthresh )
        {
        temp2[[n]] <- temp[[m]]
        names(temp2)[n]<- names(temp)[m]
        n <- (n+1)
        }
      }


    difflist[[i]]   <- lapply( c(1:length(temp2)), function(x)
                               ( (mean(temp2[[x]]$event_mean) -
                                    mean(GRcontrolL_seqsplit[[ names(temp2)[x] ]]$event_mean))
     #                            /(pore_model_list[[ names(temp2)[x]  ]]$std_dev)  )
                                )
                              )
  }

  return(difflist)
}

# ==================================================================
# --- Get extended sequence from pore data
expand_range <- function ( GRin       = stop("GRmod must be provided"),
                           refGenome  = stop("GRcontrol must be provided"),
                                    k = 5 )
{

  GR_justlocs = GRanges( seqnames = seqnames(GRin),
                         ranges   = IRanges(start = start(GRin), ))



  if ( k <= 0)
  {stop("k must be a positive integer")}
  # --------

  if ( as.character(strand(GRin) ) == "+")
  {
    lo_end   <- start(GRin) ;
    hi_end   <- start(GRin) -1 + k;
    flipseq = FALSE;

  }else if( as.character(strand(ROI_GR) )  == "-") {

    lo_end   <- start(ROI_GR) - trail;
    hi_end   <- end(ROI_GR)   + lead;    # it will be flipped later.
    flipseq = TRUE;

  }

  result = eval (
    parse(
      text=paste0(
        "refGenome$",seqnames(GRin),"[",as.character(lo_end),":",as.character(hi_end),"]"
      )
    )
  )

 return(result)
}


# ==================================================================
# --- compare current distributions between two data sets for a given sequence.

compare_current_byseq <- function( control   = stop("control GR must be provided"),
                                   treatment = stop("treatment GR must be provided"),
                                   seq       = stop("sequence must be provided"),
                                   breaks_in = seq(50,150,1),
                                   color_control   = rgb(0,0,1,alpha=.4),
                                   color_treatment = rgb(1,0,0,alpha=.4)
)
{
  base_set = c("A","C","G","T")

  # see how many bases are unconstrained:
  Num_Unconstrained = length( gregexpr("N", seq)[[1]] )
  if( Num_Unconstrained >= 1 )
    {
    seqlist = subs_Nvals( seq, Num_Unconstrained )
    }

  for ( x in c(1:Num_Unconstrained) )
    {


    }


  control_seq_subset   = control[   control$model_kmer   == seq ]
  treatment_seq_subset = treatment[ treatment$model_kmer == seq ]

  hist( control_seq_subset$event_mean,
        breaks = breaks_in,
        freq   = FALSE,
        col    = color_control,
        add    = F,
        main   = paste("Current distribution for sequence", seq) )

  hist( treatment_seq_subset$event_mean,
        breaks   = breaks_in,
        freq     = FALSE,
        col      = color_treatment,
        add      = T)

  legend("topright",
         legend = c("control", "treatment"),
         fill   = c( color_control, color_treatment))

}

# ==================================================================
# --- return normalized histogram for a given set of breaks and GR
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

}


