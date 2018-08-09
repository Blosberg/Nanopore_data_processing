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
  
breakset = seq(from = mincurrent, 
               to   = maxcurrent, 
               by   = res ) 

if( scale ) # scale by (difference from model_mean)/model_stddev
  {
  stat_params = pore_model[ seq_spec_list$seq, ]
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