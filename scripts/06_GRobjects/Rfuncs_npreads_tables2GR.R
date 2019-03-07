# ========================================================================
# --- Extract strand information based on relationship between event_indices and position. 

assign_strand <-  function( Datin = stop("Datin must be provided") )
{ 
  
  Datin$strand = NULL
  
  dtest = Datin[c("read_index", "position", "event_index") ] %>%
    #  subset(read_index %in% 1:10) %>%              # minimal check that things are working.
    group_by(read_index) %>%                         # First partition by read_index
    # mutate(position_index = order(position)) %>%     # Set position_index variable defined as the order of 
    # "position" (lowest to highest, with equal values taken 
    # in order of row number -- all within each read_index)

    mutate(  diff_event    = c(0, diff( event_index ) ) ) %>%
    mutate(  diff_position = c(0, diff( position    ) ) ) %>%   # adds diff_index col which is order of event - order of position 

    mutate(strand = case_when( 
      all( diff_event * diff_position >= 0 ) ~ '+',
      all( diff_event * diff_position <= 0 ) ~ '-',
      TRUE ~ '*'
    ) )
  
  dtest$diff_event    <- NULL
  dtest$diff_position <- NULL
  
  ftest <- left_join(Datin, dtest, by=c('read_index','position','event_index') )
  
  return(ftest)
  
}
# ========================================================================
# import table-event data & make sure there's a single-header, and each read has a sequential, unique read_index.

get_event_dat  <-  function( Event_file_list = stop("Datin must be provided") 
                             )
{
  colnames <- c( "contig", "position", "reference_kmer",  "read_index", "strand", "event_index", 
                 "event_level_mean", "event_stdv", "event_length", "model_kmer", "model_mean", 
                 "model_stdv", "standardized_level" )
  dat_all = data.frame( row.names = colnames)
  readcount_offset = 0;
  
  for ( i in c(1:length(Event_file_list) ) )  {
  
    # Event_file      = paste0(    "Ealign_", as.character(i),".cvs")
    fin               = as.character( file.path(  Event_file_list[i] ) )

    write(  paste( "---Reading csv file from: ", fin ), 
            file   = logFile,
            append = TRUE )

    dat_temp        = read.csv(  file = fin, 
                                 sep  = '\t', 
                                 stringsAsFactors=FALSE, 
                                 header = TRUE
                              ) #--- read in the current "chunk" of np data
    
    dat_temp$read_index <- dat_temp$read_index + readcount_offset  #--- offset the read_index values to ensure uniqueness
    
    dat_all           =    rbind(dat_all, dat_temp )               #--- compile the chunks together
    readcount_offset  =  ( max( dat_all$read_index)  +1  )         #--- re-calculate the necessary offset
    
  }
  
  Nreads = length(unique(dat_all$read_index))
  
  # cast the read_indices as factors (to make unique values sequentially increasing), and then 
  # cast them back into integers (since some read_indices have been removed for quality/alignment/etc. reasons
  dat_all$read_index <- as.numeric( as.factor( dat_all$read_index))
  
  
  return(dat_all)
}


# -------------------------------------------------------------------------------
plot_signal <-  function( Datin = stop("Datin must be provided") )
{ 
N=length( Datin )
  
xvals   = c( 0,  sort( c( c(1:(N-1)) ,c(1:(N-1)) ) ), N )
ymean   = unlist ( lapply( Datin , function(x) c(x$event_level_mean, x$event_level_mean ) ) )  
ystddev = unlist ( lapply( Datin , function(x) c(x$event_stdv,       x$event_stdv ) ) )  


plot(  xvals, ymean, 
       col="blue", 
       type="l", 
       lwd=2, 
       main = "current vs. event",
       xlab = "event index",
       ylab = "current" )

lines( xvals, ymean-ystddev, col="blue", type="l", lwd=0.5 )
lines( xvals, ymean+ystddev, col="blue", type="l", lwd=0.5 )

}
