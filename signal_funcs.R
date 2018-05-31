# ======== DEFINE FUNCS ======================================
group_by_base_position <-  function( Datin = stop("Datin must be provided"), 
                                     sequence_length = 5, 
                                     index = stop("index must be provided"),
                                     base_set = c("A","C","G","T")) 
{ 
  groupings=list()
  
  groupings[[1]] = Datin[ which( substr( Datin$reference_kmer, index,index) == "A"   ) , ]
  groupings[[2]] = Datin[ which( substr( Datin$reference_kmer, index,index) == "C"   ) , ]
  groupings[[3]] = Datin[ which( substr( Datin$reference_kmer, index,index) == "G"   ) , ]
  groupings[[4]] = Datin[ which( substr( Datin$reference_kmer, index,index) == "T"   ) , ]
  
  return(groupings)
}
# ======== DEFINE FUNCS ======================================
group_by_freq <-  function( Datin = stop("Datin must be provided"), 
                            sequence_length = 5, 
                            seq       = stop("target sequence to look for must be provided"),
                            xrange    = c(50,150),
                            break_set = seq(from = 49, to = 151, by = 2)
)
{ 
  if ( str_count(seq,"A") + str_count(seq,"C") + str_count(seq,"G") + str_count(seq,"T")  != nchar(seq) )
    stop("ERROR in group_by_freq: invalid seq input; Base sequence must be made up of A,C,G, and T")
  
  groupings=list()
  
  freq_array = str_count( Datin$reference_kmer, seq)
  freqmax    = length( table(freq_array) )-1
  
  for ( fr in c(0:freqmax)   ) 
  {
    groupings[[fr+1]] = Datin[ which( freq_array == fr ), ]  
  }
  
  par(mfrow=c(3,2))

  for ( fr in c(0: freqmax) )
  {
    hist( groupings[[fr+1]]$event_level_mean, 
          main   = paste(as.character(fr), "x", seq), 
          xlab   = "current mean", 
          ylab   = "", 
          xlim   = xrange, 
          breaks = break_set)
    
    # , ylim = yrange
  }
  
  
  return(groupings)
}


# -------------------------------------------------------------------------------
take_midbase_averages <-  function( Datin    = stop("Datin must be provided"), 
                                    midpoint = 3 
                                    )
{ 
RAW = list()
RAW[[1]] = Datin[ which( substr( Datin$reference_kmer, midpoint, midpoint) == "A"   ) , ]
RAW[[2]] = Datin[ which( substr( Datin$reference_kmer, midpoint, midpoint) == "C"   ) , ]
RAW[[3]] = Datin[ which( substr( Datin$reference_kmer, midpoint, midpoint) == "G"   ) , ]
RAW[[4]] = Datin[ which( substr( Datin$reference_kmer, midpoint, midpoint) == "T"   ) , ]
names(RAW) <- c("mid_A", "mid_C", "mid_G", "mid_T")

NORMED= list()
for ( i in c(1:4))
  {
  NORMED[[i]] =  ( as.numeric(RAW[[i]]$event_level_mean) - as.numeric( RAW[[i]]$model_mean) ) / as.numeric( RAW[[i]]$model_stdv ) 
  }
names(NORMED) <- c("mid_A", "mid_C", "mid_G", "mid_T")

return ( list("RAW"=RAW, "NORMED"=NORMED) )

}
