# ============================================================
# NOW LOOK FOR MODIFICATIONS:

mid_groups = take_midbase_averages( Datin = dat_win_finite )

subset_allA = dat_win_finite[ which( dat_win_finite$reference_kmer == "AAAAA") ,]
subset_allC = dat_win_finite[ which( dat_win_finite$reference_kmer == "CCCCC") ,]
subset_allG = dat_win_finite[ which( dat_win_finite$reference_kmer == "GGGGG") ,]
subset_allT = dat_win_finite[ which( dat_win_finite$reference_kmer == "TTTTT") ,]
par( mfrow = c(2,2) )
hist(subset_allA$event_level_mean, breaks = seq(from = mincurrent, to = maxcurrent, by = 1 ) , xlab = "current", main = "AAAAA")
hist(subset_allC$event_level_mean, breaks = seq(from = mincurrent, to = maxcurrent, by = 1 ) , xlab = "current", main = "CCCCC")
hist(subset_allG$event_level_mean, breaks = seq(from = mincurrent, to = maxcurrent, by = 1 ) , xlab = "current", main = "GGGGG")
hist(subset_allT$event_level_mean, breaks = seq(from = mincurrent, to = maxcurrent, by = 1 ) , xlab = "current", main = "TTTTT")


freq_group_A = group_by_freq ( Datin           = dat_win_finite, 
                               seq             = "A",
                               xrange          = current_window,
                               break_set       = seq(from = mincurrent, to = maxcurrent, by = 1 ),
                               plot_panel_dim  = c(3,2)
)


m6A_motif="GGACT"
m6A_motif_subset  = dat_win_finite[ which( dat_win_finite$reference_kmer == m6A_motif),  ]
m6A_SNP_subset_C = dat_win_finite[ which( dat_win_finite$reference_kmer == "GGCCT") ,  ]
m6A_SNP_subset_G = dat_win_finite[ which( dat_win_finite$reference_kmer == "GGGCT") ,  ]
m6A_SNP_subset_T = dat_win_finite[ which( dat_win_finite$reference_kmer == "GGTCT") ,  ]
par( mfrow = c(2,2) )

hist(m6A_motif_subset$event_level_mean, 
     breaks = seq(from = mincurrent, to = maxcurrent, by = 1 ) , 
     xlab = "current", 
     main = m6A_motif )

hist(m6A_SNP_subset_C$event_level_mean, 
     breaks = seq(from = mincurrent, to = maxcurrent, by = 1 ) , 
     xlab = "current", 
     main = "CCCGA")

hist(m6A_SNP_subset_G$event_level_mean, 
     breaks = seq(from = mincurrent, to = maxcurrent, by = 1 ) , 
     xlab = "current", 
     main = "CCGGA")

hist(m6A_SNP_subset_T$event_level_mean, 
     breaks = seq(from = mincurrent, to = maxcurrent, by = 1 ) , 
     xlab   = "current", 
     main   = "CCTGA")


#==================================================================
# ------  SAVED R-WORKSPACE AFTER EXECUTING DOWN TO HERE ----
load("/home/bosberg/projects/nanopore/signal_processing_Rworkspace.RData") 

k=5
dat_win_finite_GR = GRanges( seqnames = dat_win_finite$contig,
                             ranges      = IRanges (start = dat_win_finite$position,
                                                    end   = (dat_win_finite$position+k) 
                             ),
                             read_index        = dat_win_finite$read_index,
                             reference_kmer    = dat_win_finite$reference_kmer,
                             event_level_mean  = dat_win_finite$event_level_mean,
                             event_stdv        = dat_win_finite$event_stdv,
                             model_kmer        = dat_win_finite$model_kmer,
                             event_index       = dat_win_finite$event_index
)

save.image(file="/home/bosberg/projects/nanopore/signal_processing_Rworkspace.RData")


i_asc       = which( dat_win_finite_GR$reference_kmer == dat_win_finite_GR$model_kmer )

i_ref_eq_mod = which( dat_win_finite_GR$reference_kmer == dat_win_finite_GR$model_kmer )
dat_ref_eq_mod_GR  =  dat_win_finite_GR [ i_ref_eq_mod ]

i_ref_neq_mod = which( dat_win_finite_GR$reference_kmer != dat_win_finite_GR$model_kmer )
dat_ref_neq_mod_GR = dat_win_finite_GR [ i_ref_neq_mod ]

reads_GRL = split( dat_win_finite_GR, dat_win_finite_GR$read_index )


# ===================================================
# -- Plotting:

r <- hist(  m6A_motif_subset$event_level_mean )
r$counts= r$counts+1

par( mfrow = c(1,1) )
hist( m6A_motif_subset$event_level_mean, 
      main        ="m6A subset", 
      xlab        ="current mean", 
      xlim        = current_window, 
      ylab        = "", 
      breaks      = seq(from = mincurrent, to = maxcurrent, by = 1 ) )



# ========================================================================
# given input position (1 through 5), separate data based on which base is in that position.
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
# ========================================================================
# Group by how many times a sequence appears within a 5-bp window
group_by_freq <-  function( Datin = stop("Datin must be provided"), 
                            sequence_length = 5, 
                            seq       = stop("target sequence to look for must be provided"),
                            xrange    = c(50,150),
                            break_set = seq(from = 49, to = 151, by = 2),
                            plot_panel_dim = c(3,2)
)
{ 
  if ( str_count(seq,"A") + str_count(seq,"C") + str_count(seq,"G") + str_count(seq,"T")  != nchar(seq) )
    stop("ERROR in group_by_freq: invalid seq input; Base sequence must be made up of A,C,G, and T")
  
  groupings=list()
  
  #COLLECT AN ARRAY SHOWING THE FREQUENCY WITH WHICH THE TARGET SEQUENCE APPEARS
  freq_array = str_count( Datin$reference_kmer, seq)
  freqmax    = length( table(freq_array) )-1
  
  # BUILD GROUPINGS AROUNG THESE FREQUENCIES.
  for ( fr in c(0:freqmax)   ) 
  {
    groupings[[fr+1]] = Datin[ which( freq_array == fr ), ]  
  }
  
  par( mfrow = plot_panel_dim )
  
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


# ========================================================================
# Partition based on which base is in the middle.
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

