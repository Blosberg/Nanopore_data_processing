
histlist_4SU_IAA = readRDS( "/scratch/AG_Akalin/bosberg/nanopore/pipeline_output/20190130_1420_293_4Su_IAA/06_GRobjects/6c2fe901_kmer_histlist.rds")
histlist_HEK_WT  = readRDS("/scratch/AG_Akalin/bosberg/nanopore/pipeline_output/20180417_1233_HEK293_polyA_RNA/06_GRobjects/HEK293_polyA_kmer_histlist.rds")



kmer_olap_4SU_WT     <- unlist( lapply( seqs, function(seq)  calculate_histogram_overlap (  
                                            hist1    = histlist_4SU_IAA[[seq]],
                                            hist2    = histlist_HEK_WT[[seq]],
                                            plot     = FALSE,
                                            col_in1  = rgb( 1,      0,     0,   alpha=0.5 ),
                                            col_in2  = rgb( 0,      0,     1,   alpha=0.5 )
                                            )
                  ) )
names(kmer_olap_4SU_WT ) <- seqs

# gather the distributions with olap below threshold.
diff_observed <- kmer_olap_4SU_WT[ kmer_olap_4SU_WT < 0.85 ] 

#take the central 3 bases:
centres <- substr( names( diff_observed ), 2,4 )
blob    <- paste( centres, collapse = '')


hist( kmer_olap_4SU_WT, 100, 
      main = "Current profile overlap: 4SU_IAA vs. WT",
      xlab = "overlap [0:1]"
      )



plot_histogram_overlap ( hist1    = histlist_4SU_IAA  ,
                         hist2    = histlist_HEK_WT,
                         seq      = "GGGAA" )


plot_single_histogram  ( hist     = histlist_4SU_IAA  ,
                         seq      = "GGTAC" )


hist_mean


means_4SU     <- unlist( lapply( seqs, function(seq)  hist_mean (   hist_in    = histlist_4SU_IAA[[seq]],
                                                                    seq_in      = seq
                                                                )
                                           ) )
                                        
means_WT      <- unlist( lapply( seqs, function(seq)  hist_mean (   hist_in    = histlist_HEK_WT[[seq]],
                                                                    seq_in      = seq
                                                                )
                                           ) )
                                        

seq_selfoverlap = matrix(0, length(seqs) , length(seqs) )
colnames(seq_selfoverlap)  <- seqs
row.names(seq_selfoverlap) <-seqs

seq_selfoverlap[1:4,1:6]

for ( seq1 in seqs )
{ 
  for (seq2 in seqs )
  {
    seq_selfoverlap [seq1, seq2]  <- calculate_histogram_overlap(  hist1    = histlist_HEK_WT[[seq1]],
                                                                   hist2    = histlist_HEK_WT[[seq2]]
                                                                 )
      
      } 
  }
  

trace_histogram <- function( histlist_in  = stop("reads must be provided"),
                             sequence_in  = stop("sequence must be provided")

plot_background_hist( xdat     = xdat,
                      hist_in  = trace_histogram( histlist_in  = histlist_4SU_IAA,
                                                  sequence_in  = "AANAA") 
                     )
  