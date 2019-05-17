histlist_4SU_IAA         = readRDS( "/scratch/AG_Akalin/bosberg/nanopore/pipeline_output/20190130_1420_293_4Su_IAA/06_GRobjects/6c2fe901_kmer_histlist.rds")
histlist_HEK_unmod       = readRDS("/scratch/AG_Akalin/bosberg/nanopore/pipeline_output/20180417_1233_HEK293_polyA_RNA/06_GRobjects/HEK293_polyA_kmer_histlist.rds")

reads_HEK_unmod = readRDS("/scratch/AG_Akalin/bosberg/nanopore/pipeline_output/20180417_1233_HEK293_polyA_RNA/06_GRobjects/HEK293_polyA_reads_GRL.rds")
reads_4SU_IAA   = readRDS("/scratch/AG_Akalin/bosberg/nanopore/pipeline_output/20190130_1420_293_4Su_IAA/06_GRobjects/6c2fe901_reads_GRL.rds")
  

reads_flattened_HEK_unmod = lapply( reads_HEK_unmod$Events_GRL_splitbyread,         function(r) flatten_read( read_GR_in = r) )
reads_flattened_4SU_IAA   = lapply( reads_flattened_4SU_IAA$Events_GRL_splitbyread, function(r) flatten_read( read_GR_in = r) )

                                    
                                    
kmer_olap_4SU_WT     <- unlist( lapply( seqs, function(seq)  calculate_histogram_overlap (  
                                            hist1    = histlist_4SU_IAA[[seq]],
                                            hist2    = histlist_HEK_unmod[[seq]],
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
                         hist2    = histlist_HEK_unmod,
                         seq      = "GGGAA" )


plot_single_histogram  ( hist     = histlist_4SU_IAA  ,
                         seq      = "GGTAC" )


means_4SU     <- unlist( lapply( seqs, function(seq)  hist_mean (   hist_in    = histlist_4SU_IAA[[seq]],
                                                                    seq_in      = seq
                                                                )
                                           ) )
                                        
means_WT      <- unlist( lapply( seqs, function(seq)  hist_mean (   hist_in    = histlist_HEK_unmod[[seq]],
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
    seq_selfoverlap [seq1, seq2]  <- calculate_histogram_overlap(  hist1    = histlist_HEK_unmod[[seq1]],
                                                                   hist2    = histlist_HEK_unmod[[seq2]]
                                                                 )
      
      } 
  }
  

trace_histogram( histlist_in  = stop("reads must be provided"),
                 sequence_in  = stop("sequence must be provided")
)

plot_background_hist( xdat     = xdat,
                      hist_in  = trace_histogram( histlist_in  = histlist_4SU_IAA,
                                                  sequence_in  = "AANAA") 
                     )
# ================================
# now look at the middle trimers:

middle_trimers_tracenames <- sort( paste0("N", sequence_trace("NNN") , "N") )
names(middle_trimers_tracenames) <- middle_trimers_tracenames

middle_trimers_hists  <-  lapply( middle_trimers_tracenames, function(seq)  trace_histogram ( histlist_in  = histlist_HEK_unmod,
                                                                                              sequence_in  = seq)
                                 )
names(middle_trimers_hists) <- middle_trimers_tracenames
                                 
nseqs = length(middle_trimers_hists)
trimer_selfoverlap = matrix( 0, nseqs, nseqs) 

colnames(trimer_selfoverlap) <- middle_trimers_tracenames
rownames(trimer_selfoverlap) <- middle_trimers_tracenames

for ( seq1 in middle_trimers_tracenames )
{ 
  for (seq2 in middle_trimers_tracenames )
  {
    trimer_selfoverlap [seq1, seq2]  <-  sum( rowMins(  cbind(middle_trimers_hists[[seq1]],
                                                              middle_trimers_hists[[seq2]])  
                                                      )
                                            )
      } 
  }
  
# --- check monotonicity: (i.e. how well this trace corresponds to its components)

subtrace_overlap = list()
  
for ( seq1 in middle_trimers_tracenames )
  { 
  seqs <- sequence_trace( sequence_in = seq1)

  tracehist = trace_histogram ( histlist_in = histlist_HEK_unmod,
                                sequence_in = seq1,
                                return_as_histlist_element = TRUE)
  
  subtrace_overlap = lapply( seqs, function(x) 
                             calculate_histogram_overlap (  hist1    = tracehist,
                                                            hist2    = x
                                                          )
                             )
  }
  

lapply( c(1:length(temp2)), function(x)
                               ( (mean(temp2[[x]]$event_mean) -
                                    mean(GRcontrolL_seqsplit[[ names(temp2)[x] ]]$event_mean))
     #                            /(pore_model_list[[ names(temp2)[x]  ]]$std_dev)  )
                                )
                              )