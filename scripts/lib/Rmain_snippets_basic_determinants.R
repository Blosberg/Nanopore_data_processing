olaps_4SU_shortlived_allevents = unlist( olaps_4SU$aligned_reads$shortlived )
olaps_4SU_shortlived_kmers     = split( olaps_4SU_shortlived_allevents, olaps_4SU_shortlived_allevents$model_kmer )

olaps_4SU_longlived_allevents  = unlist( olaps_4SU$aligned_reads$longlived  )
olaps_4SU_longlived_kmers      = split( olaps_4SU_longlived_allevents, olaps_4SU_longlived_allevents$model_kmer )


#===========================



weight <- unlist( lapply( seqlist, function(x) calc_seq_property( sequence_in        = x,
                                                          Base_properties_in = Base_properties,
                                                          property           = "weight")
))

wpa <- unlist( lapply( seqlist, function(x) calc_seq_property( sequence_in   = x,
                                                          Base_properties_in = Base_properties,
                                                          property           = "wpa")
))

EA_t <- unlist( lapply( seqlist, function(x) calc_seq_property( sequence_in   = x,
                                                          Base_properties_in = Base_properties,
                                                          property           = "EA_theory")
))

EA_exp <- unlist( lapply( seqlist, function(x) calc_seq_property( sequence_in   = x,
                                                          Base_properties_in = Base_properties,
                                                          property           = "EA_exp")
))

tpass_unmod <-  unlist( lapply( seqlist, function(seq) get_mean_passage_time( GR_kmer_in = unmod_longlived_kmers[[seq]] )
))
Nsample_untreated <-  unlist( lapply( seqlist, function(seq) length( unmod_longlived_kmers[[seq]] )
));
names( Nsample_untreated ) <-seqlist


tpass_4SU_shortlived <-  unlist( lapply( seqlist, function(seq) get_passage_time_statistic( GR_kmer_in = olaps_4SU_shortlived_kmers[[seq]],
                                                                                            moment = 1)
))


tpass_4SU_shortlived_sd <-  unlist( lapply( seqlist, function(seq) get_passage_time_statistic( GR_kmer_in = olaps_4SU_shortlived_kmers[[seq]],
                                                                                            moment = 2 )
))
names( tpass_4SU_shortlived_sd ) <- seqlist

# second moment defined relative to the magnitude of the mean.
moment2 = tpass_4SU_shortlived_sd/ tpass_4SU_shortlived
hist( moment2, 100, xlab ="standard_dev[time_passage]", main="kmer time passage st. dev/mean")
large_sigmaset = moment2[ moment2 > 2.5 ]
large_sigmaset_seqs = names( large_sigmaset )


# =================================
plot( weight, poremodel_RNA$mean, 
      xlab = "<mol. weight/base>", 
      ylab = "mean current", 
      main = "Sequence current dependence on mol. weight", 
      col = "blue")
abline(lm( poremodel_RNA$mean ~ weight ), lty = "dashed" , lwd = 3)

plot( EA_t, poremodel_RNA$mean, 
      xlab = "<E-affin, Theory>", 
      ylab = "mean current", 
      main = "Sequence current dependence on <e affinity>", 
      col = "blue")
abline(lm( poremodel_RNA$mean ~ EA_t ), lty = "dashed" , lwd = 3)


plot( EA_exp, poremodel_RNA$mean, 
      xlab = "<E-affin, Exp>", 
      ylab = "mean current", 
      main = "Sequence current dependence on e affin (exp)", 
      col = "blue")
abline(lm( poremodel_RNA$mean ~ EA_exp ), lty = "dashed" , lwd = 3)


plot( weight, tpass, 
      xlab = "<mol. weight/base>", 
      ylab = "mean passage time", 
      main = "mass dependence of passage time\n Untreated HEK", 
      col = "blue")
abline(lm( tpass ~ weight ), lty = "dashed" , lwd = 3)




plot( tpass_unmod, tpass_4SU_shortlived, 
      xlab = "<passage_time> HEK untreated", 
      ylab = "<passage_time> 4SU_IAA_shortlived", 
      main = "comparison of passage time by kmer", 
      col = "blue",
      log="xy")

points( tpass_unmod[ large_sigmaset_seqs ], tpass_4SU_shortlived[ large_sigmaset_seqs ], 
      xlab = "<passage_time> HEK untreated", 
      ylab = "<passage_time> 4SU_IAA_shortlived", 
      main = "comparison of passage time by kmer", 
      col = "red",
      lwd = 3)

