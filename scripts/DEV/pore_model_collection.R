#--- example script of how I got the nanopore model data.
#--- here reads_GR is just a GRanges object of a large data set that I collected
#--- (i.e. large enough to be sure that every 5-mer occured at least once)

seqlist=which( !duplicated( reads_GR$reference_kmer) ) 
undupped = reads_GR[seqlist]

seqs_sorted = sort( reads_GR$reference_kmer[seqlist] )
indices_sorted <- match( seqs_sorted , undupped$reference_kmer)
undupped_sorted <- undupped[ indices_sorted ]


pore_model <- cbind( as.numeric(undupped_sorted$model_mean), 
                     as.numeric(undupped_sorted$model_stdv) 
)
row.names(pore_model) <- undupped_sorted$reference_kmer
colnames(pore_model) <- c("mean", "std_dev")
write.csv(pore_model, file = "~/projects/nanopore/scripts/ref/pore_model_table.csv")
temp = read.csv(file = "~/projects/nanopore/scripts/ref/pore_model_table.csv", row.names = 1)

