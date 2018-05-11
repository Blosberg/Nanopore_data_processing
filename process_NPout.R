sourcefile="/Users/blosberg/postdoc_work/projects/Nanopore_upstream/output/reads-ref_0.eventalign.txt"

Nano0=read.csv(sourcefile, sep="\t", header = TRUE )


first_occurance_index         = which(!duplicated(Nano0$read_index))

N_alignments = length( first_occurance_index )

last_index= c( (first_occurance_index-1)[2:N_alignments], dim(Nano0)[1] )

Reads_GR = GRanges( seqnames = Nano0$contig[first_occurance_index], 
                    ranges   = IRanges( start=Nano0$position[first_occurance_index], 
                                     end=Nano0$position[last_index] ))


ovlps = findOverlaps( Reads_GR, Reads_GR)
ovlps = ovlps[ which( !( queryHits(ovlps)== subjectHits(ovlps) )  ) ]
