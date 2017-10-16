#---- here is the printed history of Vedran's "magic" from the terminal. Do this again to get actual gtf files from sam once you've got minimap output with spliced reads.

rm(list=ls()) # CLEAN UP EVERYTHING

library(GenomicAlignments)
library(data.table)


PATH_in="/Users/blosberg/postdoc_work/Nanopore/ref/"
fin_Transcript=paste0(PATH_in,"170317_Homo_sapiens.GRCh37.87.chr.gtf")
fin_readalignment=paste0(PATH_in,"/barcode05_trimmed_aligned_spliceaware_mm2.sorted.bam")

algn = readGAlignments(fin_readalignment, use.names=TRUE)
GRalgn = GRanges(algn)

# rtracklayer::export.gff(gtf, 'test.gtf')
# algn
gtf = granges(algn)
gtfl = grglist(algn)

# head(gtf)
# 
# gtf$read_id = names(gtf)
# head(gtf)
# algn
# table(table(names(algn)))
# 
# #=============================
# 
# dt = data.table(rname = names(algn))
# head(dt)
# dt$chr = as.character(seqnames(algn))
# 
# # CANT DO THIS BECAUSE OF DUPLICATE ROW NAMES (ALIGNMENTS IN MULTIPLE POSITIONS)
# # dt = as.data.table(as.data.frame(algn))
# # dt = as.data.table((algn))
# 
# # remove the names from dt:
# blgn = algn
# names(blgn) = NULL
# dt = as.data.table(as.data.frame(blgn))
# 
# # create a new string with the names from align (read name) plus the alignment position, etc. that specifies each alignment precisely.
# dt$rname = names(algn)
# head(dt)
# dt[,rid := paste(rname,seqnames,start,end,width,strand,sep='_')]
# 
# # consolodate it all into a single GRanges object blgn
# head(dt)
# names(blgn) = dt$rid
# gtf = granges(blgn)
# 
# # hold the read_id itself alone as another metadata column.
# gtf$read_id = names(gtf)
# 
# #=-------------------------------------------------------
# rtracklayer::export.gff(gtf, 'test.gff')
# 
# GRLalgn=GRangesList(GRalgn)

#=-------------------------------------------------------
# on beast:
# refFile="/data/akalin/Base/Annotation/hg19/EnsemblGenes/170317_Homo_sapiens.GRCh37.87.chr.gtf"
ref_GR    = rtracklayer::import.gff(fin_Transcript) #ref    ==> reference set of transcript ranges
ref_exons_GR               =  subset(ref_GR,    type == 'exon') 
ref_exons_by_transcript    =  split(ref_exons_GR, ref_exons_GR$transcript_id )

sample_read_cum_width      = sum(  width(gtfl) )
ref_transcript_cum_width   = sum( width(ref_exons_by_transcript ) )


# N.B. gtfl contains all possible alignments --some of which are repeated for a single read_id
olap_pairs  = get_olap_pairs_n_scores( gtfl,
                                       ref_exons_by_transcript, 
                                       sample_read_cum_width, 
                                       ref_transcript_cum_width )

olaps_scordered               = olap_pairs[ order( olap_pairs$score, decreasing=TRUE ) ]

first_occurance_bool          = !duplicated(olaps_scordered$sample_transcript_id)

best_pairs                    = olaps_scordered [first_occurance_bool ]

# select a unique entry for each read_id that represents the best possible pairing 
#--among all pairs of alignments and reference transcripts-- and compile a GRlist of such best matches

transcript_overlaps_best           = best_pairs
transcript_overlaps_best$N_readexons    = elementNROWS( gtfl[ best_pairs$queryHits ] ) 
transcript_overlaps_best$cum_width = sum(           width( gtfl[ best_pairs$queryHits ] ) )


# collect info on number of exons in each transcript:

M=20
histcount = matrix(0,M,1)
mean_score_v_Nex=matrix(0,M,1)

for ( i in c(1:M) )
  {
  histcount[i]         = length(which( transcript_overlaps_best$Nexons ==i ) )
  mean_score_v_Nex[i]  = mean( transcript_overlaps_best$score[ which ( transcript_overlaps_best$Nexons == i ) ] )
  }

#=======================================================================
# the fraction of reads that overlap with transcripts:
# take the number of unique occurances of hits for the overlap check, and divide by the number of unique reads:
reads_hitting_transcripts = length( unique(transcript_overlaps_best$queryHits) )

# total number of reads that aligned somewhere
tot_reads = length( unique(names(algn) ) )

# the fraction of successfully-aligned reads that overlapped at least partially with a transcript from the reference data
overlap_frac = reads_hitting_transcripts/tot_reads

#====================================================
# LOOK FOR "NOVEL" TRANSCRIPTS:

#-- find the reads that have no hits with existing transcripts
reads_no_reftrans_hits =  which( is.na(match( names(gtfl), transcript_overlaps_best$sample_transcript_id)) )

#-- take the subset of reads that had no hits:
gtfl_norefhit= gtfl[reads_no_reftrans_hits]

#-- reduce:
gtfl_norefhit_reduced = reduce( unlist( gtfl_norefhit) )

#-- look for histogram of overlaps:
hitcount_pairs = findOverlaps( gtfl_norefhit_reduced, gtfl_norefhit) 

#--count how many times each reduced GRange was hit:
gtfl_norefhit_reduced$hits = table( queryHits(hitcount_pairs) )

yplotmax=10000
readcov_THRESH=100
maxindex = length(gtfl_norefhit_reduced)
#-- filter only the ranges with significant coverage:
gtfl_norefhit_reduced_sighit =  gtfl_norefhit_reduced[ which(gtfl_norefhit_reduced$hits > readcov_THRESH )   ]

plot(gtfl_norefhit_reduced$hits, ylim=c(0,yplotmax), axes=FALSE, ylab="coverage", xlab="index")
# points(maxindex,yplotmax, axes=TRUE)
lines(c(1,15091),c(readcov_THRESH,readcov_THRESH), col="red", lwd=2 )

#-- take the subset with significant hits:
gtfl_NOVEL_reduced =gtfl_norefhit_reduced[ which(gtfl_norefhit_reduced$hits > readcov_THRESH ) ]

#-- now take the overlaps of the original GRanges list (without hits) against the reduced set
novel_pairs = findOverlaps( gtfl_norefhit , gtfl_NOVEL_reduced )

#-- HERE are the GRANGES of the novel transcripts with significant coverage that were not included in the reference set:
gtfl_NOVEL  = gtfl_norefhit[ queryHits(    novel_pairs ) ]

