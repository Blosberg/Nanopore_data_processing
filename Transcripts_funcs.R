# get_exon_list <-  function( Target_in, DB_in) 
#   {
#   GenomicRanges::intersect(  Target_in, DB_in,
#                              ignore.strand=TRUE) ) ) #-- total length shared (or "common") to both vecs
# 
#   
#   }

check_files <-  function( InFile, RefFile) 
{ 
  if( !file.exists(InFile))
  {  stop(paste("ERROR: supplied input file: ",  InFile, " does not exist. Exiting") )  }
  if( !file.exists(RefFile))
  {  stop(paste("ERROR: supplied reference file: ",  RefFile, " does not exist. Exiting") )  }
}

#=========================================================================================
get_single_overlap <-  function( t1, t2) 
{ 
  LA = sum(width(t1))
  LB = sum(width(t2))
  
  LC = 2*sum(
    width( 
      reduce(
        intersect(  t1,
                    t2,
                    ignore.strand=TRUE
        )
      )
    )
  )
  score = LC/(LA + LB)      
  
  return(score)
}
#=========================================================================================
# break this up into pairs:
get_olap_pairs_n_scores <-  function( sample_exons_by_transcript, ref_exons_by_transcript, sample_transcript_cum_width, ref_transcript_cum_width ) 
{
  transcript_overlaps = as.data.table( 
    findOverlaps( sample_exons_by_transcript,ref_exons_by_transcript) 
  )
  #=== create an array of these lengths (some of which repeat) for every pair sample/reference overlaps. 
  LA = sample_transcript_cum_width[transcript_overlaps$queryHits]  
  LB = ref_transcript_cum_width[transcript_overlaps$subjectHits]
  
  #=== now create an array (which does not repeat) of the intersection lengths for each pair
  #=== note the factor of 2, since these lengths are found on _both_ transcripts.
  LC = 2*sum(
    width( 
      reduce(
        intersect(   sample_exons_by_transcript[transcript_overlaps$queryHits], 
                     ref_exons_by_transcript[transcript_overlaps$subjectHits],
                     ignore.strand=TRUE
        )
      )
    )
  )
  
  # compute the score for each pair set
  transcript_overlaps$score = LC/(LA+LB)
  
  # add transcript_id labels to the pair sets
  transcript_overlaps$sample_transcript_id = names( sample_exons_by_transcript)[transcript_overlaps$queryHits]
  transcript_overlaps$ref_transcript_id    = names( ref_exons_by_transcript)[transcript_overlaps$subjectHits]
  
  # sanity check:
  if(min(transcript_overlaps$score) < 0 || max(transcript_overlaps$score) >1 )
  { stop("ERROR: calculated overlap scores are not bounded by [0,1].") }
  
  return( transcript_overlaps )
}

#=========================================================================================
select_best_olap_score <- function (transcript_overlaps)
{
  #--- now order the results by the intersection score (which we calculated)
  #--- "scordered" ==> "ordered by score"
  transcript_overlaps_scordered = transcript_overlaps[order(transcript_overlaps$score , decreasing=TRUE)]
  
  
  #--- Find the first occurance of each transcript_id by decreasing score 
  first_occurance_index         = !duplicated(transcript_overlaps_scordered$sample_transcript_id)
  
  #--- Take this highest-scoring (best) result for each sample_transcript_id
  transcript_overlaps_bestscore = transcript_overlaps_scordered[first_occurance_index]
  
  print("take_best is FALSE in get_overlaps; returning all overlaps (multiple copies per sample value)")
  return(transcript_overlaps_bestscore)
}

#=====================================================================================
#----- in case you want to run the above with take_best==FALSE, you can then just do this filter after:
filter_best_olap <- function(transcript_overlaps)
{
  #--- now order the results by the intersection score (which we calculated)
  #--- "scordered" ==> "ordered by score"
  transcript_overlaps_scordered = transcript_overlaps[order(transcript_overlaps$score , decreasing=TRUE)]
  
  
  #--- Find the first occurance of each transcript_id by decreasing score 
  first_occurance_index         = !duplicated(transcript_overlaps_scordered$sample_transcript_id)
  
  #--- Take this highest-scoring (best) result for each sample_transcript_id
  transcript_overlaps_bestscore = transcript_overlaps_scordered[first_occurance_index]
  
  print("take_best is FALSE in get_overlaps; returning all overlaps (multiple copies per sample value)")
  return(transcript_overlaps_bestscore)
}

#=====================================================================================

plot_transcript_scores <- function( transcripts_GR, scorethresh=0.90, alpha_weight=0.05, alpha_weight_bar=0.9, cov_thresh=1e+1, mincov=1.5, maxcov=1E4)
{
  thickness    = 4
  # maxcov       = max(transcripts_GR$gene_cov)
  # mincov       = min(transcripts_GR$gene_cov)
  
  Ntranscripts = length(sample_transcripts_GR) 
  
  plot( sample_transcripts_GR$gene_cov[1:Ntranscripts], sample_transcripts_GR$score[1:Ntranscripts], log="x", xlim=c( mincov, maxcov), ylim=c(0,1), xlab="gene coverage", ylab = "max olap score",
        pch=20, col=rgb(0,0,0,alpha=alpha_weight))
  lines( c(mincov, maxcov), c(scorethresh, scorethresh) ,col=rgb(1,0,0,alpha=alpha_weight_bar), lwd=thickness)
  lines( c(cov_thresh, cov_thresh), c(0, scorethresh) , col=rgb(0,0,1,alpha=alpha_weight_bar), lwd=thickness)
  
  #  hist(temp$intersect_score[1:L],100, xlab="best overlap score found", main=paste("first ", as.character(L), " ranked by gene coverage")  )
  #---- this is only if you just want to order each transcript on its own by coverage
  # overlapping_trans_sel = overlapping_trans_sel[order(overlapping_trans_sel$intersect_score, 
  #                                                    decreasing=TRUE)]
  
  # overlaping_trans_intersect = sum(width(pintersect(strT_dat_exons, ref_dat_exons)))
  # save.image(file=paste0("workspaces/",dat_StrT_in_str,"_out.RData"))
  
  
}