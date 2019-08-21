
# MMe2p0
# PATH_putmod_NC50up    = paste0( sample_plotdat_DIR ,"HEK_untreated-olap-Me2p0_NC50up-plotdat.rds" )
# CITS:
PATH_putmod_NC50up    = paste0( sample_plotdat_DIR ,"HEK_untreated-olap-m6a_CITS_NC_50up-plotdat.rds" )

plotdat_putmod_NC50up = readRDS( PATH_putmod_NC50up ) # length: 1584
putmod_plotdat_covered_NC50up <- filter_for_coverage( plotdat_putmod_NC50up$ROIdat, 
                                                      mincov = mincov ) # length: 974

putmod_clustdat_covered_NC50up <- lapply( c( 1: length(putmod_plotdat_covered_NC50up) ), function(locus)
                                  iterate_pca_clusters ( sampleROI_dat_in = putmod_plotdat_covered_NC50up[[locus]],
                                                         shouldplot       = FALSE,
                                                         k_in = 5,
                                                         mincov = mincov )
                          ) # length: 974

putmod_metaplotdat_NC50up <- gather_metaplot_clusterdat (  putmod_clustdat_covered_NC50up )

# -------------------------------------------------
 plot( unlist( putmod_metaplotdat_NC50up$c2_means), 
       unlist( putmod_metaplotdat_NC50up$Silh ) )

# outlier_set = which( putmod_metaplotdat_NC50up$c2_means > 5)
# mapping <- which(nchar(putmod_metaplotdat_NC50up$Flags) ==0 )[outlier_set]
# plot_pca_clusters ( pca     = putmod_clustdat_covered_NC50up[[p]]$pca,
#                     cluster = putmod_clustdat_covered_NC50up[[p]]$clusterdat$cluster
#                    )
# 
# plot_samplesignal_over_ROI( SampleName     = "MM_2pO", 
#                             sampleROI_dat  = putmod_plotdat_covered_NC50up[[p]],
#                             refgen         = ref_Genome,
#                             squiggle_type  = "none",
#                             line_darkness  = 0.5,
#                             normed         = F   )

## ----------------------------------------------------------------------------