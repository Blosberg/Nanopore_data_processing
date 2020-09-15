
# MMe2p0
PATH_putmod_NC50dn    = paste0( sample_plotdat_DIR ,"HEK_untreated-olap-Me2p0_NC50dn-plotdat.rds" )
# CITS:
PATH_putmod_NC50dn    = paste0( sample_plotdat_DIR ,"HEK_untreated-olap-m6a_CITS_NC_50dn-plotdat.rds" )

plotdat_putmod_NC50dn = readRDS( PATH_putmod_NC50dn )
putmod_plotdat_covered_NC50dn <- filter_for_coverage( plotdat_putmod_NC50dn$ROIdat, 
                                                    mincov = mincov )

putmod_clustdat_covered_NC50dn <- lapply( c( 1: length(putmod_plotdat_covered_NC50dn) ), function(locus)
                                  iterate_pca_clusters ( sampleROI_dat_in = putmod_plotdat_covered_NC50dn[[locus]],
                                                         shouldplot       = FALSE,
                                                         k_in = 5,
                                                         mincov = mincov )
                       
                          )
putmod_metaplotdat_NC50dn <- gather_metaplot_clusterdat (  putmod_clustdat_covered_NC50dn )

## ----------------------------------------------------------------------------

 plot( unlist( putmod_metaplotdat_NC50dn$c2_means), 
       unlist( putmod_metaplotdat_NC50dn$Silh ) )
## ----------------------------------------------------------------------------