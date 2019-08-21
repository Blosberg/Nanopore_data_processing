# Everything below here is buffer text to use as templates for function execution:
# group="CITS"
# group="Ill_2pO"
# i=1

suppressPackageStartupMessages( library(GenomicRanges) )
suppressPackageStartupMessages( library(Biostrings) )
suppressPackageStartupMessages( library(cluster) )
source('~/projects/nanopiper/scripts/08_GRdata_vis/Rfuncs_plot_squiggles_pca.R')

mincov = 10
SampleName = "HEK293_untreated"
Path_refgen = "/scratch/AG_Akalin/refGenomes/hg19_canon/hg19_canon.fa"
ref_Genome    <- readDNAStringSet( Path_refgen )

reads_PATH     = "/fast/AG_Akalin/bosberg/nanopore/pipeline_output/20180417_1233_HEK293_polyA_RNA/08_plotdat_vis/"

## ----------------------------------------------------------------------------
PATH_CITS    = paste0( reads_PATH, "HEK_untreated-olap-m6a_CITS-plotdat.rds")
PATH_Me2p0   = paste0( reads_PATH, "HEK_untreated-olap-Me2p0-plotdat.rds" )

plotdat_Me2p0 = readRDS( PATH_Me2p0 )

Me2p0_plotdat_covered <- filter_for_coverage( plotdat_Me2p0$ROIdat, 
                                              mincov = mincov )

CITS_clustdat_covered <- lapply( c(1:length(CITS_plotdat_covered) ), function(locus)
                                    iterate_pca_clusters ( sampleROI_dat_in = CITS_plotdat_covered[[locus]],
                                                           shouldplot       = FALSE,
                                                           k_in             = 4,
                                                           mincov           = mincov 
                                                          )
                                )
CITS_metaplotdat <- gather_metaplot_clusterdat (  CITS_clustdat_covered )

## ----------------------------------------------------------------------------
PATH_CITS_NC50up    = paste0( CITS_DIR ,"HEK_untreated-olap-m6a_CITS_NC_50up-plotdat.rds" )
plotdat_CITS_NC50up = readRDS( PATH_CITS_NC50up )
CITS_plotdat_covered_NC50up <- filter_for_coverage( plotdat_CITS_NC50up$ROIdat, 
                                                    mincov = mincov )

CITS_clustdat_covered_NC50up <- lapply( c( 1: length(CITS_plotdat_covered_NC50up) ), function(locus)
                                  iterate_pca_clusters ( sampleROI_dat_in = CITS_plotdat_covered_NC50up[[locus]],
                                                         shouldplot       = FALSE,
                                                         k=4,
                                                         mincov = mincov )
                       
                          )
CITS_metaplotdat_NC50up <- gather_metaplot_clusterdat (  CITS_clustdat_covered_NC50up )

## ----------------------------------------------------------------------------

PATH_CITS_NC50dn    = paste0( CITS_DIR ,"HEK_untreated-olap-m6a_CITS_NC_50dn-plotdat.rds" )
plotdat_CITS_NC50dn = readRDS( PATH_CITS_NC50dn )
CITS_plotdat_covered_NC50dn <- filter_for_coverage( plotdat_CITS_NC50dn$ROIdat, 
                                                    mincov = mincov )

CITS_clustdat_covered_NC50dn <- lapply( c( 1: length(CITS_plotdat_covered_NC50dn) ), function(locus)
                                  iterate_pca_clusters ( sampleROI_dat_in = CITS_plotdat_covered_NC50dn[[locus]],
                                                         shouldplot       = FALSE,
                                                         k=4,
                                                         mincov = mincov )
                       
                          )
CITS_metaplotdat_NC50dn <- gather_metaplot_clusterdat (  CITS_clustdat_covered_NC50dn )

## ----------------------------------------------------------------------------



# which c2 > 10:
# 39  76  77 103 147 168 190 202 213
p=76 
# =======================================
plot_samplesignal_over_ROI( SampleName     = SampleName, 
                            sampleROI_dat  = CITS_plotdat_covered[[p]],
                            refgen         = ref_Genome,
                            squiggle_type  = "none",
                            line_darkness  = 0.5,
                            normed         = F   )
# =======================================
plot_dwelltime_over_ROI( SampleName    = SampleName,
                         sampleROI_dat = CITS_plotdat_covered[[p]],
                         refgen        = ref_Genome,
                         plot_logdwell = T )

# =======================================
pca_normdiff <- get_pca_clusters ( sampleROI_dat_in  = CITS_plotdat_covered[[p]],
                                   should_zeropadd   = TRUE,
                                   k=4
                                   )
# =======================================
locusdat <- iterate_pca_clusters ( sampleROI_dat_in = CITS_plotdat_NC_50dn_covered[[p]],
                                   shouldplot       = TRUE,
                                   k_in = 4,
                                   mincov = mincov
                                  )

# =======================================



CITS_metaplotdat_NC50up <- gather_metaplot_clusterdat ( clusterdat_in = CITS_clustdat_covered_NC50up )

CITS_clustdat_covered_50dn = lapply( c( 1: length(CITS_plotdat_NC_50dn_covered) ), function(locus)
                                       iterate_pca_clusters ( sampleROI_dat_in = CITS_plotdat_NC_50dn_covered[[locus]],
                                                              shouldplot       = FALSE,
                                                              k                = 4,
                                                              mincov           = mincov
                                                            )
                                    )
CITS_metaplotdat_50dn <- gather_metaplot_clusterdat ( clusterdat_in = CITS_clustdat_covered_50dn )





                       



# extract important observables:
c2_means <- lapply( c( 1: Ngroups ), function(group)
                      unlist( lapply( c(1: length(sampleROI_dat_by_group[[group]] ) ), function(loci)
                             rowSums( pca_clust_dat[[group]][[loci]]$clust_kmeans$centers^2 )[2]  ) )
                   )
names( c2_means ) <- names( sampleROI_dat_by_group )

Sil_scores <- lapply( c( 1: Ngroups ), function(group)
                     unlist( lapply( c(1: length(sampleROI_dat_by_group[[group]] ) ), function(loci)
                              mean( pca_clust_dat[[group]][[loci]]$Silh[,3] )  ) )
                     )
names( Sil_scores ) <- names( sampleROI_dat_by_group )


pca_dist <- dist( pca_normdiff$pca$x )
readraw_dist <- dist( sampleROI_dat_by_group[[group]][[i]]$read_normdiff )


rowSums( pca_normdiff$pca$x * pca_normdiff$pca$x )
