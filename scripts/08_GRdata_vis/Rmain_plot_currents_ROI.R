# Plot current from reads along RsoI:

## Collect arguments
args <- commandArgs(TRUE)

## Default setting when no arguments passed
if(length(args) < 1) {
  args <- c("--help")
}

## Help section
if("--help" %in% args) {
  cat("
      Render to report

      Arguments:
      pathin_reads         --path to Granges List structure with reads aligned to reference genome,
      pathin_RsoI          --path to Granges list of regions of interest in the current study,
      pathout_ROIdat       --where should we send the plots

      logFile              -- self-explanatory

      Example:
      $Rscript ./Rmain_olap_reads_w_RsoI.R --arg1=1 --arg2='output.txt' --arg3=TRUE \n\n")

  q(save="no")
}

## Parse arguments (we expect the form --arg=value)
parseArgs <- function(x) strsplit(sub("^--", "", x), "=")
argsDF    <- as.data.frame(do.call("rbind", parseArgs(args)))
argsL     <- as.list(as.character(argsDF$V2))

names(argsL) <- argsDF$V1
# BiocManager::install("Biostrings")
# BiocManager::install("BSgenome.Hsapiens.UCSC.hg19")
suppressPackageStartupMessages( library(GenomicRanges) )
suppressPackageStartupMessages( library(Biostrings) )

# ======= DEBUGGING: DELETE THIS: ========
argsL=list(
"pathin_reads"   = "/scratch/AG_Akalin/bosberg/nanopore/pipeline_output/testset_0/07_GRprocessing/TESTSET0_read_ROIolap_m6A_put.rds",  # path to Granges List structure with reads aligned to reference genome
"pathin_RsoI"    = "/fast/AG_Akalin/bosberg/nanopore/ref/Regions_of_interest/m6A_putlocs_Linder.rds",                 # path to Granges list of regions of interest in the current study

"path_funcdefs"  = "/clusterhome/bosberg/projects/nanopiper/scripts/08_GRdata_vis/Rfuncs_plot_currents_ROI.R",
"pathin_refGen"  = "/fast/AG_Akalin/refGenomes/hg19_canon/hg19_canon.fa",                         # path to Granges list of regions of interest in the current study

"pathout_ROIplotdat" = "/scratch/AG_Akalin/bosberg/nanopore/pipeline_output/20180417_1233_HEK293_polyA_RNA/08_testing/TESTSET0_m6A_plotdat.rds",                             # where should we send the plots
"sampleName"     = "HEK293_untreated",
"poremodel_ref"  = "/clusterhome/bosberg/projects/nanopiper/dev/ref/poremodel_RNA.csv",

"mincov_in"      = 10,
"plotrange_in"   = 10,

"logFile"        = "/scratch/AG_Akalin/bosberg/nanopore/pipeline_output/20180417_1233_HEK293_polyA_RNA/08_testing/test.log"                             # self-explanatory
)

source( argsL$path_funcdefs )
# ======= DEBUGGING: DOWN TO HERE: ========

# overlaps CONVENTION: findoverlaps ( region of interest, reads        )
#                                      < queryHits >    | < subjectHits >
#                                     __________________|______________
#                                       locus_i         |   read_i
#                                       ...             |   ...

# ========================================
# Read in inputs:
ref_Genome    <- readDNAStringSet( argsL$pathin_refGen )
readdat       <- readRDS( argsL$pathin_reads ) # these are just all the reads with at least one overlap on a ROI:
                                               # no structural organization in place (as some reads might overlap multiple RsOI)
putloci       <- readRDS( argsL$pathin_RsoI  )
poremodel     <- read.table( file = argsL$poremodel_ref,
                             sep = "\t",
                             header = TRUE
                             )
# ========================================
# Read in inputs:

OLAP_skip_TOL = 3

# writeLines( "Filtering loci for coverage.",
#             argsL$logFile )
#
# reduce list of RsoI to those with sufficient coverage.
loci_filtered_for_coverage <- filter_loci_for_coverage (  loci   = putloci$Region_groups,
                                                          reads  = readdat$aligned_reads,
                                                          mincov = argsL$mincov_in,
                                                          OLAP_skip_TOL = OLAP_skip_TOL
                                                        )

# expand loci slightly to ensure overlap is registered,
# even if the exact base is skipped.
loci_filtered_for_coverage_expanded <- loci_filtered_for_coverage
for ( group in  names( loci_filtered_for_coverage )  )
  {
  if ( identical( all( width ( loci_filtered_for_coverage_expanded[[group]] )  < OLAP_skip_TOL ) , TRUE ) )
    {
    start( loci_filtered_for_coverage_expanded[[group]] ) <-  ( start( loci_filtered_for_coverage[[group]] ) -OLAP_skip_TOL )
    end(   loci_filtered_for_coverage_expanded[[group]] ) <-  ( end(   loci_filtered_for_coverage[[group]] ) +OLAP_skip_TOL )
    }
  }

# NOW check overlaps of reads with only the _covered_ loci.
#TODO: We run this overlaps more than once; if speedup is required later, eliminate this repetition.
overlaps_by_group = lapply ( names( loci_filtered_for_coverage ), function(group)
                               findOverlaps(  loci_filtered_for_coverage_expanded[[group]],
                               readdat$aligned_reads[[group]] )
                               )
names( overlaps_by_group ) <- names( loci_filtered_for_coverage )

# Don't need this anymore. The loci should refer directly to the exact ROI
rm(loci_filtered_for_coverage_expanded)

# overlaps_by_group queryHits now references the indices of COVERED loci.

# TODO: add row/col names for the output data structures.
sampleROI_dat_by_group <- lapply( names(loci_filtered_for_coverage),
                                      function(group) collect_group_sample_dat_over_ROIs(
                                                          SampleName           = argsL$sampleName,
                                                          ref_Genome           = ref_Genome,
                                                          aligned_group_reads  = readdat$aligned_reads[[group]],
                                                          group_overlaps       = overlaps_by_group[[group]],
                                                          group_loci_covered   = loci_filtered_for_coverage[[group]],
                                                          poremodel_in         = poremodel,
                                                          plotrange_in         = plotrange_in )
                                    )
names(sampleROI_dat_by_group) <- names(loci_filtered_for_coverage)

saveRDS( object = sampleROI_dat_by_group,
         file   = pathout_ROIplotdat      )

# writeLines( "Finished exporting RDS file for plot data.",
#             argsL$logFile )
#
# Everything below here is buffer text to use as templates for function execution:
# group="CITS"
# group="Ill_2pO"
# i=1
plot_samplesignal_over_ROI( sampleROI_dat  = sampleROI_dat_by_group[[group]][[i]],
                            refgen         = ref_Genome,
                            squiggle_type  = "event",
                            line_darkness  = 0.1,
                            normed         =  F )
# =======================================
plot_dwelltime_over_ROI( sampleROI_dat = sampleROI_dat_by_group[[group]][[i]],
                         refgen        = ref_Genome,
                         plot_logdwell = T )

# =======================================

i=0;
i=i+1; pca_normdiff <- get_pca_normdiff (  sampleROI_dat_in =  sampleROI_dat_by_group[[group]][[i]]  )


group="CITS"
pca_clust_dat=list()
Ngroups = length( names( sampleROI_dat_by_group ) )

pca_clust_dat= lapply( c( 1: Ngroups ), function(group)
                       lapply( c(1:length(sampleROI_dat_by_group[[group]] ) ),
                              function(x) get_pca_clusters (  sampleROI_dat_in =  sampleROI_dat_by_group[[group]][[x]],
                                                               shouldplot      = FALSE) )  )
names( pca_clust_dat ) <- names( sampleROI_dat_by_group )


c2_means   = list()
Sil_scores = list()

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
