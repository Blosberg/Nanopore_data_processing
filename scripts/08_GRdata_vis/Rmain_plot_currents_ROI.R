# Plot current from reads along RsoI:

BLAMAROONY

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
      pathout_plot         --where should we send the plots
      assembly             --genome
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
BiocManager::install("Biostrings")
BiocManager::install("BSgenome.Hsapiens.UCSC.hg19")
suppressPackageStartupMessages( library(GenomicRanges) )
suppressPackageStartupMessages( library(Biostrings) )

##==============================================================


# ======= DEBUGGING: DELETE THIS: ========
argsL=list(
"pathin_reads"   = "/scratch/AG_Akalin/bosberg/nanopore/pipeline_output/20180417_1233_HEK293_polyA_RNA/07_GRprocessing/HEK_untreated_read_ROIolap_m6A_put.rds",  # path to Granges List structure with reads aligned to reference genome
"pathin_RsoI"    = "/fast/AG_Akalin/bosberg/nanopore/ref/Regions_of_interest/m6A_putlocs_Linder.rds",                 # path to Granges list of regions of interest in the current study
"pathin_refGen"  = "/fast/AG_Akalin/refGenomes/hg19_canon/hg19_canon.fa",                         # path to Granges list of regions of interest in the current study
"pathout_plot"   = "/scratch/AG_Akalin/bosberg/nanopore/pipeline_output/20180417_1233_HEK293_polyA_RNA/08_testing/testplot",                             # where should we send the plots
"logFile"        = "/scratch/AG_Akalin/bosberg/nanopore/pipeline_output/20180417_1233_HEK293_polyA_RNA/08_testing/test.log",                             # self-explanatory
"assembly"       = "hg19",
"poremodel_ref"  = "/clusterhome/bosberg/projects/nanopiper/dev/ref/pore_model_table.csv",
"sampleName"     = "HEK293_untreated_testsession",
"mincov_in"      = 10,
"plotrange_in"   = 10
)
group="CIMS"

# ======= DEBUGGING: DOWN TO HERE: ========

# overlaps CONVENTION: findoverlaps ( region of interest, reads        )
#                                      < queryHits >    | < subjectHits >
#                                     __________________|______________
#                                       locus_i         |   read_i
#                                       ...             |   ...

# ========================================
# Read in inputs:
ref_Genome    <- readDNAStringSet( argsL$pathin_refGen )
readdat       <- readRDS( argsL$pathin_reads )
putloci       <- readRDS( argsL$pathin_RsoI  )
poremodel     <- read.table( file = argsL$poremodel_ref,
                             sep = "\t",
                             header = TRUE
                             )
# ========================================
# Read in inputs:

# reduce list of RsoI to those with sufficient coverage.
loci_filtered_for_coverage <- filter_loci_for_coverage (  loci   = putloci$Region_groups,
                                                          reads  = readdat$aligned_reads,
                                                          mincov = argsL$mincov_in
                                                        )

# NOW check overlaps of reads with only the _covered_ loci.
#TODO: We run this overlaps more than once; if speedup is required later, eliminate this repetition.
overlaps_by_group = lapply ( names( loci_filtered_for_coverage ), function(group)
                               findOverlaps(  loci_filtered_for_coverage[[group]],
                               readdat$aligned_reads[[group]] )
                               )
names( overlaps_by_group ) <- names( loci_filtered_for_coverage )
# overlaps_by_group queryHits now references the indices of COVERED loci.

i=9
sampleROI_dat <-  get_sampledat_over_ROI ( SampleName         = argsL$sampleName,
                                           overlapping_reads  = readdat$aligned_reads[[group]][
                                                                 subjectHits( overlaps_by_group[[group]][
                                                                     queryHits( overlaps_by_group[[group]] ) == i ] ) ],
                                           # Comment explaining the above lines.
                                           # for plotting, we want to collect:
                                           #                     ^ The subset of the reads
                                           #                        ^ That are referenced as the subject
                                           #                           ^ in an overlap-pair for which the query is the one we are looking at (i.e. ="i").

                                           refgen             = ref_Genome,
                                           ROI_raw            = loci_filtered_for_coverage[[group]][i],
                                           plotrange          = argsL$plotrange_in,
                                           poremodel_ref      = poremodel
                                         )

plot_samplesignal_over_ROI( sampleROI_dat  = sampleROI_dat,
                            refgen         = ref_Genome)

plot_dwelltime_over_ROI( sampleROI_dat = sampleROI_dat,
                         refgen        = ref_Genome,
                         log           = TRUE )

# =======================================
i=0
i=i+1; plot_lt_passage_over_ROI   ( SampleName = argsL$sampleName,
                             overlapping_reads  = readdat$aligned_reads[[group]][
                                                                 subjectHits( overlaps_by_group[[group]][
                                                                     queryHits( overlaps_by_group[[group]] ) == i ] ) ],
                                         # Comment explaining the above lines.
                                         # for plotting, we want to collect:
                                         #                    ^ The subset of the reads
                                         #                       ^ That are referenced as the subject
                                         #                          ^ in an overlap-pair for which the query is the one we are looking at (i.e. ="i").

                             refgen        = ref_Genome,
                             ROI_raw       = loci_filtered_for_coverage[[group]][i],
                             plotrange =  10
                           )

