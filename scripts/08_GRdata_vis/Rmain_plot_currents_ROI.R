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
"pathin_reads"   = "/Users/lana/Desktop/train_dat/olap_files/HEK_untreated_read_ROIolap_m6A_put.rds",  # path to Granges List structure with reads aligned to reference genome
"pathin_RsoI"    = "/Users/lana/Desktop/train_dat/putmod_locs/m6A_putlocs_Linder.rds",                 # path to Granges list of regions of interest in the current study
"pathin_refGen"  = "/Users/lana/Desktop/train_dat/refGenome/hg19_canon.fa.gz",                         # path to Granges list of regions of interest in the current study
"pathout_plot"   = "/Users/lana/Desktop/train_dat/test_out/localtest.out",                             # where should we send the plots
"logFile"        = "/Users/lana/Desktop/train_dat/test_out/localtest.log",                             # self-explanatory
"assembly"       = "hg19",
"sampleName"     = "HEK293_untreated",
"mincov_in"      = 10
)
suppressPackageStartupMessages( library(BSgenome.Hsapiens.UCSC.hg19 ) )
ref_Genome    <- BSgenome.Hsapiens.UCSC.hg19
group="CIMS"

# ======= DEBUGGING: DOWN TO HERE: ========

# ========================================
# Read in inputs:
ref_Genome    <- readDNAStringSet( argsL$pathin_refGen )
readdat       <- readRDS( argsL$pathin_reads )
putloci       <- readRDS( argsL$pathin_RsoI  )

# ========================================
# Read in inputs:

# TODO: add sanity check that reference names agree.

read_RoI_filtered_for_coverage <- filter_loci_coverage (  loci   = putloci$Region_groups,
                                                          reads  = readdat$aligned_reads,
                                                          mincov = mincov_in
                                                        )
readdat$aligned_reads


chr_Locus = as.character( seqnames( Locus) )
strand_Locus = as.character( unique( strand( read) ) )

plotframe = GRanges( seqnames = chr_Locus,
                     strand   = strand_Locus,
                     ranges   = IRanges( c( start(Locus):end(Locus) ) ,
                                         end =  c( start(Locus):end(Locus) )  )
                     )

# take reads that were aligned and putative modification sites:
reads_of_interest <- findOverlaps( HEK293_MM2po_reads$aligned_reads$Ill_2pO, Locus )

# ================= ???? ======================

i=1
read = reads_of_interest[[i]]
end(read) <- start(read)

olaps = findOverlaps( plotframe, read )

# plot a single read:
plot( start(read[ subjectHits(olaps)] ),
                read[ subjectHits( olaps )  ]$event_mean,
      type = "l",
      xlim = c( start(Locus),
                end(Locus)),
      main= paste(  SampleName, "\n",
                    as.character( unique( seqnames(Locus) ) ),
                    ":",
                    as.character(start(Locus) ),
                    "-",
                    as.character( end(Locus)  ) ),
      xlab = paste( "+ strand sequence." ),
      xaxt="n",

      ylim = c(mincurrent, maxcurrent),
      ylab = "current[pA]",
      col  = rgb( 0, 0, 1, 0.5 )
      )

# prepare the sequence tick-mark labels:
xcoords = seq( start(Locus), end(Locus), 1 )
xlabels = unlist( strsplit( as.character( ref_Genome$chr1[ start(Locus): end(Locus) ] ),
                   split = "" ) )

# apply the sequence tick-mark labels:
for ( i  in c(1:length(xcoords)) )
  { axis(1,
       at = xcoords[i],
       xlabels[i]
       ) }

#-----------add lines for subsequent reads (keeping the same boundaries, etc.):

i=100
read = reads_of_interest[[i]]
end(read) <- start(read)

olaps = findOverlaps( plotframe, read )

# plot a single read:
lines( start(read[ subjectHits(olaps)] ),
                read[ subjectHits( olaps )  ]$event_mean,
      col  = rgb( 0, 0, 1, 0.5 )
      )




