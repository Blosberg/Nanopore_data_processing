# modification_analysis.R
# ---Driver script for mod analysis.

library(Gviz)
library(GenomicRanges)
library(BSgenome.Hsapiens.UCSC.hg19 )

assembly   = "hg19"
mincurrent = 50
maxcurrent = 250

SampleName = "HEK293_unmod"

Refgenome_path <- "/scratch/AG_Akalin/refGenomes/hg19_canon/hg19_canon.fa"
hg19_ref       <-  readDNAStringSet( Refgenome_path )
hg19_strack    <-  SequenceTrack(hg19_ref)

HEK293_MM2po_reads <- readRDS( "/fast/AG_Akalin/bosberg/nanopore/pipeline_output/20180417_1233_HEK293_polyA_RNA/07_GRprocessing/HEK_untreated_read_ROIolap_MihaM_2pO.rds" )
Miha_2po_sites <- readRDS("/fast/AG_Akalin/bosberg/nanopore/ref/Regions_of_interest/MihaM_2po_Ill.rds")


Ni = 63 # arbitrarily chosen: just need to plot something
Locus <- Miha_2po_sites$Region_groups$Ill_2pO[Ni]
start( Locus  ) <- 26227566
end(   Locus  ) <- 26227586 
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
xlabels = unlist( strsplit( as.character( hg19_ref$chr1[ start(Locus): end(Locus) ] ), 
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




