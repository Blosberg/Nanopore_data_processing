# modification_analysis.R
# ---Driver script for mod analysis.

library(Gviz)
library(GenomicRanges)
library(BSgenome.Hsapiens.UCSC.hg19 )

assembly       <- "hg19"

Refgenome_path <- "/scratch/AG_Akalin/refGenomes/hg19_canon/hg19_canon.fa"
hg19_ref       <-  readDNAStringSet( Refgenome_path )
hg19_strack    <-  SequenceTrack(hg19_ref)

HEK293_MM2po_reads <- readRDS( "/fast/AG_Akalin/bosberg/nanopore/pipeline_output/20180417_1233_HEK293_polyA_RNA/07_GRprocessing/HEK_untreated_read_ROIolap_MihaM_2pO.rds" )
Miha_2po_sites <- readRDS("/fast/AG_Akalin/bosberg/nanopore/ref/Regions_of_interest/MihaM_2po_Ill.rds")

Ni = 63 # arbitrarily chosen: just need to plot something

# take reads that were aligned and putative modification sites:
untreated_olaps_w_2po <- findOverlaps( HEK293_MM2po_reads$aligned_reads$Ill_2pO, Miha_2po_sites$Region_groups$Ill_2pO )

locus_i_olaps <- untreated_olaps_w_2po[ subjectHits(untreated_olaps_w_2po) == Ni ]
# take just the reads that overlap:

plot_leading  = 10
plot_trailing = 10

Locus <- Miha_2po_sites$Region_groups$Ill_2pO[Ni]
start( Locus  ) <- start( Locus ) -10
end(   Locus  ) <- end(   Locus ) +10

reads_of_interest <- HEK293_MM2po_reads$aligned_reads$Ill_2pO[ queryHits( locus_i_olaps )  ]

# Shows the chromosome context overall. Don't really need this
# itrack <- IdeogramTrack(genome=genome, chromosome=chr)
# plotTracks(list(itrack, gtrack, atrack) )

# establish direction (based on strand info)
atrack <- AnnotationTrack(Locus, name="MM2po")
plotTracks(atrack)

# get Genome coordinates:
gtrack <- GenomeAxisTrack()
plotTracks(list(gtrack, atrack))

# get the sequence reference:
strack <- SequenceTrack(Hsapiens, chromosome=chr)

# now get the current data:
readofInterest <- reads_of_interest[[1]]

test = readofInterest[1:25]

dtrack <- DataTrack( data=readofInterest $event_mean,
                     start=start(readofInterest ),
                     end= end(readofInterest ),
                     chromosome= seqnames( readofInterest ),
                     genome=assembly,
                     name="ReadofInterest" )

plotTracks( c(gtrack, atrack, dtrack, strack),
            from = start( Locus  ) ,
            to = end(   Locus  ),
            type = "a" )




#======================================


HEK_atrack  <- AnnotationTrack(HEK_rps16_signal, name = "HEK")

genome( signal) <- "hg19"

IVum_atrack <- AnnotationTrack(IVum_rps16_signal, name = "IVum")
genome(IVum_rps16_signal) <- "hg19"

chr <- as.character(unique(seqnames(HEK_rps16_signal)))
gen <- unique( genome(HEK_rps16_signal) )

itrack <- IdeogramTrack(genome = gen, chromosome = chr)

a = 39923778
b = 39923790


plotTracks(list(itrack, gtrack, HEK_atrack, grtrack, IVum_atrack),
            from = a, to = b,
           cex = 0.8)


HEK_locsignal  = HEK_rps16_signal[ start(HEK_rps16_signal) > a  ]
HEK_locsignal  = HEK_locsignal[ start(HEK_locsignal)       < b  ]

IVum_locsignal  = IVum_rps16_signal[ start(IVum_rps16_signal) > a  ]
IVum_locsignal  = IVum_locsignal[ start(IVum_locsignal)       < b  ]



plot( start(IVum_locsignal),
      IVum_locsignal$event_mean,
      main ="Current vs. sequence",
      col    = "red",
      type   = "l",
      lwd    = 2,
      xaxt   = "n",
      xlab   = " chr19: 39925572 - 39925660 ",
      ylab   = "current [pA]",
      ps =0.5
     )

axis(1, at=start(IVum_locsignal),
     labels=substr( IVum_locsignal$model_kmer, 1, 1),
     ps = 0.5)

lines( start(IVum_locsignal), IVum_locsignal$event_mean - IVum_locsignal$event_stdv,
      col=rgb(1,0,0,alpha=0.25), lwd = 0.5)
lines( start(IVum_locsignal), IVum_locsignal$event_mean + IVum_locsignal$event_stdv,
      col=rgb(1,0,0,alpha=0.25), lwd = 0.5)

lines( start(HEK_locsignal), HEK_locsignal$event_mean,
      col=rgb(0,0,1,alpha=1), lwd = 2)
lines( start(HEK_locsignal), HEK_locsignal$event_mean - HEK_locsignal$event_stdv,
      col=rgb(0,0,1,alpha=0.25), lwd = 0.5)
lines( start(HEK_locsignal), HEK_locsignal$event_mean + HEK_locsignal$event_stdv,
      col=rgb(0,0,1,alpha=0.25), lwd = 0.5)

legend( "topright",
        legend=c("IV signal", "HEK cells"),
        col   =c("red","blue"),lwd = 2

        )
#======================================


HEK_rps16_flattened_reads <- lapply( HEK_rps16_reads, function(x) flatten_read( read_GR_in = x ) )




reduced_IVum_chr19_splitby_read <- lapply( IVum_chr19_splitby_read,
                                          function(x) reduce(x) )

unknown_mark = "N"

poremodel_statconsts = read.csv( "scripts/ref/pore_model_table.csv",
                                  stringsAsFactors = F,
                                  header           = TRUE,
                                  row.names        = 1,
                                  sep= "\t"
                                 )



plotTracks(DataTrack(twoGroups, name = "a"), type="a")


data(twoGroups)

# ===========================================

HEK_splitby_kmer     = split( HEK_RDSdat$allevents_GR,     HEK_RDSdat$allevents_GR$model_kmer     )
IVunmod_splitby_kmer = split( IVunmod_RDSdat$Events_all_GR, IVunmod_RDSdat$Events_all_GR$model_kmer )

mean_current_HEK     <- lapply( HEK_splitby_kmer,
                                function(x) mean(x$event_mean) )
current_kmer_HEK = unlist( mean_current_HEK)

mean_current_IVum    <- lapply( IVunmod_splitby_kmer,
                                function(x) mean(x$event_mean) )
current_kmer_IVum = unlist( mean_current_IVum )

IVkmers = names( mean_current_IVum )

xmin <- min( unlist(mean_current_IVum) )
xmax <- max( unlist(mean_current_IVum) )

plot( mean_current_HEK[IVkmers],
      mean_current_IVum[IVkmers],
      xlab = "HEK293 (kmer-subset)",
      ylab = "In-vitro, unmodified (all kmers)",
      col  = "blue",
      main = "Direct current comparison: HEK293 vs. in-vitro"
      )
lines( c(xmin, xmax),
       c(xmin, xmax),
       col= rgb(0,0,0,alpha=0.75),
       lwd=5)





HEK_splitby_kmer[[1]]$model_kmer

# ===========================================


allevents_GR_uncalled <-  RDSdat$allevents_GR[ which(   grepl( unknown_mark,
                                                               RDSdat$allevents_GR$model_kmer) )  ]
allevents_GR_called   <-  RDSdat$allevents_GR[ which(  !grepl( unknown_mark,
                                                               RDSdat$allevents_GR$model_kmer) )  ]
poremodel_statconsts = read.csv( "scripts/ref/pore_model_table.csv",
                                  stringsAsFactors = F,
                                  header           = TRUE,
                                  row.names        = 1,
                                  sep= "\t"
                                 )


# --- NOW CONSIDER PUTATIVE MODIFICATION SITES  -----

PATH_putloc="/scratch/AG_Akalin/bosberg/nanopore/ref/Linder_2015_nature_supp/put_m6A_locs_CITS.csv"
CITS_put <- get_putlocs( PATH_putloc )


# just check for overlaps now to see which ones have the highest coverage
calledreads_CITS_pairs = findOverlaps(allevents_GR_called, CITS_put )
most_covered_CITSput_indices = as.numeric(
                                   names( sort(
                                               table(
                                                      subjectHits( calledreads_CITS_pairs )
                                                      ),decreasing=TRUE
                                               ) )
                                          )
CITSput_rankedbycoverage = CITS_put[ most_covered_CITSput_indices]

# ========

pi     = 3.141592358979
angles = seq(0, 2*pi, 0.001)
xcirc  = cos(angles)
ycirc  = sin(angles)
lines( xcirc, ycirc, col="black")


kspace_normed_current_matrix <- get_kspace_normed_current_vectors ( readsdat_in         = allevents_GR_called ,
                                                                    putlocs_GR_in       = CITS_put[ most_covered_CITSput_indices ],
                                                                    k                   = 5
                                                                  )
