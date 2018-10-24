# modification_analysis.R
# ---Driver script for mod analysis.

dataset = "/scratch/AG_Akalin/bosberg/nanopore/pipeline_output/20180417_1233_HEK293_polyA_RNA/GRdat_allevents.RDS"
RDSdat <- readRDS(dataset)

unknown_mark = "N"

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

hg19_ref <- readDNAStringSet("/scratch/AG_Akalin/refGenomes/hg19_canon/hg19_canon.fa")


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
xcirc  = cos(angle)
ycirc  = sin(angle)
lines( xcirc, ycirc, col="black")


kspace_normed_current_matrix <- get_kspace_normed_current_vectors ( readsdat_in         = allevents_GR_called ,
                                                                    putlocs_GR_in       = CITS_put[ most_covered_CITSput_indices ],
                                                                    k                   = 5,
                                                                    breaks_in           = current_hist_breaks
                                                                   )
