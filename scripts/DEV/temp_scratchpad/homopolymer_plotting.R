# load the HEK data
dataset = "/scratch/AG_Akalin/bosberg/nanopore/pipeline_output/20180417_1233_HEK293_polyA_RNA/GRdat_allevents.RDS"
RDSdat <- readRDS(dataset)


xmin = 50
xmax = 150
res  = 0.5
current_hist_breaks  <- seq( xmin, xmax, res)
# breaks in current level that define our "binning" for histogram plots.

unknown_mark = "N"

# ===========================================

allevents_GR_uncalled <-  RDSdat$allevents_GR[ which(   grepl( unknown_mark, 
                                                               RDSdat$allevents_GR$model_kmer) )  ]
allevents_GR_called   <-  RDSdat$allevents_GR[ which(  !grepl( unknown_mark, 
                                                               RDSdat$allevents_GR$model_kmer) )  ] 

#--- select only model kmers that don't contain "N"

GRL_splitbymodelkmer <- split(  allevents_GR_called, 
                                allevents_GR_called$model_kmer )  #--- create a list separated by model_kmer values.
                                # kmer5_freq <- table( reads_GR_called$model_kmer ) 

current_histlist_by_kmer <- lapply( GRL_splitbymodelkmer, function(x) get_normalized_current_hist( GR_input = x,
                                                                                                   breaks   = current_hist_breaks ) ) 
# list of 4^k = 1024 histograms of current values --one for each sequence.

temp <- hist( GRL_splitbymodelkmer[["AAAAA"]]$event_mean, 
              breaks=current_hist_breaks,
              plot = FALSE)
xdat <- temp$mids # --- the positions along the x-axis at which we plot the binned current values.


poremodel_statconsts = read.csv( "scripts/ref/pore_model_table.csv",
                                  stringsAsFactors = F,
                                  header           = TRUE,
                                  row.names        = 1,
                                  sep= "\t"
                                 )
                                   
# ====================================================================

seq = "GGATC"

poremodel_statconsts[seq, ]
poremodel_statconsts["GGATC", ]

plot_hist_against_normal ( xdat_in     = xdat,
                           hist_in     =  current_histlist_by_kmer[[seq]],
                           mean_in     = poremodel_statconsts[seq, 1],
                           stddev_in   = poremodel_statconsts[seq, 2],
                           col_in      = rgb( 0,      0,     1,   alpha=0.5 )
                          )


cumulative_NULLseq_hist =   trace_over_unknown_base ( histlist_in  = current_histlist_by_kmer,
                                                      sequence_in  = "NNNNN"  ) 
#--- this histogram is taken over all 1024 k-mer sequences, equally weighted
  
col_A = rgb( 0,           0,         1,   alpha=0.5  )
col_C = rgb( 0.3548, 0.4484, 0.1967742,   alpha=0.5  )
col_G = rgb( 1,         0.5,         0,   alpha=0.6 )
col_T = rgb( 1,           0,         0,   alpha=0.5  )

plot_background_hist( xdat, 
                      cumulative_NULLseq_hist )

plot_homopolymers( xdat_in     = xdat,
                   histlist_in = current_histlist_by_kmer,
                   pos_N       = 5,
                   scale_HP    = 0.1, 
                   col_A       = col_A,
                   col_C       = col_C,
                   col_G       = col_G,
                   col_T       = col_T,
                   isRNA       = TRUE   )

#---- this function can take any "seq_in" and plot the corresponding histogram 
#---- (also accepts N-values which are then traced over)
plot_specific_sequence_currenthist ( xdat_in      = xdat,
                                     histlist_in  = current_histlist_by_kmer,
                                     seq_in       = "GTATC",
                                     scale_HP     = 0.15,
                                     col_in       = rgb( 1,      0,     0,   alpha=0.5 ),
                                     add          = TRUE
                                     )
  

# ===================================================


  
# ====== PLOTTING =======
scale_HP = 0.15
i=4 

seq = High_coverage_reads_plist[[i]][1]$model_kmer

plot_background_hist( xdat, 
                      cumulative_NULLseq_hist )
plot_specific_sequence_currenthist ( xdat_in      = xdat,
                                     histlist_in  = current_histlist_by_kmer,
                                     seq_in       = seq,
                                     scale_HP     =scale_HP,
                                     col_in       = rgb( 1,      0,     0,   alpha=0.5 ),
                                     add          = TRUE
)

polygon( c( min(xdat_in),
            xdat_in,
            max(xdat_in)
           ),
          scale_HP * c( 0,
                        High_coverage_reads_phistlist[[i]],
                        0 ),
        col  = rgb( 0,      0,     1,   alpha=0.5 ),
        lty  = "blank"
)

