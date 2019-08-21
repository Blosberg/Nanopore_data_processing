suppressPackageStartupMessages( library(GenomicRanges) )
suppressPackageStartupMessages( library(Biostrings) )
suppressPackageStartupMessages( library(cluster) )


                  
# --------------------------------------------------
SampleName = "HEK293_untreated"
allReadDat_PATH = "/fast/AG_Akalin/bosberg/nanopore/pipeline_output/20180417_1233_HEK293_polyA_RNA/06_GRobjects/HEK_untreated_reads_GRL.rds"
allReadDat = readRDS( allReadDat_PATH )

Path_refgen = "/scratch/AG_Akalin/refGenomes/hg19_canon/hg19_canon.fa"
ref_Genome    <- readDNAStringSet( Path_refgen )


ROI_DIR="/fast/AG_Akalin/bosberg/nanopore/ref/Regions_of_interest/hg19/"
putloc_CITS  = readRDS( paste0( ROI_DIR, "ROI_m6A-CITS.rds" ) )
putloc_2pomG = readRDS( paste0( ROI_DIR, "ROI_Me2p0.rds") )
putloc_2pomG_NC50up = readRDS( paste0( ROI_DIR,  "ROI_Me2p0-NC50up.rds") )
putloc_2pomG_NC50dn = readRDS( paste0( ROI_DIR,  "ROI_Me2p0-NC50dn.rds") )


# putloc_2pomG_NC50dn$loci$seq_context <- unlist( lapply( c(1:length(putloc_2pomG_NC50dn$loci)), 
#                                                        function(locus) 
#                                                            get_sequence_at_reference( GRange_ob_in = putloc_2pomG_NC50dn$loci[locus], 
#                                                                                       refGen_in = ref_Genome
#                                                                                       )
#                                                 ) )
# Find the most common sequence:
putloc_CITS_loci_GGACT = putloc_CITS$loci[ putloc_CITS$loci$seq_context == "GGACT" ] 
# These data already have the "A" in the middle (3rd) position.
# --------------------------------------------------

allReadDat_allevents_3shifted <- unlist( allReadDat$Events_GRL_splitbyread ) 
end( allReadDat_allevents_3shifted )   <- end(   allReadDat_allevents_3shifted ) + 3
start( allReadDat_allevents_3shifted ) <- start( allReadDat_allevents_3shifted ) + 3
# allevents_GGACT <- allReadDat_allevents[ which( allReadDat_allevents$model_kmer == "GGACT" ) ]
# allevents_3shifted_subset_GGACT = allReadDat_allevents_3shifted[ allReadDat_allevents_3shifted$model_kmer == "GGACT" ]
# GGACTEvents_3shifted_allputloci_overlaps = findOverlaps ( allevents_3shifted_subset_GGACT,
#                                                           putloc_CITS$loci )

Events_3shifted_2pOputloci_overlaps     = findOverlaps ( allReadDat_allevents_3shifted,
                                                         putloc_2pomG$loci )
Events_3shifted_2pO_NC50uploci_overlaps = findOverlaps ( allReadDat_allevents_3shifted,
                                                         putloc_2pomG_NC50up$loci )


GGACTEvents_3shifted_onputloc = allevents_3shifted_subset_GGACT[ queryHits( GGACTEvents_3shifted_allputloci_overlaps ) ]

GGACTEvents_onputlocs_add1downstream_seq = unlist( lapply( c(1:length(GGACTEvents_3shifted_onputloc)), 
                                                       function(x)  
                                                         get_sequence_at_reference( GRange_ob_in = GGACTEvents_3shifted_onputloc[x], 
                                                                                    refGen_in = ref_Genome,
                                                                                    upstream_incl = 1,
                                                                                    dnstream_incl = 3)
                                                  ))


# this is now the GRanges object of all GGACT events on a locus that is followed by a "G" making a 6-mer of GGACTG
GGACT_onputloci_part_of_6mers <- GGACTEvents_3shifted_onputloc[  which( GGACTEvents_onputlocs_add1downstream_seq  == "GACTG" )   ]

temp_olap = findOverlaps( GGACT_onputloci_part_of_6mers,
                          putloc_CITS$loci  )
# These are the unique sites that are putative modifications and have a GGACTG footprint.
GGACTG_putloc_CITS = putloc_CITS$loci[ unique( subjectHits( temp_olap ) )  ]

GGACTG_putloc_downstream_by1 = unique( shift_GR_streamwise( GRange_ob_in = GGACTG_putloc_CITS, 
                                                            n            = 1, 
                                                            direction    = "down" 
                                                           )
                                     )

GGACTG_putloc_downstream_by1$seq_context = unlist( lapply( c(1:length(GGACTG_putloc_downstream_by1)), 
                                                       function(x)
                                                         get_sequence_at_reference( GRange_ob_in = GGACTG_putloc_downstream_by1[x],
                                                                                    refGen_in = ref_Genome,
                                                                                    upstream_incl = 2,
                                                                                    dnstream_incl = 2)
                                                  ))




olap_temp = findOverlaps( allReadDat_allevents_3shifted, 
                          GGACTG_putloc_downstream_by1
                         )
GACTGEvents_onputloci_part_of_6mers =   allReadDat_allevents_3shifted[ queryHits( olap_temp) ]                                                  


# ---- split up by read and location
GGACT_onputloc_splitby_readi_and_loc =  split( GGACT_onputloci_part_of_6mers, 
                                               (GGACT_onputloci_part_of_6mers$read_index)*start(GGACT_onputloci_part_of_6mers) 
                                               )
GACTG_events_ds1_from_putloc_splitby_readi_and_loc <- split ( GACTGEvents_onputloci_part_of_6mers, 
                                                              GACTGEvents_onputloci_part_of_6mers$read_index * start(GACTGEvents_onputloci_part_of_6mers) 
                                                             )

# ---- flatten the current average for each read-location set
currentvals_GGACT_onputloc_flattened <- lapply( 
                                                 GGACT_onputloc_splitby_readi_and_loc, 
                                                     function(posn)
                                                        sum( posn$event_mean * posn$event_length)/sum(posn$event_length) 
                                               )
currentvals_GACTG_dsby1_from_putloc <- lapply( 
                                               GACTG_events_ds1_from_putloc_splitby_readi_and_loc, 
                                                   function(posn)
                                                        sum( posn$event_mean * posn$event_length)/sum(posn$event_length) 
                                              )
# ----- extract read indices and look for matches
temp = unlist( GRangesList( lapply(GGACT_onputloc_splitby_readi_and_loc, 
                  function(ro) shift_GR_streamwise( GRange_ob_in = ro[1], 
                                                    n            = 1, 
                                                    direction    = "down" 
                                                    )  
              )))
temp$read_index = unlist( lapply( GGACT_onputloc_splitby_readi_and_loc, 
                              function(ro) 
                                 ro[1]$read_index
                              ))

GGACT_onputloc_shifted_tag <- start(temp)*temp$read_index


GACTG_dsby1_fromputloc_readloc_tags <- unlist( lapply(GACTG_events_ds1_from_putloc_splitby_readi_and_loc, 
                                                  function(ro) 
                                                      unique( ro$read_index * start(ro) ) 
                                           ) 
                                        )
                                 
match_set <- match( GACTG_dsby1_fromputloc_readloc_tags,
                    GGACT_onputloc_shifted_tag)


ypoints <- unlist( currentvals_GACTG_dsby1_from_putloc[  which(! is.na(match_set)) ] )
xpoints <- unlist( currentvals_GGACT_onputloc_flattened[ match_set[ which(! is.na(match_set)) ] ] )

normed_data_x <- ( xpoints-poremodel["GGAUC",1])/poremodel["GGAUC",2] 
normed_data_y <- ( ypoints-poremodel["GAUCG",1])/poremodel["GGAUC",2] 

 plot( normed_data_x, 
       normed_data_y,
       xlim = c(-4,3),
       ylim = c(-4,3),
       col = rgb(0,0,1,0.5),
       main = "CITS modification sites")
 
   pi=3.14159265358979;
    angles = seq(0, 2*pi, 0.001);
    xcirc  = cos(angles);
    ycirc  = sin(angles);
    lines(   xcirc,    ycirc, col="black", lty = 1)
    lines( 2*xcirc,  2*ycirc, col="black", lty = 2)
    lines( 3*xcirc,  3*ycirc, col="black", lty = 3)

    scatter_dat <- cbind( normed_data_x, normed_data_y )
 
     cluster_dat  <-       kmeans( x        = scatter_dat,
                                   centers  = 2 )

     points( normed_data_x[which(cluster_dat$cluster==1)],
             normed_data_y[which(cluster_dat$cluster==1)],
             col = rgb(1,0,0,0.5)
            )

z=sqrt( normed_data_x^2 + normed_data_y^2)

GGACT_shifted_onedown <- lapply( GGACT_splitby_readi_and_loc,
                                    function( pos )
                                        shift_GR_streamwise( GRange_ob_in = pos[1], direction = "down", n = 1 )
                                 )
GGACT_shifted_onedown <- unique( GGACT_shifted_onedown )


allReadDat_allevents_3shifted[ query]            
GGACT_shifted_onedown
            
#@@@@ 


             
for( posn in c(1: length(GGACT_splitby_readi_and_loc ))) 
  {
  
  
  }
  
twostep_splitby_read <- split( two_step_event_chain, 
                               two_step_event_chain$read_index )



#@@@@@


test <- unlist( lapply(c(1:length(putloc_CITS_loci_GGACT)),
                        function(x)
                          get_sequence_at_reference( GRange_ob_in = putloc_CITS_loci_GGACT[x],
                                                     upstream_incl = 1, 
                                                     dnstream_incl = 3, 
                                                     refGen_in = ref_Genome
                                                   )
                        ) )
putloc_6mers = putloc_CITS_loci_GGACT[ which( test == "GACTG") ]



GGACTG_p1_positions = shift_GR_streamwise( putloc_6mers,
                                          n = 1 ,
                                          direction = "down")

findOverlaps( allReadDat_allevents_3shifted, GGACTG_p1_positions )

#### @@@@ =============




olap_temp = findOverlaps( allReadDat$Events_GRL_splitbyread,  putloc_GGACTG_6mers )


reads_overlapping_GGACTG_CITS_putloc_6mers = allReadDat$Events_GRL_splitbyread[ queryHits( olap_temp  ) ]

 
# --------------------------------------------------
allReadDat_allevents_windowshifted = allReadDat_allevents
end( allReadDat_allevents_windowshifted )   <- end( allReadDat_allevents_windowshifted ) + 5
start( allReadDat_allevents_windowshifted ) <- start( allReadDat_allevents_windowshifted ) + 1

allEvents_windowshifted_allputloci_overlaps = findOverlaps ( allReadDat_allevents_windowshifted, 
                                                             putloc_CITS$loci )

All_events_windoshifted_that_do_overlap <- allReadDat_allevents_windowshifted[ queryHits( allEvents_windowshifted_allputloci_overlaps ) ]


# --------------------------------------------------
GGACT_putloci = putloc_CITS$loci[ putloc_CITS$loci$seq_context == "GGACT" ]


GGACT_putloc_overlaps = findOverlaps(
                                     allevents_GGACT,
                                     GGACT_putloci 
                                     )
shift=3
allevents_GGACT_shifted          <- allevents_GGACT
start( allevents_GGACT_shifted ) <- start( allevents_GGACT_shifted ) - shift
end( allevents_GGACT_shifted )   <- end(   allevents_GGACT_shifted ) - shift

GGACT_shifted_putloc_overlaps =
findOverlaps( allevents_GGACT_shifted,
              GGACT_putloci )

allevents_GGACT_shifted[ queryHits( GGACT_shifted_putloc_overlaps )   ] 

# --------------------------------------------------
j=4


test <- lapply( c(1:length(GGACT_putloci)),
                   function(x) get_sequence_at_reference( GRange_ob_in = GGACT_putloci[x], 
                                                          refGen_in    = ref_Genome
                                                         ) 
              )
        
# --------------------------------------------------


ROI_i = 1
locus = putloc_CITS$loci[ ROI_i ]

Nloci          = length( putloc_CITS$loci )
context_list_bioS <- lapply( c(1:Nloci), 
                  function ( li )
                       eval(parse( text= paste0(  "ref_Genome$", 
                                                  seqnames( putloc_CITS$loci[li] ),
                                                  "[", as.character( start(putloc_CITS$loci[li]) -2 ),
                                                  ":", as.character( end(putloc_CITS$loci[li]) +  2 ), "]"
                                                )
                                 )
                          )
                       )

strand_integer = rep( 1, Nloci )
neg_strands_inds =  which( as.character( strand( putloc_CITS$loci ) ) == "-" )   

revcomped_list <- lapply( neg_strands_inds, 
                            function( i )
                              reverseComplement( context_list_bioS[[i]] )
                          )
context_list_bioS[ neg_strands_inds ] <- revcomped_list

context_array = unlist(  lapply( c(1:Nloci), 
                                 function( li )
                                    as.character( context_list_bioS[[li]] )
                                   )   )

