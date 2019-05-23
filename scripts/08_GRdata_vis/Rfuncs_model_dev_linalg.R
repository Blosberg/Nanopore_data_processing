# ==== Script functions for analysing a RoI, and the degree to which overlapping
# reads deviate from the standard model
# ================================================================================

# --- filter a GRanges object for overlap with a putative location of intererest
# --- long-term goal: make this take sets of locations
get_kspace_normed_current_vectors <- function( readsdat_in           = stop("GRanges must be provided"),
                                               putlocs_GR_in         = stop("putative location of interest must be provided"),
                                               poremodel_statconsts  = stop(" matrix of stat params for the model must be provided"),
                                               k                     = 5,
                                               breaks_in             = seq(50,150,1),
                                               method                = "pca",
                                               add_synthdat          = FALSE
                                               )
{
  check_sanity( putlocs_GR_in )

  #### LAPPLY LEVEL HERE TO DO EVERY putloc ######
  overlaps    = findOverlaps( readsdat_in,
                              putlocs_GR_in  )       # ignore.strand is FALSE by default, so strands are matched.
#  @@@ execution per i-th ROI can be done from here

  pcalist       = list()

  ROIi= 1
  pcalist[[ROIi]] <- get_normed_vector_data( readsdat_loc_in = readsdat_in[ queryHits( overlaps [ subjectHits(overlaps) == ROIi ] ) ],
                                             putloc_GR_in    = putlocs_GR_in[ ROIi ] )


  testROIs = GRanges( seqnames = seqnames(putlocs_GR_in),
                      ranges = IRanges( start = start(putlocs_GR_in)+20, end = end(putlocs_GR_in)+20 ),
                      strand = strand(putlocs_GR_in) )

  testoverlaps = findOverlaps( readsdat_in,
                               testROIs   )
  testpcalist = list()

  testROIi=3
  testpcalist[[ROIi]] <- get_normed_vector_data( readsdat_loc_in = readsdat_in[ queryHits( testoverlaps [ subjectHits(testoverlaps) == testROIi ] ) ]  ,
                                                 putloc_GR_in    = testROIs[ testROIi ] )



}

# ==================================================================
# --- get scatter plot deviation for a particular ROI:

get_normed_vector_data <-function ( readsdat_loc_in = stop("reads overlapping this location must be provided"),
                                    putloc_GR_in    = stop("single specific ROI must be provided") )
{


  # FOR NOW JUST DO THE FIRST LOCATION: @@@
  # Reads_for_this_loc        = readsdat_in[ queryHits( overlaps[ subjectHits(overlaps)==1 ] ) ]


  region_dat <- get_starting_positions_and_sequences ( ROI_GR_in      = putloc_GR_in,
                                                       covering_reads = readsdat_loc_in,
                                                       k              = k
                                                      )
  sequence_list <- unlist( region_dat$model_sequences[ as.character(region_dat$pos_array ) ] )

  # for a given location of interest, split up the overlaps by read:
  Reads_on_loc_splitby_read = split( readsdat_loc_in,
                                     readsdat_loc_in$read_index )

  # grab list of vectors of mean current values
  # (Each entry in the list corresponds to a read)
  all_overlapping_reads_as_kvecs     <- lapply( Reads_on_loc_splitby_read,
                                                function(x) convert_single_read_overlap_to_kdimvec( readloc_GR_in  = x,
                                                                                                    startpos_array = region_dat$pos_array )
                                                )


  # convert that list of vectors into a matrix
  # each column is a vector of mean current values
  all_overlapping_reads_as_kvecs_mat              <- simplify2array( all_overlapping_reads_as_kvecs)

  # rename the rows according to the model_kmer that corresponds to this window
  row.names( all_overlapping_reads_as_kvecs_mat ) <- sequence_list

  #@@@@@
  synthdat_raw <- get_synthdat ( N                   =  50,
                                 seqs                =  sequence_list,
                                 events_splitby_kmer =  readsdat_splitby_kmer )
  synthdat_normed <-  (all_overlapping_reads_as_kvecs_mat-ROI_statconsts[,1] )/( ROI_statconsts[,2] )
  # all_overlapping_reads_as_kvecs_mat <- cbind(all_overlapping_reads_as_kvecs_mat, synthdat)
  #@@@@@@

  # import the statistical constants for the relevant sequences:
  ROI_statconsts <- poremodel_statconsts[ sequence_list ,]

  all_olapreads_normedkvecs <- ( all_overlapping_reads_as_kvecs_mat-ROI_statconsts[,1])/( ROI_statconsts[,2] )

  # For the moment, treat NA's as zero deviation from the model.
  all_olapreads_normedkvecs[ is.na( all_olapreads_normedkvecs) ] <- 0

  # grab the principle components
  pca <- prcomp( t(all_olapreads_normedkvecs),
                 center = FALSE,
                 scale. = FALSE)
  # pca$rotation %*% t(pca$x) -> gives you back the original all_olapreads_normedkvecs

  rot_in2_pca_mat = solve ( pca$rotation )
  # "solve" in R means "invert", because R is terrible.

  synthdat_normed_rot2pca = rot_in2_pca_mat %*% synthdat_normed
  # these are now _row_ -based vectors

  plot( pca$x[,1], pca$x[,2],
        col  = "red",
        xlab = "PC 1",
        ylab = "PC 2",
        main = "Deviation from model" )

  points( synthdat_normed_rot2pca[1,],
          synthdat_normed_rot2pca[2,],
          col  = "blue",
          pch  = 2 )

  pi=3.141592358979;  angles = seq(0, 2*pi, 0.001); xcirc  = cos(angles); ycirc  = sin(angles)
  lines( xcirc, ycirc, col="black")

  return(pca)
}
# ===============================================================
# create a set of vectors orthonormal to an input:
gen_orthonorm_vecs <- function(vec_in)
{
  temp=dim(vec_in)

  if(temp[1] <= 1 || temp[2] != 1)
  {
    stop(paste("vector input of dimension",temp[1]," x ", temp[2],"; incompatible with orthogonal generation") )
  }
  N=temp[1]

  result     = matrix(0,N,N)
  result[,1] = vec_in/(norm(vec_in,type="F"))

  for(i in c(2:N))
  {
    temp_a   = matrix(0,N,1)
    temp_b   = matrix(0,N,1)

    temp_a[1:(i-1)] = 1
    temp_a[i]       = 2

    for( j in c(1:(i-1) ) )
    {
      temp_a = temp_a - proj(result[, j, drop = FALSE],temp_a) #-- subtract projections onto each other dimension
    }

    result[,i] = temp_a/(norm(temp_a,type="F"))
  }

  return(result)
}

