# HANDLE THE CASE WHEN THE OVERLAP SUBSET IS EMPTY:
# ==================================================================

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
# --- grab the "starting positions" for each kmer when you want the pos of interest to be in the 'ith' position.
# --- N.B. "Start" means lowest valued x-position in the range --> *REGARDLESS* of strand or direction of transcription.
get_starting_positions_and_sequences <- function( ROI_GR_in      = stop("GR_in must be provided"), 
                                                  covering_reads = stop("covering_reads must be provided"),
                                                  k              = 5               )
{
# We want the position array in order such that the 'i'th 
# entry corresponds to whenever the modification is in position i
# e.g. pos[1] is the starting position of whatever window leads to 
# the modification of interest being in position 1 of its respective strand.
    
  if( as.character( strand(ROI_GR_in)) == "+")  
    { pos_offset = c(0:(k-1)) } else if( as.character( strand( ROI_GR_in) ) == "-" ){
      pos_offset = c((k-1):0) } else { 
      stop("strand ill-defined.")
    }
  pos_array = start(ROI_GR_in) - pos_offset

  reads_splitby_startpos <- split( covering_reads, start(covering_reads)  )
  if( length(reads_splitby_startpos) < k ) # Must have at least one read cover every position.
    { return(NULL) } else{
      model_sequences        <- lapply( reads_splitby_startpos, function(x) x[1]$model_kmer ) 
    }
  #------------
  
  return( list("pos_array" = pos_array, 
               "model_sequences" = model_sequences) 
          )
}

# ===============================================================
# --- given a histogram of data, and the associated bins, plot it alongside a normal distribution
plot_hist_against_normal <- function( xdat_in     = stop("xdat_in must be provided"),
                                      hist_in     = stop("hist_in must be provided"),
                                      mean_in     = 0,
                                      stddev_in   = 1,
                                      col_in      = rgb( 0,      0,     1,   alpha=0.5 ),
                                      scale_HP    = 1
                                     )
{

  xnormed = ( (xdat_in - mean_in) /stddev_in );
  dX = xnormed[2]-xnormed[1]
  
  if ( max( abs(diff(xnormed) -dX)/dX  ) > 0.01)
    {stop("non-uniform spacing in histogram x-data")}
  
  plot( xnormed, 
        dnorm( xnormed, mean=0, sd=1), 
        type="l", 
        lwd=1,
        ylab = " probability density ",
        xlab = "normalized deviation",
        ylim = c(0,0.5)
        )
  
    polygon( c( min(xnormed),
                xnormed,
                max(xnormed)
    ),
    scale_HP * (1/dX) * c( 0,
                           hist_in, # N.B. this need not _necessarily_ contain unknown bases (i.e. "N"'s -- but it CAN handle them.)
                           0),
    col  = col_in,
    lty  = "blank"
    )
      
}
# ===============================================================
# --- take k-dimentional data and rotate it into PC 1, 2 
# Rotate_kD_data_to_pc12_hist  <- function( k=5 )
# {
#   ...
#   }
# # ===============================================================
# # --- convert each overlapping read into a k-dim vector of normalized deviations
# # --- FOR each read, get a k-dimensional vector describing its position. 
# convert_GRL_to_readsegged_kdim  <- function( GR_in = stop("GR_in must be provided"),
#                                              k=5, 
#                                              )
# {
#   ...
#   }

# --- Given a specific read overlapping with a particular position, 
# --- construct a vector of normalized deviations from the model current  
convert_single_read_overlap_to_kdimvec  <- function( readloc_GR_in  = stop("GR_in must be provided"),
                                                     startpos_array = stop("start array must be provided"),
                                                     k              = 5
                                                    )
{ # --- startpos_array is a k-dimensional array of starting positions[i] such that
  # --- the "i"-th position is the putative region of interest
  # --- it can be increasing or decreasing depending on the strand, 
  # --- but the "start" position is always the point of reference.

   readoverlap_splitby_startpos <- split(  readloc_GR_in, 
                                           start(readloc_GR_in) )  
  
#   if ( length(readoverlap_splitby_startpos) < k )
#     {print("WARNING: not all positions covered.") } # turn this into a "return(na)" statement.
   
   read_pos_current_mean <- lapply( readoverlap_splitby_startpos, 
                                    function(x) mean(x$event_mean) ) 
   
   return(  unlist(read_pos_current_mean)[ as.character(startpos_array) ] )
   # ---- AVERAGE CURRENT VALUES ARE RETURNED AS A VECTOR WITH ELEMENT NAMES CORRESPONDING TO POSITION.
   # ---- at this stage, normalization wrt. the model is *NOT YET* performed 
   }
# # ===============================================================
# # --- just take a histogram of data, and plot it alongside a normal distribution
#   
#   if( as.character( strand(GR_in)) == "-")  
#     { pos_offset = c(0:(k-1)) } else if( as.character( strand(GR_in) ) == "+" ){
#       pos_offset = c((k-1):0) } else { 
#       stop("strand ill-defined.")
#     }
#     
#   pos_array = start(GR_in) - pos_offset 
#   return(pos_array)
# }

# ===============================================================
check_sanity <- function( putlocs_GR_in = stop("RsoI must be provided")  )  
{
  if( max( abs( start( putlocs_GR_in ) - end( putlocs_GR_in )  ) ) !=0 )
    { stop(" putative locations should be single-base positions, not ranges.")}
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

# ===============================================================
# Import putative modification sites: 
get_putlocs <- function( PATHin = stop("PATHin must be provided"))
{
  
  CITS_rawdat=read.csv( PATHin, sep = "\t", header = TRUE )
  
  CITS_put =  GRanges( seqnames = CITS_rawdat$Chr ,
                       ranges   = IRanges ( start  = CITS_rawdat$Start,
                                            end    = CITS_rawdat$End ),
                       strand   = CITS_rawdat$Strand
  )
  
  return( sort( CITS_put[ which( seqnames(CITS_put) != "chrM") ] )  ) # 6543 locations
}


# ===============================================================
# look for the bin with the highest counts within a certain range of a given a hist( ..., plot = false) dataset, .
find_maxfreq_length <- function(  histdat   = stop("histdat must be provided"),
                                  a = stop("a must be provided"),
                                  b = stop("b must be provided")
                           )  
{

subset_1 <- which(histdat$mids > a)
subset_2 <- subset_1[ which(histdat$mids[subset_1] < b) ]
return( histdat$mids[ subset_2[ which.max( histdat$counts[ subset_2 ]  )  ] ] )
}
