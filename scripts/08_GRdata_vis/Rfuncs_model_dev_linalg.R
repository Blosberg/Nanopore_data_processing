# ==== Script functions for analysing a RoI, and the degree to which overlapping
# reads deviate from the standard model
# ================================================================================

# --- filter a GRanges object for overlap with a putative location of intererest
# --- long-term goal: make this take sets of locations
get_pca_normdiff  <- function( sampleROI_dat_in  = stop("GRanges must be provided"),
                               k                 = 5,
                               method            = "pca", # TODO: add a KL-div based method.
                               shouldplot        = TRUE
                              )
{

  Nxpos = dim( sampleROI_dat_in$read_normdiff )[2]

  # get the subset of data points near the centre of the ROI:
  ROI_dat_lin = sampleROI_dat_in$read_normdiff[ , c( (ceiling(Nxpos/2)-k+1) : (ceiling(Nxpos/2)+k-1) ) ]


  ROI_dat_lin[ is.na( ROI_dat_lin ) ] <- 0

  # grab the principle components
  pca <- prcomp( (ROI_dat_lin),
                 center = FALSE,
                 scale  = FALSE)

  # if you do pca <- prcomp ( A )
  # A_copy <- pca$rotation %*% t(pca$x)
  # ^ gives you back the original value "A"

  # to get the matrix that performs the forward operation (i.e. rotates from original space _into_ PCA space):
  # rot_in2_pca_mat = solve ( pca$rotation )
  # ("solve" in R means "invert", because R is terrible.)

  # pcax_copy <- t( rot_in2_pca_mat %*% t( ROI_dat_lin ) )
  # ^ This will be equal to pca$x (within machine tolerance)
  # synthdat_normed_rot2pca = rot_in2_pca_mat %*% synthdat_normed
  # these are now _row_ -based vectors


  clust_kmeans <- kmeans( x = ROI_dat_lin,
                          centers = 2,
                          iter.max = 10 )
  dev = sqrt( rowSums( clust_kmeans$centers * clust_kmeans$centers  )  )


  if( shouldplot )
  {
    plot( # pca$x[,1] * (1/sqrt(Nxpos)),
          # pca$x[,2] * (1/sqrt(Nxpos)),
          pca$x[,1],
          pca$x[,2],
          col  = "red",
          xlab = "PC 1",
          ylab = "PC 2",
          main = "Deviation from model" )

    # draw a unit circle to represent a single std. dev.
    pi=3.14159265358979;
    angles = seq(0, 2*pi, 0.001);
    xcirc  = cos(angles);
    ycirc  = sin(angles);
    lines( xcirc, ycirc, col="black")

   points( # pca$x[,1] * (1/sqrt(Nxpos)),
           # pca$x[,2] * (1/sqrt(Nxpos)),
           pca$x[ clust_kmeans$cluster == 1 ,1],
           pca$x[ clust_kmeans$cluster == 1 ,2],
           col  = "blue",
           lw   = 3 )

  }

  return(pca)
}

# ===============================================================
# Check for clustering:
# cluster_normdiff_reads <- function ( ... )

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

