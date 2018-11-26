# ===============================================================
# look for another read which, together with rps16, will provide maximal kmer coverage:
get_targetcov <- function(  read_GR      = stop("reads_GRL must be provided"),
                            target_kmers = stop("RsoI must be provided"),
                            return_seqs  = FALSE
                          )  
{
  L= unique(  nchar( target_kmers ))
  if ( length( L ) != 1 )
       { stop("ERROR: non-uniform length of target_kmers") }
  # ----------
  
  target_kmers_found <- target_kmers [ !is.na( match(  target_kmers ,  substr( read_GR$model_kmer, 1,  L) )  ) ]
  
  if ( return_seqs )
  { return(target_kmers_found) }
  
  complimentarity    <- length( target_kmers_found )
  
  # location <- reduce(read_GR)
  return( complimentarity )
}
# ===============================================================
# --- sort through a sequence and check its diversity:
get_kmer_divers <- function( seq_in = stop("GRanges must be provided"),
                             n      = stop("n must be provided") 
                            )
{
N=nchar(seq_in)  
  
kmerset = list()

for ( p in c(1:N-n+1))
  {
  kmerset[p]  <- substr( seq_in, p, p+n-1 )
  }

kmerset_array <- unique( unlist( kmerset ) )

return( kmerset_array)
}
