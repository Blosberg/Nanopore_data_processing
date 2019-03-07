# These are functions for handling the conversion of bases to more generalized constructs:
# e.g., if you want to input ACGNW, and you want the "N" to trace over all of ACGT, and the "W" to trace over the "Weak" bases: AT, etc. See IUPAC naming convention for a comprehensive list.

# ==================================================================
# --- produce list of sequences that trace over unknown bases (default mapping: IUPAC)
sequence_trace <- function( sequence_in       = stop("sequence     must be provided"),
                            base_tracemap_in  = list(  
                              # --- Standard IUPAC notation: ----
                              "A" = "A",
                              "C" = "C",
                              "G" = "G",
                              "T" = "T",
                              "W" = c("A","T"),       # Strong
                              "S" = c("C","G"),       # Weak
                              "M" = c("A","C"),       # aMino
                              "K" = c("G","T"),       # Keto
                              "R" = c("A","G"),       # puRine
                              "Y" = c("C","T"),       # pYrimidine
                              "B" = c("C","G","T"),   # not-A
                              "D"  = c("A","G","T"),  # not-C
                              "H" = c("A","C","T"),   # not-G
                              "V" = c("A","C","G"),   # not-(T or U)
                              "N" = c("A","C","G","T")
                            )
)
{
  # how long is the sequence? 
  k=nchar( sequence_in )
  
  # break it up into individual characters
  seq_in_singlechars = unlist( strsplit(sequence_in, split="") )
  
  # Check how many map to explicit and trace bases:
  Charmap_trace    = match( seq_in_singlechars, unlist( names(base_tracemap_in) ) )
  
  # check that all base positions are recognized and accounted for.
  if ( length( Charmap_trace [!is.na(Charmap_trace ) ] ) != k )
  { stop("sequence contains unrecognized bases") }
  
  
  # now go through all k's 
  seqs <- c()
  seq_mapped  <- as.vector( base_tracemap_in[ seq_in_singlechars]  )
  
  for (i in 1:k) {
    
    if( i == 1) {
      
      seqs <- unlist(strsplit( seq_mapped[[i]], 
                               split = NULL)) # first entry will not need appending
      
    } else {
      
      # now to that previous vector we add all permutations of the current character [i]
      seqs_new <- expand.grid(seqs, 
                              unlist(strsplit(seq_mapped[[i]], 
                                              split = NULL))) 
      
      #collapse into strings
      seqs_new <- apply(seqs_new,1,paste0, collapse = "")
      seqs <- seqs_new #update the seqs vector
      
    }
  }
  
return(seqs)
  
}
# ==================================================================
# --- take average of bases in a given position
trace_histogram <- function( histlist_in  = stop("reads must be provided"),
                             sequence_in  = stop("sequence must be provided"),
                             return_as_histlist_element = FALSE
                             )
{
  k             = nchar(sequence_in) # number of bases under consideration.

  sequence_list = sequence_trace( sequence_in  = sequence_in  )
  nseqs = length( sequence_list)

  # TODO: add a sanity cheeck to make sure all histlists are using the same x-binning.
  dx_arr = diff( histlist_in[[ 1 ]]$breaks )
  
  cumulative_hist=matrix( 0, length(dx_arr ), 1)
  for( seq in sequence_list )
    {
    cumulative_hist = cumulative_hist + (histlist_in[[seq]]$density * dx_arr )
  }
  
  NORM = sum(cumulative_hist);
  if( abs((NORM - nseqs)/nseqs) > 0.0001 )
    {stop("ERROR: normalization violation in sequence trace function.")}

  if( return_as_histlist_element )
  {
    result = list( "breaks"    = histlist_in[[ 1 ]]$breaks,
                   "density"   = (cumulative_hist/nseqs)/dx_arr,
                   "mids"      = histlist_in[[ 1 ]]$mids,
                   "xname"     = paste0(histlist_in$xname, "-->seq: ", sequence_in),
                   "equidist"  = histlist_in[[ 1 ]]$equidist
                   )
    
  }else{
  result =  (1/NORM)*cumulative_hist; 
  }
  
  return( result )
}

