# ===============================================
# checking if some really basic properties of DNA sequence might actually determine current flow.
# e.g. molecular weight, # of atoms,  electronegativity, stiffness, etc.


Base_properties=list( "A" = list("weight" = as.numeric( 135.13),    "Natoms" = as.numeric(15), "wpa" = 9.008667, "EA_theory"= -2.54, "EA_exp"= -0.74 ),
                      "C" = list("weight" = as.numeric( 111.1),     "Natoms" = as.numeric(13), "wpa" = 8.546154, "EA_theory"= -1.97, "EA_exp"= -0.40 ),
                      "G" = list("weight" = as.numeric( 151.13),    "Natoms" = as.numeric(16), "wpa" = 9.445625, "EA_theory"= -2.82, "EA_exp"= -1.23 ),
                      "T" = list("weight" = as.numeric( 126.115),   "Natoms" = as.numeric(15), "wpa" = 8.407667, "EA_theory"= -1.85, "EA_exp"= -0.32 ),
                      "U" = list("weight" = as.numeric( 112.08676), "Natoms" = as.numeric(12), "wpa" = 9.340563, "EA_theory"= -1.77, "EA_exp"= -0.11)
                    )    

poremodel_RNA=read.csv( "/clusterhome/bosberg/projects/nanopiper/dev/ref/poremodel_RNA.csv", 
                    sep="\t")

# get the whole list of 1024 sequences
seqlist = row.names( poremodel_RNA) 

# ==================================================================
# --- calculate weight of sequence.
calc_seq_property <- function( sequence_in        = stop("sequence     must be provided"),
                               Base_properties_in = stop("base properties list must be provided"),
                               property           = stop("property of interest must be specified.")
                            )
{
  # how long is the sequence? 
  k=nchar( sequence_in )
  
  # break it up into individual characters
  seq_in_singlechars = unlist( strsplit(sequence_in, split="") )
  
  # Check how many map to explicit and trace bases:
  Charmap = match( seq_in_singlechars, unlist( names( Base_properties_in) ) )
  
  # check that all base positions are recognized and accounted for.
  if ( length( which(!is.na(Charmap) ) ) != k )
    { stop("sequence contains unrecognized bases") }
  
  result = 0
  
  for ( i in 1:k )
    {
    result = result + Base_properties_in[[ Charmap[i] ]][[ property ]] 
    if( is.null(  Base_properties_in[[ Charmap[i] ]][[ property ]] ))
      { stop("Null property submitted to base set.") }
    }
return( result/k )
}

# ===============================================================
# --- return list of mean dwell times by kmer
get_passage_time_statistics <- function( GRL_kmers_in    = stop("kmer-split event GRL obj must be provided"),
                                         moment_in       = 1
                                   )
{
  
 result = unlist(  lapply( GRL_kmers_in, function(x)   passage_time_statistic ( GR_kmer_in      = x, 
                                                                                moment  = moment_in
                                                                               )
                           )
                   )
 
 return(result)
}
# --- return mean dwell time for this specific kmer
get_passage_time_statistic <- function( GR_kmer_in    = stop("GR obj for specific kmer must be provided"),
                                        moment        = 1
                                   )
{
 if( moment == 1 ){ # 1rst moment == mean
 result = mean( GR_kmer_in$event_length )
 }else if (moment == 2) { # end moment == std. dev
 result = sd( GR_kmer_in$event_length )
 }
 
 return(result)
}
  