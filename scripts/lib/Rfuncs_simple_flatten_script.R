simple_flatten <- function( GR_in = stop("GR_in must be provided") )
  {
  strand = as.character( unique( strand( GR_in ) ))

  GR_in_splitby_startpos <- split( GR_in, start( GR_in ) )
  
  # model_kmer = lapply(
  #                     GR_in_splitby_startpos,
  #                        function(posn) 
  #                           unique(posn$model_kmer) 
  #                     )
  # 
  # for( posni in c(1:length(model_kmer)))
  #   {
  #   NNNspots = which( model_kmer[[posni]] =="NNNNN" )
  #   
  #   if( length( NNN_spots ) == length( model_kmer[[posni]] ) )
  #       { 
  #       model_kmer[[posni]] <- "ALLNs" } else{
  #       model_kmer[[posni]] <- model_kmer[[posni]][ -c(NNNspots) ]
  #       }
  #   
  #   if( length( model_kmer[[posni]] ) > 1 )
  #     { 
  #     stop("invalid name list") 
  #     }
  #   }
  
  result <- unlist( lapply( GR_in_splitby_startpos, 
                        function(posn)
                            sum(posn$event_mean * posn$event_length)/sum(posn$event_length) 
                                  )
                    )
  
  if( strand =="-" )
    {
    result    <- rev( result )
#    model_kmer <- rev( model_kmer )
    }  
  
#  model_kmer = unlist( model_kmer )
  return( result )
  
  }