# ==================================================================
# --- get list (separated by position of modification) of differences
# --- between current with modified base vs. control. Each list entry
# --- tabulates

get_moddiff_vs_modpos_seq <- function ( GRmod      = stop("GRmod must be provided"),
                                        GRcontrol  = stop("GRcontrol must be provided"),
                                        k=5,
                                        sigthresh = 5 )
{
  # --- temporary kludge fix for strand alignment. Figure out a more efficient way to get both strands.
  GRcontrol <- GRcontrol[strand(GRcontrol)=="+"]

  SEQHITS = sort( unique(GRmod$reference_kmer) )
  Nunique_sequences = length( SEQHITS )

  i=1
  result = matrix(0,Nunique_sequences, 3  )
  row.names(result) <- SEQHITS

  GRmodL_modpos_split <- split(GRmod, GRmod$modposition)
  GRcontrolL_seqsplit <- split(GRcontrol, GRcontrol$reference_kmer)

  difflist = list()

  for ( i in c (1:k))
    {
    temp = split( GRmodL_modpos_split[[i]], GRmodL_modpos_split[[i]]$reference_kmer)

    m=1;n=1;
    temp2=GRangesList()
    for ( m in c (1:length(temp)))
      {
      if ( length( temp[[m]] ) > sigthresh )
        {
        temp2[[n]] <- temp[[m]]
        names(temp2)[n]<- names(temp)[m]
        n <- (n+1)
        }
      }


    difflist[[i]]   <- lapply( c(1:length(temp2)), function(x)
                               ( (mean(temp2[[x]]$event_mean) -
                                    mean(GRcontrolL_seqsplit[[ names(temp2)[x] ]]$event_mean))
     #                            /(pore_model_list[[ names(temp2)[x]  ]]$std_dev)  )
                                )
                              )
  }

  return(difflist)
}
