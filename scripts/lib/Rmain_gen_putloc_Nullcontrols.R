putloc_dir="/fast/AG_Akalin/bosberg/nanopore/ref/Regions_of_interest/"

filenames = c("m6A_putlocs_Linder.rds",  
  "MihaM_2po_Ill.rds" )

putloc_path_list = paste0(putloc_dir, filenames )
baseshift = 50

for (putloc_path in putloc_path_list )
{
  putloc_set = readRDS( putloc_path )
  
  for ( groupname in names( putloc_set$Region_groups ) )
  {
   
   upstream_groupname   = paste0( groupname, "_NULLcontrol", as.character(baseshift),"bp_upstream" )
   downstream_groupname = paste0( groupname, "_NULLcontrol", as.character(baseshift),"bp_downstream" )
   
   temp_putlocs  = putloc_set$Region_groups[[groupname]]
   temp_downstream = temp_putlocs
   temp_upstream   = temp_putlocs
   
   multiplier = matrix ( 1, length(temp_putlocs), 1 )
   multiplier[ which (strand(temp_putlocs) == "-" )] = -1
   
   ranges( temp_downstream ) = IRanges( start = start( temp_downstream ) + (multiplier*baseshift ),
                                        end   = end( temp_downstream ) + (multiplier*baseshift ) )
   
   ranges( temp_upstream )   = IRanges( start = start( temp_upstream ) - (multiplier*baseshift ),
                                        end   = end(   temp_upstream ) - (multiplier*baseshift ) )

   putloc_set$Region_groups[[upstream_groupname]]   <- temp_upstream
   putloc_set$Region_groups[[downstream_groupname]] <- temp_downstream
      
   saveRDS( object = putloc_set, 
            file   = putloc_path );
  }
}