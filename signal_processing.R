# import nanopore data and run some statistics on them:
# rm(list=ls()) # CLEAN UP EVERYTHING

Event_folder = "/scratch/AG_Akalin/bosberg/nanopore/pipeline_output/20180417_1233_polyA_RNA/04_eventalign/"


Event_file_all   = "all_first10M.csv"
fin_all = file.path(Event_folder, Event_file_all)
dat_first10M  = read.csv( file = fin_all, 
                      sep = '\t', 
                      stringsAsFactors=FALSE, 
                      header = TRUE,
                      nrow= N
)
test_first10M  = dat_first10M$event_level_mean + 3.14
print(paste("finished with part ", as.character(i) ) )


dat_test_parts=list()
i=0
prefix_string="Ealign_splitby_10M_a_1M_c_100K_f_10K_i_1K_g_100_"

for ( char in c('a', 'b', 'c', 'd', 'e', 'f', 'g', 'h', 'i', 'j') )
  {
  
  print(  paste( " -------------------------- Beginning character ", char , " -----------------------------" ) )
  
  i=i+1
  Event_file           = paste0( prefix_string, char )
  fin                  = file.path( Event_folder, Event_file )
  print( paste("------ READING IN FILE: ", fin ) )
  
  if( i== 1)
  {should_be_header = TRUE} else{ should_be_header = FALSE }
  
  dat_test_parts[[i]]  = read.csv( file = fin, 
                                   sep = '\t', 
                                   stringsAsFactors=FALSE, 
                                   header = should_be_header
  )
  names( dat_test_parts[[i]] ) <- names( dat )
  
  print( dat_test_parts[[i]]$event_level_mean[1:10]  )
  print( dat_test_parts[[i]]$event_level_mean[1:10] + 1.01  )
  
  print(  paste("finished with part ", char ) )
}



# ========================================================================================

fin_test_inf        = file.path(Event_folder, "test_inf.csv" )
fin_test_allfinite = file.path(Event_folder, "test_allfinite.csv" )

dat_allfinite  = read.csv( file = fin_test_allfinite, 
                            sep = '\t', 
                            stringsAsFactors=FALSE, 
                            header = TRUE
                           )

dat_inf  = read.csv( file=fin_test_inf, 
                     sep='\t', 
                     stringsAsFactors=FALSE, 
                     header = TRUE
)


mincurrent=39
maxcurrent=171
current_window = c( mincurrent, maxcurrent )

dat  = read.csv( file=fin, 
                 sep='\t', 
                 stringsAsFactors=FALSE, 
                 header = TRUE
                 )

# hist(dat$model_mean)

temp           = dat[  dat$event_level_mean  > mincurrent, ]
dat_windowed   = temp[ temp$event_level_mean < maxcurrent, ] 
dat_win_finite = dat_windowed [ which (! is.na (dat_windowed$event_level_mean) ), ]
rm(temp)

mid_groups = take_midbase_averages( Datin = dat )

# ============================================================

freq_group_A = group_by_freq ( Datin           = dat_win_finite, 
                               seq             = "A",
                               xrange          = current_window,
                               break_set       = seq(from = mincurrent, to = maxcurrent, by = 2 )
)

m6A_motif_subset = dat_win_finite[ which( dat_win_finite$reference_kmer == "GGACT") ,]
hist( m6A_motif_subset$event_level_mean, main="m6A subset", xlab="current mean", xlim=xrange, ylab="", breaks=50 )

