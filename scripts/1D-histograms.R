





current_histmin=100 
current_histmax=140
current_histres=0.5
breakset = seq( current_histmin, current_histmax, current_histres)

GGACT_3shifted_middle_olap_splitby_read = split( GGACT_3shifted_middle_olap, GGACT_3shifted_middle_olap$read_index )


plot_hist_with_edge_trimming( GGACT_from_all_events$event_mean, 
                              breaks_in = breakset,
                              main = "GGACT-seq CITS putmod against bulk"
                              )
plot_hist_with_edge_trimming( GGACT_3shifted_middle_olap$event_mean, 
                              breaks_in = breakset,
                              add = T,
                              col = rgb( 1, 0, 0, 0.5 )
                              )

length( GGACT_from_all_events)
20*c(1:6000)
