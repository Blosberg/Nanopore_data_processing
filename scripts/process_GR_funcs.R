plot_compare_motif <-  function( seq  = stop("seq  must be provided"),
                                 Gro1 = stop("Gro1 must be provided"), 
                                 Gro2 = stop("Gro2 must be provided"),
                                 mincurrent = 50,
                                 maxcurrent = 150
                                 )
{ 
temp1 = Gro1[  Gro1$event_mean  > mincurrent ]
temp1 = temp1[ temp1$event_mean < maxcurrent ]

temp2 = Gro2[  Gro2$event_mean  > mincurrent ]
temp2 = temp2[ temp2$event_mean < maxcurrent ]

temp1_subset= temp1[ temp1$reference_kmer == seq ]
temp2_subset= temp2[ temp2$reference_kmer == seq ]

# current_window = c( mincurrent, maxcurrent )
par( mfrow = c(2,1) )

hist( temp1_subset$event_mean, 
      breaks = seq(from = mincurrent, to = maxcurrent, by = 1 ) , 
      xlab = "current", 
      main = m6A_motif )
hist( temp2_subset$event_mean, 
      breaks = seq(from = mincurrent, to = maxcurrent, by = 1 ) , 
      xlab = "current", 
      main = m6A_motif )

return(0)

}