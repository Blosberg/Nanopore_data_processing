
#=========================================================================================
plot_sequence_current <-  function( seq,
                                    poremodel_table = stop("table must be provided"),
                                    k=5,
                                    add= TRUE) 
{ 
  N=nchar(seq);
  
  xvals   = c(1:(N-k+1) )
  yvals   = poremodel_table[ substr(seq, 1, k ) ,1]
  
  for (i in c(1:N-k+1))
  {
    yvals[i] =  poremodel_table[ substr(seq, i, (i+k-1) ) ,1]
  }
  
  return(yvals)
}
#=========================================================================================

poremodel_statconsts[  substr(rownames(poremodel_statconsts),1,3)=="GCC", ] 

BC7 = "ACAAGACAGCAC"

a = 70; b=80
temp    <- poremodel_statconsts[  poremodel_statconsts[,1] > a, ] 
window1 <- rownames( temp[ temp[,1] < b , ] )
tail1   <- substr(window1,2,5)

a = 110; b=120
temp     <- poremodel_statconsts[  poremodel_statconsts[,1] > a, ] 
window2  <- rownames(temp[ temp[,1] < b , ] )
head2    <- substr(window2,1,4)

window2_headfilter       <- window2[ !is.na( match(head2, tail1)) ]
window2_headfilter_tail  <- substr(window2_headfilter, 4,5)

a = 110; b=120
temp     <- poremodel_statconsts[  poremodel_statconsts[,1] > a, ] 
window5  <- rownames(temp[ temp[,1] < b , ] )
head5    <- substr(window5,1,2)

window5_update = window5[ !is.na(match(head5, window2_headfilter_tail))  ]
window5_update_tail = substr(window5_update,4,5)

a = 110; b=115
temp     <- poremodel_statconsts[  poremodel_statconsts[,1] > a, ] 
window8  <- rownames(temp[ temp[,1] < b , ] )
head8    <- substr(window8,1,2)

window8_update <- window8[ !is.na(match(head8,window5_update_tail )) ]


poremodel_statconsts[ window8[ is.na( match(head8, window5_update_tail ) )] ,]


ACAAGAGCTGAT

#=====================
AVOID: GGATCC, and ATGCAT


BC1 = "CCCCGATTAGAT"
BC2 = "AGATGCGACCCC"
BC3 = "AGCCCTGATGCC"
BC4 = "CTGCCGACATTA"
BC5 = "CAAGCCCTTGAC"
BC6 = "ACGGATCGCCAG"
BC7 = "ACAAGAGCTGAT"

B1 =  plot_sequence_current( seq= BC1, poremodel_table = poremodel_statconsts )
B2 =  plot_sequence_current( seq= BC2, poremodel_table = poremodel_statconsts)
B3 =  plot_sequence_current( seq= BC3, poremodel_table = poremodel_statconsts)
B4 =  plot_sequence_current( seq= BC4, poremodel_table = poremodel_statconsts )
B5 =  plot_sequence_current( seq= BC5, poremodel_table = poremodel_statconsts )
B6 =  plot_sequence_current( seq= BC6, poremodel_table = poremodel_statconsts )
B7 =  plot_sequence_current( seq= BC7, poremodel_table = poremodel_statconsts )


plot( xvals, B1, type = "o", col = "red", xlab = "position", ylab = "current[pA]", main="barcode current levels")
lines( xvals, B2, type = "o", pch=2, col = "blue")
lines( xvals, B3, type = "o", pch=3, col = "green")
lines( xvals, B4, type = "o", pch=4, col = "black")
lines( xvals, B5, type = "o", pch=5, col = "pink")
lines( xvals, B6, type = "o", pch=6, col = "orange")
lines( xvals, B7, type = "o", pch=7, col = "brown")


