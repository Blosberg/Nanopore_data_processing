# prepare some plotting scripts:
source('scripts/lib/plot_current_distributions_funcs.R', echo=TRUE)
M=100

# import reference:
# read in reference data of gene locations
genelist_gtf = readRDS("/data/akalin/Base/Annotation/hg19/EnsemblGenes/170317_Homo_sapiens.GRCh37.87.chr.gtf.granges.rds")
# split these up by gene name:
genelist_splitbyname_GRL = split( genelist_gtf, 
                                  genelist_gtf$gene_name )

# =============== Now examing experimental data of gene life-times: ==============
Gene_lifetime_measurement = "/clusterhome/bosberg/projects/nanopiper/dev/ref/IV_incorporation_MM_calcs/readcounts_slopes.txt"

readcount_slopes = read.csv( file   = Gene_lifetime_measurement, 
                             row.names = 1,
                             header = FALSE,
                             sep    = "\t"
                            )

# col names provided in email from Emanuel
names( readcount_slopes ) <- c( "Gene_length",
                                "Readcounts_2017",
                                "Readcounts_2018",
                                "Names_again",
                                "Stab_1",
                                "Stab_2" )
readcount_slopes <- readcount_slopes [,c(1:3,5:6) ]

# get number of genes overall:
N = dim( readcount_slopes)[1]


# ====== now select based on lifetime: 

# Take the first "M" from either end.
Thresh1_short = sort( readcount_slopes$Stab_1, decreasing = FALSE )[M]
Thresh2_short = sort( readcount_slopes$Stab_2, decreasing = FALSE )[M]

list_shortlived_1 = which ( readcount_slopes[,"Stab_1"] <= Thresh1_short )
list_shortlived_2 = which ( readcount_slopes[,"Stab_2"] <= Thresh2_short )

matchlist_short =  match( list_shortlived_1,
                          list_shortlived_2 )
matchlist_short = matchlist_short[ !is.na(matchlist_short) ]

list_shortlived_both = list_shortlived_2[ matchlist_short ]

shortlived_both <- readcount_slopes[ list_shortlived_both, ]

# ---- now get the M longest-lived genes:

Thresh1_long = sort( readcount_slopes$Stab_1, decreasing = TRUE )[M]
Thresh2_long = sort( readcount_slopes$Stab_2, decreasing = TRUE )[M]

list_longlived_1 = which ( readcount_slopes[,"Stab_1"] >= Thresh1_long )
list_longlived_2 = which ( readcount_slopes[,"Stab_2"] >= Thresh2_long )

matchlist_long =  match( list_longlived_1,
                    list_longlived_2 )

matchlist_long =  matchlist_long[ !is.na(matchlist_long) ]
list_longlived_both = list_longlived_2[ matchlist_long ]

longlived_both <- readcount_slopes[ list_longlived_both, ]

# --------------  now do some plotting --------------------------

par(mfrow=c(1,1))
plot( readcount_slopes[, "Stab_1"], ylab = "Stab_1", main = "Transcript stability, 2017" )
#points( shortlived_1[, "Stab_1" ], col="red" )

points( list_shortlived_both, 
        shortlived_both$Stab_1,
        lwd = 3,
        col="red" )

points( list_longlived_both, 
        longlived_both$Stab_1, 
        lwd = 3,
        col="blue" )

legend( 200, -0.2,
        legend = c( "Transcripts", "Low stability", "High stability"),
        col=c("black", "red", "blue", lty=1:2, cex=0.8),
        lwd=c(1,3,3),
        pch=c(1,1,1)
        )

# ==================================================
assembly="hg19"

# Now find the locations of these genes on our reference assembly.

# find matching names:
found_list_shortlived_loci  <- match( row.names(shortlived_both),
                                      names(genelist_splitbyname_GRL))

found_list_longlived_loci  <- match( row.names(longlived_both),
                                     names(genelist_splitbyname_GRL))

# remove NA
found_list_shortlived_loci <- found_list_shortlived_loci[ !is.na(found_list_shortlived_loci) ]
found_list_longlived_loci  <- found_list_longlived_loci[ !is.na(found_list_longlived_loci) ]

# build an actual GRL of these loci
shortlived_GRL <-   genelist_splitbyname_GRL[ found_list_shortlived_loci ] 
longlived_GRL  <-   genelist_splitbyname_GRL[ found_list_longlived_loci ] 

# build a single structure containing a list of various regions of interest grouped
RsoI = list( "Region_groups" = list( "shortlived"    = shortlived_GRL,
                                     "longlived"     = longlived_GRL),
             "assembly"      = assembly,
             "Num_g_select"  = M )

saveRDS( RsoI, file = paste0("/scratch/AG_Akalin/bosberg/Stability/stabilityseg_hg19_M",as.character(M),".rds") )
