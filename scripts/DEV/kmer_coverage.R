# IVunmod_dataset = "/scratch/AG_Akalin/bosberg/nanopore/pipeline_output/20180913_1457_invitro_rps16_unmod/06_GRobjects/invitro_rps16_unmod_reads_GR.rds"
# IVunmod_RDSdat <- readRDS(IVunmod_dataset)
# kmers_in_rps16 <- unique( IVunmod_RDSdat$Events_all_GR$model_kmer )

deBruijn_4="CCCCAGACGAGCACAACTGGGCGTAAGGCCTATACTCAAGAACACGTCGCCCTTCGAATGCCGTTTTCACTACATCTCCTGAAATAGCCAATTGCGCTGTCCACCTCTAGTATTTGGTACCGATTAGGGGTCTGCTCGTGTTGACCCGCATTCCGGATATGAGTTCTTGTGCAGCGGGAAGCTTACGGCAAAAGTGATGTAGAGGACTTTATCGGTTAACGCGACAGGTGGAGATCAGTCATAAACCATGGCTAATCCC"

deBruijn_5="CCCCCATCCGTAAGGTCCTATAATGTCGCCCTCTCCTTCATCAAGGCTTAAAACTTAGTGCCCGACTGGTTGTTCTTGGGATCAGAGCCCAGACGGCCTACCCAAACCTCGAGACTATGTTTAGACATATACCTAACCGATGATCCCCGGTGAAATAGAGATCGACCTGCTAAGCCGCGTCGGACAATAAGTCGATCTACATCGCAACTGTGGCGACAGGCCGTTCGCTATCGTGGTACCACAGTTTCTATTGTCTTTTTCGTTAGCCAAGTTGCGCCGAGTATTTGTGTTGGTCGTCCAATGGCCCCTGGCAGATTTACTGAAGGAGGTGGAAACGTCACATTTCACTAATTGCTTTCCCTTGTACACTGCGAGGGCACTCGCGCTTGCCTTACGACGTTTTAAGACCCGCACACGTACTTGATGCCATTATAGTAACTCCGCCAGTACGTGCTCGTATCCTGAGAATACTCAATTAGGTATGCGTAGAAGCTTCTCAGGGTGCGGAGTGGGGGTAGGCATTCCAGGATATTACAGCGGGAGAGGAAGTAGTTCCGAATTTTGAATCACGGGCTGACCAACAAGAACCCTAGTCTCGGGTCTGCCGGATGGGTTTGGAGCTCTGATACAAATCTGTATATGGATTGGCTACTAGGGGCCACTTTATCTCTAGCTGGACGATTCGGCGCAGCCTGTTATTCAGCTAGATGTAGCGTTGACACCGCTGTCCGGCAAGCGATAGCACCATGGTGTCAAAGATAACAGAAAAAGCAATCGGTTCAACGCGGCTCATAAACACAACCAGCAGGTTACCGGGGACCGTCTAAATGACTTCCTCCACCTTTGCAGTGTAAAGTGATTAACGAAAGGGAACATGTGACGCTCCCGTGAGCGCGAACGGAATGCTGCAAAATTCTGGGCGTGTGCATGAACTACGGTAATATCATTGAGGCGGTCAGTCATGCACGAGCATAGGACTCTTATGAGTTAATCCATACGCATCTTCGAAGAGTCCCACGCCTCACCCC"

HEK_dataset = "/scratch/AG_Akalin/bosberg/nanopore/pipeline_output/20180417_1233_HEK293_polyA_RNA/GRdat_allevents.RDS"
HEK_RDSdat <- readRDS(HEK_dataset)


poremodel_statconsts = read.csv( "scripts/ref/pore_model_table.csv",
                                  stringsAsFactors = F,
                                  header           = TRUE,
                                  row.names        = 1,
                                  sep= "\t"
                                 )


PATH_putloc="/scratch/AG_Akalin/bosberg/nanopore/ref/Linder_2015_nature_supp/put_m6A_locs_CITS.csv"
CITS_put <- get_putlocs( PATH_putloc )
# ===========================================

all_k5mers <- row.names( poremodel_statconsts )
all_k4mers <- unique( substr( all_k5mers, 1,4) ) 

# original generation of this variable commented out above (saved in .rds file for easy retrieval) :
# saveRDS( rpsplasmid_k5merset_array, file = "/scratch/AG_Akalin/bosberg/nanopore/pipeline_output/20180913_1457_invitro_rps16_unmod/rps16_plasmid_k5mers.rds")

rps16_plasmid_k5mers  <- readRDS("/scratch/AG_Akalin/bosberg/nanopore/pipeline_output/20180913_1457_invitro_rps16_unmod/rps16_plasmid_k5mers.rds")
rps16_plasmid_k5mers <- kmers_in_rps16[ kmers_in_rps16 != "NNNNN" ]

kmers_not_in_rps16_plasmid  <- all_kmers[  which( is.na ( match( all_kmers,  rps16_plasmid_k5mers ) )) ]
k4mers_not_in_rps16_plasmid <- all_k4mers[ which( is.na ( match( all_k4mers, k4mers_in_rpsplasmid ) )) ]

# ===========================================
length_THRESH = 700


length_check       <- lapply( HEK_RDSdat$allevents_GRL_splitbyread, 
                              function(x) length( unique ( start(x) ) ) 
                             )
length_check_array <- unlist(length_check) 

HEK_short_reads_GRL <- HEK_RDSdat$allevents_GRL_splitbyread[ length_check_array < length_THRESH ]

#--------------------------------------------
# with putmod overlaps: 

HEK_short_reads_putmod_olaps                 <- findOverlaps( HEK_short_reads_GRL, CITS_put )
HEK_short_reads_AND_putmod_olaps_GRL         <- HEK_short_reads_GRL[ unique( queryHits( HEK_short_reads_putmod_olaps ) )  ]
HEK_short_reads_with_putmod_olaps_k5mercompl <- unlist( lapply( HEK_short_reads_AND_putmod_olaps_GRL, 
                                                       function(x) get_targetcov ( read_GR    = x,
                                                                                   target_kmers = k5mers_not_in_rps16
                                                                                  ) 
                                                      ))

HEK_short_reads_with_putmod_olaps_k4mercompl <- unlist( lapply( HEK_short_reads_AND_putmod_olaps_GRL, 
                                                       function(x) get_targetcov ( read_GR      = x,
                                                                                   target_kmers = k4mers_not_in_rps16_plasmid
                                                                                  ) 
                                                      ))


read_indices_over_200_k5mer_comp   <- HEK_short_reads_with_putmod_olaps_k5mercompl[ HEK_short_reads_with_putmod_olaps_k5mercompl > 200 ] 

hist(HEK_short_reads_with_putmod_olaps_k4mercompl[ names(read_indices_over_200_k5mer_comp)])

which.max(HEK_short_reads_with_putmod_olaps_k4mercompl[ names(read_indices_over_200_k5mer_comp)])

# and the best result seems to be "131710"
comp_sequence <- extract_sequence_from_read(  HEK_short_reads_AND_putmod_olaps_GRL[["131710"]] )


read_indices_satisfying_k4mer_comp <- HEK_short_reads_with_putmod_olaps_k4mercompl[ HEK_short_reads_with_putmod_olaps_k4mercompl >=25 ]

#--------------------------------------------
# with putmod overlaps: 
HEK_short_reads_GRL

HEK_short_reads_kmercompl <- unlist( lapply( HEK_short_reads_GRL, 
                                             function(x) get_targetcov ( read_GR    = x,
                                                                         target_kmers = kmers_not_in_rps16
                                                                        ) 
                                          ))

#--------------------------------------------

reads_with_putmod_olaps_kmercompl_array <- unlist(reads_with_putmod_olaps_kmercompl)

optimal_i <- which.min(unlist(length_check))

reduced_optimal <- reduce( HEK_reads_with_putmod_olaps_and_kmercompl_over_THRESH[optimal_i])
#####################

