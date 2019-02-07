# saveRDS( readlist_on_chr19, file = "/scratch/AG_Akalin/bosberg/nanopore/pipeline_output/20180913_1457_invitro_rps16_unmod/readlist_on_chr19.rds")
# readlist_on_chr19 <- readRDS( "/scratch/AG_Akalin/bosberg/nanopore/pipeline_output/20180913_1457_invitro_rps16_unmod/readlist_on_chr19.rds")
# saveRDS( IVum_c19_splitby_read_GRL, file="/scratch/AG_Akalin/bosberg/nanopore/pipeline_output/20180913_1457_invitro_rps16_unmod/06_GRobjects/reads_c19_GRL.rds")
# IVum_c19_splitby_read_GRL <- readrds("/scratch/AG_Akalin/bosberg/nanopore/pipeline_output/20180913_1457_invitro_rps16_unmod/06_GRobjects/reads_c19_GRL.rds")

rps_plasmid = toupper("ctcccgtcgcgacgcagtgctcaaggcgcctgcgcagaccctgaaaagcggccagggtggcccctagctttccttttccggttgcggcgccgcgcggtgaggttgtctagtccacgctcggagccatgccgtccaagggcccgctgcagtctgtgcaggtcttcggacgcaagaagacagcgacagctgtggcgcactgcaaacgcggcaatggtctcatcaaggtgaacgggcggcccctggagatgattgagccgcgcacgctacagtacaagctgctggagccagttctgcttctcggcaaggagcgatttgctggtgtagacatccgtgtccgtgtaaagggtggtggtcacgtggcccagatttatgctatccgtcagtccatctccaaagccctggtggcctattaccagaaatatgtggatgaggcttccaagaaggagatcaaagacatcctcatccagtatgaccggaccctgctggtagctgaccctcgtcgctgcgagtccaaaaagtttggaggccctggtgcccgcgctcgctaccagaaatcctaccgataagcccatcgtgactcaaaactcacttgtataataaacagtttttgagggattttaaagtttcaagaactgtgtgtggccttatgtgttggtattggatgtttaaccagacagaaatcagtaaacatcctggacctataaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa")
N_rps = nchar(rps_plasmid)
N_comp <- nchar(comp_sequence)

transcript_NM_002896=toupper("ttcaaggcaaacgaatgcacgtgcagttgtccaccagccggcttaggactgcgcccgggatgggagaccagagcggctgctatcggtgcgggaaagaggggcactggtccaaagagtgtccgatagatcgttcaggccgcgtggcagacttgaccgagcaatataatgagcaatacggagcagtgcgtacgccttacaccatgagctatggggattcattgtattacaacaacgcgtacggagcgctcgatgcctactacaagcgctgccgtgctgcccggtcctatgaggcagtggcagctgcagctgcctccgtgtataattacgcagagcagaccctgtcccagctgccacaagtccagaatacagccatggccagtcacctcacctccacctctctcgatccctacgatagacacctgttgccgacctcaggagctgctgccacagctgctgctgcagcagcagccgctgctgctgttactgcagcttccacttcatattacgggcgggatcggagccccctgcgtcgcgctacagccccagtccccactgttggagagggctacggttacgggcatgagagtgagttgtcccaagcttcagcagccgcgcggaattctctgtacgacatggcccggtatgagcgggagcagtatgccgatcgggcgcggtactcagccttttaaagcttga")
N_NM00=nchar(transcript_NM_002896)


# Restrictions (to be avoided): 
GGATCC
ATGCAT
# http://jgeisler0303.github.io/deBruijnDecode/

dB4 <- "CCCCAGACGAGCACAACTGGGCGTAAGGCCTATACTCAAGAACACGTCGCCCTTCGAATGCCGTTTTCACTACATCTCCTGAAATAGCCAATTGCGCTGTCCACCTCTAGTATTTGGTACCGATTAGGGGTCTGCTCGTGTTGACCCGCATTCCGGATATGAGTTCTTGTGCAGCGGGAAGCTTACGGCAAAAGTGATGTAGAGGACTTTATCGGTTAACGCGACAGGTGGAGATCAGTCATAAACCATGGCTAATCCC"
dB5 <- "CCCCCATCCGTAAGGTCCTATAATGTCGCCCTCTCCTTCATCAAGGCTTAAAACTTAGTGCCCGACTGGTTGTTCTTGGGATCAGAGCCCAGACGGCCTACCCAAACCTCGAGACTATGTTTAGACATATACCTAACCGATGATCCCCGGTGAAATAGAGATCGACCTGCTAAGCCGCGTCGGACAATAAGTCGATCTACATCGCAACTGTGGCGACAGGCCGTTCGCTATCGTGGTACCACAGTTTCTATTGTCTTTTTCGTTAGCCAAGTTGCGCCGAGTATTTGTGTTGGTCGTCCAATGGCCCCTGGCAGATTTACTGAAGGAGGTGGAAACGTCACATTTCACTAATTGCTTTCCCTTGTACACTGCGAGGGCACTCGCGCTTGCCTTACGACGTTTTAAGACCCGCACACGTACTTGATGCCATTATAGTAACTCCGCCAGTACGTGCTCGTATCCTGAGAATACTCAATTAGGTATGCGTAGAAGCTTCTCAGGGTGCGGAGTGGGGGTAGGCATTCCAGGATATTACAGCGGGAGAGGAAGTAGTTCCGAATTTTGAATCACGGGCTGACCAACAAGAACCCTAGTCTCGGGTCTGCCGGATGGGTTTGGAGCTCTGATACAAATCTGTATATGGATTGGCTACTAGGGGCCACTTTATCTCTAGCTGGACGATTCGGCGCAGCCTGTTATTCAGCTAGATGTAGCGTTGACACCGCTGTCCGGCAAGCGATAGCACCATGGTGTCAAAGATAACAGAAAAAGCAATCGGTTCAACGCGGCTCATAAACACAACCAGCAGGTTACCGGGGACCGTCTAAATGACTTCCTCCACCTTTGCAGTGTAAAGTGATTAACGAAAGGGAACATGTGACGCTCCCGTGAGCGCGAACGGAATGCTGCAAAATTCTGGGCGTGTGCATGAACTACGGTAATATCATTGAGGCGGTCAGTCATGCACGAGCATAGGACTCTTATGAGTTAATCCATACGCATCTTCGAAGAGTCCCACGCCTCACCCC"


length_check       <- lapply( IVum_c19_splitby_read_GRL, 
                           function(x) length( unique ( start(x) ) ) 
                          )
length_check_array = unlist(length_check)

hist(length_check_array, 200) 
IV_readlengths_hist <- hist(length_check_array, 200, plot = FALSE) 

find_maxfreq_length( histdat = IV_readlengths_hist, a=400, b=600 )

a=95; b = 105;
peak_subset   <- which( length_check_array ==   495 )
sequence      <- extract_sequence_from_read( IVum_c19_splitby_read_GRL[[ peak_subset[2] ]] )

sequences_peak1 <- unlist( lapply( IVum_c19_splitby_read_GRL[subset_2], 
                            function(x) extract_sequence_from_read( x ) 
                           ) ) 
# peak lengths occure at:  101,  147,  261,  495

kmer_divers       <- unlist( lapply( IVum_c19_splitby_read_GRL, 
                                    function(x) length( unique (x$model_kmer ) ) 
                                    ))

kmer_divers_array = unlist( kmer_divers )

IVum_c19_splitby_read_GRL
# -----------------------------

k5merset = list()
k4merset = list()

k5mer_compset  = list()
k4mer_compset  = list()

k5mer_NMset  = list()
k4mer_NMset  = list()


for ( p in c(1:N_rps-4))
  {
  k5merset[p]   <- substr(rps_plasmid, p, p+4 )
  }
for ( p in c(1:N_rps-3))
  {
  k4merset[p]      <- substr(rps_plasmid, p, p+3 ) 
 }

for ( p in c(1:N_NM00-4))
  {
  k5mer_NMset[p] <- substr(transcript_NM_002896, p, p+4 )
  }
for ( p in c(1:N_NM00-3))
  {
  k4mer_NMset[p] <- substr(transcript_NM_002896, p, p+3 ) 
 }


rpsplasmid_k5merset_array <- unique( unlist( k5merset ) )
rpsplasmid_k4merset_array <- unique( unlist( k4merset ) )

comp_k5merset_array <- unique( unlist( k5mer_compset ) )
comp_k4merset_array <- unique( unlist( k4mer_compset ) )

NM_k5merset_array <- unique( unlist( k5mer_NMset ) )
NM_k4merset_array <- unique( unlist( k4mer_NMset ) )


k5mer_diversity = length( unique(k5merset_array))
k4mer_diversity = length( unique(k4merset_array))

k4mers_in_rpsplasmid <- unique(k4merset_array)
