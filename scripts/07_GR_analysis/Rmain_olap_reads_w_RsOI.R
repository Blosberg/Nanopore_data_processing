# Compare observed reads with known reference locations of interest.
## Help section
if("--help" %in% args) {
  cat("
      Render to report

      Arguments:
      Rfuncs_olap_w_RsOI  -- script with function definitions used here.
      reads_GRL_in        -- read_object structure in GRL format.
      output_Olap         -- rds filename read-separated GRL output corresponding to overlap with each of the RsOI.
      samplename          -- name of sample for documentation
      RsOI_files          -- rds files corrsponding to a set of RsOI to compare the reads against.
                             There should be only one per execution of this script (the script should be called 
                             once for every ROI).
      logFile             -- filename to pipe output to

      Example:
      ./test.R --arg1=1 --arg2='output.txt' --arg3=TRUE \n\n")

  q(save="no")
}

## Parse arguments (we expect the form --arg=value)
parseArgs <- function(x) strsplit(sub("^--", "", x), "=")
argsDF    <- as.data.frame(do.call("rbind", parseArgs(args)))
argsL     <- as.list(as.character(argsDF$V2))

names(argsL) <- argsDF$V1


##########################################

# import the nanopore reads as GRL
Reads_raw  <- readRDS ( argsL$ )
reads_untreated <- readRDS("/scratch/AG_Akalin/bosberg/nanopore/pipeline_output/20180417_1233_HEK293_polyA_RNA/06_GRobjects/HEK293_polyA_reads_GRL.rds")
reads_4SU_IAA   <- readRDS("/scratch/AG_Akalin/bosberg/nanopore/pipeline_output/20190130_1420_293_4Su_IAA/06_GRobjects/6c2fe901_reads_GRL.rds")

# find overlaps with the genes we've identified as short-lived"
olaps_untreated <- findOverlaps( reads_untreated$Events_GRL_splitbyread, shortlived_alignments_GRL )
olaps_4SU_IAA   <- findOverlaps( reads_4SU_IAA$Events_GRL_splitbyread,   shortlived_alignments_GRL )

# select only those transcripts.
shortlived_transcripts_untreated_GRL <- reads_untreated$Events_GRL_splitbyread[ queryHits( olaps_untreated) ]
shortlived_transcripts_4SU_IAA_GRL   <- reads_4SU_IAA$Events_GRL_splitbyread[   queryHits( olaps_4SU_IAA)   ]

#-----------------------------------------

shortlived_allevents_untreated_GR <- unlist(shortlived_transcripts_untreated_GRL ) 
shortlived_allevents_4SU_IAA_GR   <- unlist(shortlived_transcripts_4SU_IAA_GRL ) 

# now split up by kmer:
shortlived_kmers_untreated_GR <- split( shortlived_allevents_untreated_GR, shortlived_allevents_untreated_GR$model_kmer )
shortlived_kmers_untreated_GR <- shortlived_kmers_untreated_GR[ names(shortlived_kmers_untreated_GR) != "NNNNN" ]

shortlived_kmers_4SU_IAA_GR <- split( shortlived_allevents_4SU_IAA_GR, shortlived_allevents_4SU_IAA_GR$model_kmer )
shortlived_kmers_4SU_IAA_GR <- shortlived_kmers_4SU_IAA_GR[ names(shortlived_kmers_4SU_IAA_GR) != "NNNNN" ]

allseqs <- names(shortlived_kmers_untreated_GR ) 

# build_histlist
shortlived_untreated_histlist <- lapply( shortlived_kmers_untreated_GR, function(x)  get_fullhist_data  ( GR_input = x )  )

shortlived_4SU_IAA_histlist <- lapply( shortlived_kmers_4SU_IAA_GR, function(x)  get_fullhist_data  ( GR_input = x )  )


shortlived_kmerhistolap_4SU_WT     <- unlist( lapply( allseqs, function(seq)  calculate_histogram_overlap (  
                                            hist1    = shortlived_untreated_histlist[[seq]],
                                            hist2    = shortlived_4SU_IAA_histlist[[seq]]
                                            )
                  ) )

names(shortlived_kmerhistolap_4SU_WT) <- allseqs

# take the subset with largest difference:
subset <- shortlived_kmerhistolap_4SU_WT[ shortlived_kmerhistolap_4SU_WT< 0.75]

#
mid_T_subset <- shortlived_kmerhistolap_4SU_WT[ which( substr( names(shortlived_kmerhistolap_4SU_WT), 3,3 ) == "T" ) ]



i=0
i=i+1; seq = names(mid_T_subset)[i];  plot_histogram_overlap (   hist1    = shortlived_untreated_histlist[[seq]],
                                                                 hist2    = shortlived_4SU_IAA_histlist[[seq]],
                                                                 seq_in   = seq )

  
