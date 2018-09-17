# library(Biostrings)
# read.DNAStringSet
#rm(list=ls()) # CLEAN UP EVERYTHING
# library("Biostrings")
# hg19_canon_ref = readDNAStringSet("/scratch/AG_Akalin/refGenomes/hg19_canon/hg19_canon.fa")
# negstrand_reads <- reads_GR[ strand( reads_GR)=="-",]
# source("https://bioconductor.org/biocLite.R") #--- ignore warnings telling you to upgrade
# biocLite("GenomicRanges")

# ===== LOAD GRANGES OBJECTS PRODUCED BY SNAKEMAKE ======

# @@@ TODO: make these arguments come from the command-line
# e.g. [hard-coded for now]

DIR_GR="/home/bosberg/projects/nanopore/pipeline_out/20180417_1233_polyA_RNA/07_GRobjects/"
file_GR="ec9bd1f6_GR.RData"
load( file.path( DIR_GR, file_GR ) )

hg19_ref <- readDNAStringSet("/scratch/AG_Akalin/refGenomes/hg19_canon/hg19_canon.fa")

# ======= LOAD PUTATIVE MODIFICATION LOCATIONS: =========

put_path="/scratch/AG_Akalin/bosberg/nanopore/ref/Linder_2015_nature_supp/"
CITS_filename="put_m6A_locs_CITS.csv"
CIMS_filename="put_m6A_locs_CIMS.csv"

CITS_rawdat=read.csv( file.path( put_path, CITS_filename ), sep = "\t", header = TRUE )
CIMS_rawdat=read.csv( file.path( put_path, CIMS_filename ), sep = "\t", header = TRUE )

CIMS_put =  GRanges( seqnames = CIMS_rawdat$Chr ,
                     ranges   = IRanges ( start  = CIMS_rawdat$Start,
                                          end    = CIMS_rawdat$End ),
                     strand   = CIMS_rawdat$Strand
                    )
CITS_put =  GRanges( seqnames = CITS_rawdat$Chr ,
                     ranges   = IRanges ( start  = CITS_rawdat$Start,
                                          end    = CITS_rawdat$End ),
                     strand   = CITS_rawdat$Strand
                    )

CIMS_put = sort( CIMS_put[ which( seqnames(CIMS_put) != "chrM") ] ) # 9411 locations
CITS_put = sort( CITS_put[ which( seqnames(CITS_put) != "chrM") ] ) # 6543 locations

#  ================================================
# convert GRanges object into something with/without mod.

Nreads   = length(reads_GR)

reads_cims_putmod_pairs   = findOverlaps( reads_GR, CIMS_put )
reads_cims_putmod_indices = which( !is.na( match( c(1:Nreads), queryHits(reads_cims_putmod_pairs)) )  )
reads_cims_nonmod_indices = which(  is.na( match( c(1:Nreads), queryHits(reads_cims_putmod_pairs)) )  )

reads_cims_putmod         = reads_GR[ reads_cims_putmod_indices ]
reads_cims_nonmod         = reads_GR[ reads_cims_nonmod_indices ]


### @@@ this is more complicated with CIMS where there are repeats

# -------

reads_cits_putmod_pairs   = findOverlaps( reads_GR, CITS_put )
# find the indices of reads (from 1 to Nreads) that have at least one match in queryHits
reads_cits_putmod_indices = which( !is.na( match( c(1:Nreads), queryHits(reads_cits_putmod_pairs)) )  )

# Likewise, find the indices that have zero matches
reads_cits_nonmod_indices = which(  is.na( match( c(1:Nreads), queryHits(reads_cits_putmod_pairs)) )  )

reads_cits_putmod         = reads_GR[ reads_cits_putmod_indices ]
reads_cits_nonmod         = reads_GR[ reads_cits_nonmod_indices ]


modPairs = findOverlaps(reads_cits_putmod, CITS_put )
if (  !identical( queryHits( modPairs), c(1:length(modPairs)) )
      ){ stop("non-consecutive modcorrespondencePairs") }

# Create a meta-data column tracking the position of the modification within the kmer 
reads_cits_putmod$modposition   = start( CITS_put[ subjectHits(modPairs) ] ) - start( reads_cits_putmod[queryHits(modPairs)] ) +1

# Now "align" the strands, in the sense that we set reference_kmer to the actual kmer
# that went through the pore (i.e. take the reverse-compliment for all reads on the "-" strand)
reads_cits_putmod_strandAligned = strand_align(reads_cits_putmod)


# result = eval ( 
#  parse( 
#    text=paste0(
#      "refgen$",seqnames(ROI_GR),"[",as.character(lo_end),":",as.character(hi_end),"]" 
#    ) 
#  ) 
# )

#  ================================================
# write something specific for the m6A motif

m6A_motif  = "GGACT"
pore_model = read.csv("scripts/ref/pore_model_table.csv", row.names = 1)
pore_model_list <- setNames(split(pore_model, seq(nrow(pore_model))), rownames(pore_model))

test  <- seq_spec_compare ( seq        = m6A_motif,
                                         criterion  = "CITS",
                                         SOI_GR     = reads_cits_putmod, 
                                         control_GR = reads_cits_nonmod 
)

CITS_m6amotif_part <- seq_spec_compare ( seq        = m6A_motif,
                                         criterion  = "CITS",
                                         SOI_GR     = reads_cits_putmod_strandAligned, 
                                         control_GR = reads_cits_nonmod 
)




pos_seq_diffs = get_moddiff_vs_modpos_seq ( reads_cits_putmod_strandAligned, reads_cits_nonmod   )

#=================================


plot_seq_spec_comparison  ( CITS_m6amotif_part,
                            pore_model = pore_model,
                            mincurrent = 100,
                            maxcurrent = 140,
                            res = 0.5,
                            scale=TRUE )    

plot_seq_spec_comparison  ( test,
                            pore_model = pore_model,
                            mincurrent = 100,
                            maxcurrent = 140,
                            res = 0.5,
                            scale=TRUE )    

# ===================================================================
# now generalize the above to arbitrary sequence:

# GR2= expand_range(GRin, k=4, hg19_ref)

as.character(seqnames(CIMS_put[1]))
# This is how to get the sequence from the reference of a single GRange
get_refgen_seqs (  refgen = hg19_ref,
                   ROI_GR = CIMS_put[2], 
                   lead       = 2,
                   trail      = 2,
                   RNAstrand  = FALSE
)
# so now do it in a batch.
CITs_sequences = sapply( 1:length(CITS_put), function(x) get_refgen_seqs( refgen = hg19_ref, 
                                                                ROI_GR = CITS_put[x], 
                                                                lead=2, 
                                                                trail=2, 
                                                                RNAstrand=FALSE)  
              )
CIMs_sequences = sapply( 1:length(CIMS_put), function(x) get_refgen_seqs( refgen = hg19_ref, 
                                                                          ROI_GR = CIMS_put[x], 
                                                                          lead=2, 
                                                                          trail=2, 
                                                                          RNAstrand=FALSE)  
)

#====

CIMS_put_df = as.data.frame(CIMS_put) 
short_df=CIMS_put_df[1:10,]

command = cbind("hg19_ref$", paste0(short_df$seqnames,"[" , as.character(short_df$start),"]" ))
# collect "command strings", using the values from the database

apply( command, 1, function(x) eval( paste0(command[x,1], command[x,2])) )
# execute those strings as commands, and do it for each of the putative positions

allAtest= apply( c(1: length(command)), FUN=function(x) as.character( eval(parse(text=command[x])) ) )
                 )
#===

m6A_motif_subset  = dat_win_finite[ which( dat_win_finite$reference_kmer == m6A_motif),  ]

start( reads_GR[1:N]  )


