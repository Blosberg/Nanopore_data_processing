# library(Biostrings)
# read.DNAStringSet
# library("Biostrings")
# hg19_canon_ref = readDNAStringSet("/scratch/AG_Akalin/refGenomes/hg19_canon/hg19_canon.fa")
# negstrand_reads <- reads_GR[ strand( reads_GR)=="-",]
# source("https://bioconductor.org/biocLite.R") #--- ignore warnings telling you to upgrade
# biocLite("GenomicRanges")

# ===== LOAD GRANGES OBJECTS PRODUCED BY SNAKEMAKE ======

DIR_GR="/home/bosberg/projects/nanopore/pipeline_out/20180417_1233_polyA_RNA/06_GRobjects/"
file_GR="ec9bd1f6f626d73709089ba0bfc63929581e9ede_GR.RData"
load( file.path( DIR_GR, file_GR ) )

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

#  ================================================

# start importing lines from "analysis.R" here.
m6A_motif="GGACT"
Nreads   = length(reads_GR)

reads_cims_putmod_pairs   = findOverlaps( reads_GR, CIMS_put )
reads_cims_putmod_indices = which( !is.na( match( c(1:Nreads), queryHits(reads_cims_putmod_pairs)) )  )
reads_cims_nonmod_indices = which(  is.na( match( c(1:Nreads), queryHits(reads_cims_putmod_pairs)) )  )
reads_cims_putmod         = reads_GR[ reads_cims_putmod_indices ]
reads_cims_putNonmod      = reads_GR[ reads_cims_nonmod_indices ]

reads_cits_putmod_pairs   = findOverlaps( reads_GR, CITS_put )
reads_cits_putmod_indices = which( !is.na( match( c(1:Nreads), queryHits(reads_cits_putmod_pairs)) )  )
reads_cits_nonmod_indices = which(  is.na( match( c(1:Nreads), queryHits(reads_cits_putmod_pairs)) )  )
reads_cits_putmod         = reads_GR[ reads_cits_putmod_indices ]
reads_cits_putNonmod      = reads_GR[ reads_cits_nonmod_indices ]

# ===================================================================


plot_compare_motif(seq = m6A_motif, 
                   Gro1= reads_cits_putmod, 
                   Gro2= reads_cits_putNonmod,
                   mincurrent  = 110,
                   maxcurrent  = 140)

m6A_motif_subset  = dat_win_finite[ which( dat_win_finite$reference_kmer == m6A_motif),  ]

start( reads_GR[1:N]  )


