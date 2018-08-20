# Environment from this execution is stored in : "/home/bosberg/projects/nanopore/signal_processing_Rworkspace.RData"
# import nanopore data and run some statistics on them:
# rm(list=ls()) # CLEAN UP EVERYTHING

## Collect arguments
args <- commandArgs(TRUE)

## Default setting when no arguments passed
if(length(args) < 1) {
  args <- c("--help")
}

## Help section
if("--help" %in% args) {
  cat("
      Render to report
      
      Arguments:
      Rfuncs_tableGRconv_file --script with function definitions used here.
      output      --filename for the output GRanges .RData 
      Ealign_files -- array of files to use as input
      logFile      -- filename to pipe output to 
     
      Example:
      ./test.R --arg1=1 --arg2='output.txt' --arg3=TRUE \n\n")
  
  q(save="no")
}

## Parse arguments (we expect the form --arg=value)
parseArgs <- function(x) strsplit(sub("^--", "", x), "=")
argsDF    <- as.data.frame(do.call("rbind", parseArgs(args)))
argsL     <- as.list(as.character(argsDF$V2))

names(argsL) <- argsDF$V1

# catch output and messages into log file
# out <- file(argsL$logFile, open = "wt")
# sink(out, type = "output")
# sink(out, type = "message")

# Run Functions -----------------------------------------------------------

# e.g. (replace this list with actual arguments)
Rfuncs_tableGRconv_file  <- argsL$Rfuncs_tableGRconv_file 
output                   <- argsL$output
logfile                  <- argsL$logFile
Ealign_files             <- unlist( strsplit(argsL$Ealign_files,",")  )

#===============================================================
suppressPackageStartupMessages( library(GenomicRanges) )
suppressPackageStartupMessages( library(dplyr)         )
source(Rfuncs_tableGRconv_file)

dat_all         = get_event_dat( Ealign_files )

read_list_final = unique( dat_all$read_index  )
Nreads          = length( read_list_final     )

if( !( identical( as.numeric(c(1:Nreads)) , as.numeric(read_list_final) )  ))
  {  writeLines("ERROR: final read list is not step-wise increasing", logfile );
     exit(1)

  }


# ========================================================================================
# dat_all$event_level_mean[1:20]
# hist(dat_all$event_level_mean, breaks=500)
mincurrent=50
maxcurrent=150

current_window = c( mincurrent, maxcurrent )

temp           = dat_all[  dat_all$event_level_mean  > mincurrent, ]
dat_windowed   = temp[ temp$event_level_mean < maxcurrent, ] 
dat_win_finite = dat_windowed [ which (! is.na (dat_windowed$event_level_mean) ), ]
rm(temp)


# =======================
# Add strand information:
dat_win_finite_stranded = assign_strand( dat_win_finite )

# Check that they are all attributed correctly
Watson = dat_win_finite_stranded[ dat_win_finite_stranded$strand =="+", ];
Crick  = dat_win_finite_stranded[ dat_win_finite_stranded$strand =="-", ];
Unk    = dat_win_finite_stranded[ dat_win_finite_stranded$strand =="*", ];

if( dim(Unk)[1] > 0 )
  { writeLines("ERROR: reads being assigned unknown strand", logfile );
    stop(paste("ERROR: reads being assigned unknown strand")) }

# ================================================
# BUILD GRanges OBJECT TO Process

reads_GR = GRanges( seqnames   = dat_win_finite_stranded$contig, 
                    strand     = dat_win_finite_stranded$strand, 
                    IRanges(   start  =  dat_win_finite_stranded$position+1,
                               end    =  dat_win_finite_stranded$position+5 ), # --- nopolish provides position before beginning of 5bp window.
                    
                    reference_kmer = dat_win_finite_stranded$reference_kmer,
                    read_index     = dat_win_finite_stranded$read_index,
                    event_index    = dat_win_finite_stranded$event_index,
                    event_mean     = dat_win_finite_stranded$event_level_mean,
                    event_stdv     = dat_win_finite_stranded$event_stdv,
                    event_length   = dat_win_finite_stranded$event_length,
                    model_mean     = dat_win_finite_stranded$model_mean,
                    model_stdv     = dat_win_finite_stranded$model_stdv
)



reads_GRL = split( reads_GR, 
                   reads_GR$read_index )

# rm(dat_all)
# rm( dat_win_finite)
# rm( dat_win_finite_stranded)
# rm( dat_stranded)
# rm( Watson)
# rm( Crick )
# rm( Unk )
# rm( dat_windowed )
# save.image(file="/home/bosberg/projects/nanopore/signal_processing_Rworkspace.RData")

save( file= output, reads_GR, reads_GRL )
