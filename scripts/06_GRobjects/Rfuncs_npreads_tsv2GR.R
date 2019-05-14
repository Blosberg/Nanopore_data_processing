# ========================================================================
# --- Extract strand information based on relationship between event_indices and position.

assign_strand <-  function( Datin = stop("Datin must be provided"),
                            perform_sanity_checks = FALSE )
{

  Datin$strand = NULL

  dtest = Datin[c("read_index", "position", "event_index") ] %>%
    #  subset(read_index %in% 1:10) %>%              # minimal check that things are working.
    group_by(read_index) %>%                         # First partition by read_index
    # mutate(position_index = order(position)) %>%     # Set position_index variable defined as the order of
    # "position" (lowest to highest, with equal values taken
    # in order of row number -- all within each read_index)

    mutate(  diff_event    = c(0, diff( event_index ) ) ) %>%
    mutate(  diff_position = c(0, diff( position    ) ) ) %>%   # adds diff_index col which is order of event - order of position

    mutate(strand = case_when(
      all( diff_event * diff_position >= 0 ) ~ '+',
      all( diff_event * diff_position <= 0 ) ~ '-',
      TRUE ~ '*'
    ) )

  dtest$diff_event    <- NULL
  dtest$diff_position <- NULL

  ftest <- left_join(Datin, dtest, by=c('read_index','position','event_index') )

  # Check that they are all attributed correctly

  if( perform_sanity_checks )
    {
    Watson = ftest[ ftest$strand =="+", ];
    Crick  = ftest[ ftest$strand =="-", ];
    Unk    = ftest[ ftest$strand =="*", ];

    if( dim(Unk)[1] > 0 )
      {
      writeLines("ERROR: reads being assigned unknown strand", logFile );
      stop(paste("ERROR: reads being assigned unknown strand"))
      }
    }

  return(ftest)
}

# ===============================================================
# make sure the specific position parameters make sense.
sanity_check_tbl_spec_position  <- function( readposn_tbl_in        = stop("read_GR_in must be provided")
                                       )
{
    posn            = unique( readposn_tbl_in$position )
    reference_kmer  = unique( readposn_tbl_in$reference_kmer )
    model_kmer      = unique( readposn_tbl_in$model_kmer     )

      if ( max ( length(posn),
                 length(reference_kmer),
                 length(model_kmer)   ) > 1
          )
    {
      stop("non-unique values for position-constant parameters")
      }
    else{
      return( 0 )
      }
}
# ===============================================================
# Take a read in table format and eliminate redundancy (average overlapping positions)

flatten_read_tbl  <- function( read_tbl_in         = stop("read_GR_in must be provided"),
                               perform_sanity_checks = FALSE
                             )
{
  # --- takes a specific read and combines all events associated with the same position
  # --- producing a single value for each position.
  # --- skipped positions are left as "NA"s

  if ( perform_sanity_checks )
    { #--- if we want to be safe.
    chr         = unique( read_tbl_in$contig )
    read_index  = unique( read_tbl_in$read_index )
    if( length( chr) != 1 || length( read_index ) != 1 )
      { stop("irregular seqnames or read_indices in read") }
  } else { #--- if we are confident, and want to be fast.
    chr         =  read_tbl_in$contig[1]
    read_index  =  read_tbl_in$read_index[1]
    }

  read_tbl_in_Ns_removed <- read_tbl_in[ read_tbl_in$model_kmer != "NNNNN" ,  ]
  # ----------

  # create a list, with each element containing groups of GR objects with the same start position
  read_tbl_in_splitby_pos <- split(   read_tbl_in_Ns_removed,
                                      read_tbl_in_Ns_removed$position )

    if( perform_sanity_checks )
    {
    # If all is ok, it should return 0 everywhere.
    sanity  <- unlist( lapply( read_tbl_in_splitby_pos,
                       function(x) sanity_check_tbl_spec_position ( readposn_tbl_in = x)
                              )
                       ) # length( which( sanity != 0 ) )
    }


    # for each starting position, compact the events into a single observation.
    #, so run lapply over the list of positions.
    return( data.table::rbindlist( lapply( read_tbl_in_splitby_pos,
                                           function(x) flatten_tbl_at_spec_position( readposn_tbl_in = x,
                                                                                     chr_in          = chr
                                                                                    )
                                          )
                                )
            )
    # ... and then return the whole r-binded list as the flattened read.

}

# ===============================================================
# Define, precisely, the one-element GRanges object for each position in the above function.
flatten_tbl_at_spec_position  <- function( readposn_tbl_in        = stop("read_GR_in must be provided"),
                                           chr_in                = stop("chr_in must be provided"),
                                           perform_sanity_checks = FALSE
                                          )
{
if( dim( readposn_tbl_in )[1] == 1 )
  { return( readposn_tbl_in ) }


    event_index_str  = paste0( as.character( min ( readposn_tbl_in$event_index ) ),
                               "-",
                               as.character(max ( readposn_tbl_in$event_index ) ) )

#    model_kmer       = readposn_tbl_in$model_kmer[1]
#    posn             = readposn_tbl_in$position[1]
#    read_index       = readposn_tbl_in$read_index[1]
   # ^ These have already been checked to be unique assuming sanity_checks were activated.
#    reference_kmer  = readposn_tbl_in[1]$reference_kmer
#    model_mean      = readposn_tbl_in[1]$model_mean
#    model_stdv      = readposn_tbl_in[1]$model_stdv
   # ^ redundant metadata columns omitted for efficiency.

  # duration of overlapping events:
  passage_time = sum(  readposn_tbl_in$event_length)

  event_mean   = sum(  readposn_tbl_in$event_level_mean * readposn_tbl_in$event_length ) /passage_time
  event_stdv   = sqrt( sum( readposn_tbl_in$event_stdv*readposn_tbl_in$event_stdv * readposn_tbl_in$event_length ) / passage_time )

  output  <- readposn_tbl_in[1,]

  output$event_index      <- event_index_str
  output$event_level_mean <- event_mean
  output$event_stdv       <- event_stdv
  output$event_length     <- passage_time

  return(output)
}

# ===============================================================
# Inefficiently double-check that the specific position parameters make sense.
sanity_check_tbl_spec_position  <- function( readposn_tbl_in        = stop("read_GR_in must be provided")
                                       )
{
    S               = unique( readposn_tbl_in$position )
    reference_kmer  = unique( readposn_tbl_in$reference_kmer )
    model_kmer      = unique( readposn_tbl_in$model_kmer     )

      if ( max ( length(S),
                 length(reference_kmer),
                 length(model_kmer)   ) > 1
          )
    {
      stop("non-unique values for position-constant parameters")
      }
    else{
      return( 0 )
      }
}

# ===============================================================
# Save the observed model_parameters:
writeout_pore_model  <- function( event_dat   = stop("read_GR_in must be provided"),
                                  fout        = stop("output must be provided"),
                                  strand_type = "RNA"
                                )
{
model_dat = split( event_dat,
                   event_dat$model_kmer )

model_mean_list <- lapply( model_dat, function(x) unique(x$model_mean) )
model_stdv_list <- lapply( model_dat, function(x) unique(x$model_stdv) )

test_length     <- unlist( lapply( 1:length(model_mean_list), function(x) length(model_mean_list[[x]]) + length(model_stdv_list[[x]]) ) )

if( min(test_length != 2 ) || max (test_length !=2 ))
  { stop("ERROR: irregular model values found in table data. Terminating.") }

model_dat <- cbind( unlist(model_mean_list), unlist(model_stdv_list))
colnames(model_dat) <- c("mean", "stdv")

model_dat <- model_dat[ rownames(model_dat) != "NNNNN", ]

if( strand_type == "RNA" )
  { rownames(model_dat) <- gsub("T", "U", rownames(model_dat) ) }

write.table( model_dat,
	     file = fout,
	     sep="\t")
}

# ===============================================================
# Extract Event information from the nanopolish output:
get_event_dat  <- function( Event_file_list = stop("Datin must be provided"),
                            logFile         = stop("logFile must be provided"))
  {
    colnames <-
      c(
        "contig",
        "position",
        "reference_kmer",
        "read_index",
        "strand",
        "event_index",
        "event_level_mean",
        "event_stdv",
        "event_length",
        "model_kmer",
        "model_mean",
        "model_stdv",
        "standardized_level"
      )
    dat_all = data.frame()
    readcount_offset = 0

    for (i in c(1:length(Event_file_list)))  {
      # Event_file      = paste0(    "Ealign_", as.character(i),".tsv")
      fin               = as.character(file.path(Event_file_list[i]))

      write(paste("---Reading tsv file from: ", fin),
            file   = logFile,
            append = TRUE)

      dat_temp        = read.csv(
        file = fin,
        sep  = '\t',
        stringsAsFactors = FALSE,
        header = TRUE
      ) #--- read in the current "chunk" of np data

      dat_temp$read_index <- dat_temp$read_index + readcount_offset  #--- offset the read_index values to ensure uniqueness

      dat_all           =   rbind( dat_all, dat_temp )               #--- compile the chunks together
      readcount_offset  =  (max(dat_all$read_index)  + 1)            #--- re-calculate the necessary offset

    }

    Nreads = length(unique(dat_all$read_index))

    # cast the read_indices as factors (to make unique values sequentially increasing), and then
    # cast them back into integers (since some read_indices have been removed for quality/alignment/etc. reasons
    dat_all$read_index <- as.numeric(as.factor(dat_all$read_index))


    return(dat_all)
  }

# ===============================================================
# Convert list of reads in table format to GRL:
convert_tbl_readlist_to_GRL  <- function( Readlist_in   = stop("Readlist must be provided"),
                                          Flatten_reads = TRUE )
{
  if(Flatten_reads)
  {
    # tic <- Sys.time()
    Flatreadlist <-
      lapply( Readlist_in, function(x)
              flatten_read_tbl ( read_tbl_in   = x,
                                 perform_sanity_checks = FALSE))
    # toc <- Sys.time()

    result <- GRangesList(lapply(Flatreadlist, function(x)
      GRanges(
        seqnames     = x$contig,
        strand       = x$strand,
        range        = IRanges( start = x$position,
                                end   = x$position),
        read_index   = x$read_index,
        event_index  = x$event_index,
        event_mean   = x$event_level_mean,
        event_stdv   = x$event_stdv,
        event_length = x$event_length,
        model_kmer   = x$model_kmer
      )))
  } else {

    result <- GRangesList( lapply( Readlist_in, function(x)
      GRanges(
        seqnames     = x$contig,
        strand       = x$strand,
        range        = IRanges( start = x$position,
                                end   = x$position),
        read_index   = x$read_index,
        event_index  = x$event_index,
        event_mean   = x$event_level_mean,
        event_stdv   = x$event_stdv,
        event_length = x$event_length,
        model_kmer   = x$model_kmer
      )))

  }

  return ( result )
}
