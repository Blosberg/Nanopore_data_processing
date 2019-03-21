#!/home/bosberg/projects/nanopore/dev/guix/.guix-profile/bin/python3
import os, sys, json, csv, yaml
import argparse

# import IPython;

#-----------------------------------------------------------------
#--- Define settings for config file -----------------------------

# import default settings:
config_defaults="dev/config_defaults.json"
config = yaml.safe_load(open(config_defaults, 'r'))

include: os.path.join( config["scripts"]["script_folder"], config["scripts"]["pyfunc_defs"] )

# import user-supplied settings, over-writing as necessary
config_userin="config.json"
update_config( config, yaml.safe_load( open( config_userin, 'r')))

# look in the path provided for each sample and build a list of "chunks" in
# each case.  these are just directory names and should be simply numbers
# (0,1,2,...)
for sLoop_countchunks in config["samplelist"]:
    # set path within this sample's data set to look for chunks
    DATPATH = os.path.join(config["PATHIN"], config["samplelist"][sLoop_countchunks]["subdir"], "fast5", "pass" )

    # define a list of these subentries:
    unsorted_stringlist =  [ entry.name for entry in os.scandir( DATPATH ) if entry.is_dir() ]
    intlist = list( map(int, unsorted_stringlist) )
    intlist.sort()   

    config["samplelist"][sLoop_countchunks]["chunkdirlist"] = list( map(str, intlist)) 


# store a log of the config file we just built:
config_log="configlog_out.json"

with open(config_log, 'w') as outfile:
    dumps = json.dumps(config,
                       indent=4, 
                       sort_keys=True,
                       separators=(",",": "), 
                       ensure_ascii=True)
    outfile.write(dumps)


# ------------------------------------------------------
#--- Define Dependencies:


MM2        = config["progs"]["minimap"]
SAMTOOLS   = config["progs"]["SAMTOOLS"]
nanopolish = config["progs"]["nanopolish"]

RefTranscriptome   = config["ref"]["Transcriptome"]
GENOME_VERSION     = config["ref"]["Genome_version"]
RmdReportScript    = os.path.join(config["scripts"]["script_folder"],"final_report","Nanopore_report.Rmd")
input_data_type    = config["input_data_type"]

R_tables2GR_main     = os.path.join( config["scripts"]["script_folder"], config["scripts"]["Rmain_table2GRconv"] )
R_tables2GR_funcs    = os.path.join( config[ "scripts"]["script_folder"], config[ "scripts"]["Rfuncs_table2GRconv"] )

R_flattenreads_main  = os.path.join( config["scripts"]["script_folder"], config["scripts"]["Rmain_flattenreads"] )
R_flattenreads_funcs = os.path.join( config["scripts"]["script_folder"], config["scripts"]["Rfuncs_flattenreads"] )

R_build_histlist_main   = os.path.join( config["scripts"]["script_folder"], config["scripts"]["Rmain_build_histlist"] )

#------------------------------------------------------
#--- Define output directories

SUBDIR_SYMLINKS           = "01_symlinks_fastqindex_chunks/"
# -- all the rest will be created by snakemake automatically

SUBDIR_ALIGNED_MINIMAP    = "02_MM_aligned_chunks/"
SUBDIR_FILTERED_MINIMAP   = "03_MM_filtered_chunks/"
SUBDIR_SORTED_MINIMAPPED  = "04_MM_sortedbam/"
SUBDIR_EVENTALIGN         = "05_eventalign/"
SUBDIR_GR                 = "06_GRobjects"
SUBDIR_REPORT             = "Final_report/"

DIR_REFGENOME             = config['ref']['Genome_DIR']

#------------------------------------------------------
#--- Check that the pipeline can be executed:

# check for write access to refgenome dir
if ( not os.access(DIR_REFGENOME, os.W_OK) ):
   print("Write access to refgenome folder is denied. Checking if necessary indexing files already exist: ... ")

   if( not os.path.isfile(DIR_REFGENOME , config['ref']['Genome_version'] + ".mmi")):
      bail("minimap index files not found, and cannot be created. Aborting")

   else:
      print("Refgenome index files are present. Continuing... ")

#--- Create symbolic links to PATHIN so that indexing/etc can be performed in
# written pathout

if ( input_data_type == "raw_minION"):
   # the snakemake rules defined in the following included script assume
   # that the minION output has been "chunked" into sequential files, and
   # will have to be reassembled.
   include   : os.path.join( config["scripts"]["script_folder"], config["scripts"]["rules_chunks"] )

   for sLooplinki in config["samplelist"]:

       os.makedirs( os.path.join( config["PATHOUT"], config["samplelist"][sLooplinki]["subdir"], SUBDIR_SYMLINKS),  exist_ok=True)

       # get subdir for this sample:
       samplePATH = os.path.join( config["PATHIN"], config["samplelist"][sLooplinki]["subdir"] )

       # linke to each chunk directly:
       for linkindex in config["samplelist"][sLooplinki]["chunkdirlist"]:
          linkname = sLooplinki + "_" + str(linkindex) + "." + config["fastq_suffix"]

          target   = getPathCase( samplePATH, "fastq", "pass",config["samplelist"][sLooplinki]["fastq_prefix"] + str(linkindex) + "." +  config["fastq_suffix"], "raw_minION")

          linkloc  = os.path.join( config["PATHOUT"], config["samplelist"][sLooplinki]["subdir"], SUBDIR_SYMLINKS, linkname)

          makelink(  target, linkloc )


elif( input_data_type == "fastq"):
   # Do some other stuff
   include   : os.path.join( config["scripts"]["script_folder"], config["scripts"]["rules_wholefastq"] )

   for sLoop2linki in config["samplelist"]:

       # get subdir for this sample:
       samplePATH = os.path.join(config["PATHIN"], config["samplelist"][sLoop2linki]["subdir"] )

       linkname = sLoop2linki + config['samplelist'][sLoop2linki]["fastq_suffix"]
       makelink( os.path.join(samplePATH, config["samplelist"][sLoop2linki]["fastq_prefix"] + config["samplelist"][sLoop2linki]["fastq_suffix"] ), os.path.join( config["PATHOUT"], config["samplelist"][sLoop2linki]["subdir"], SUBDIR_SYMLINKS, linkname))

  
else:
   print("Unrecognized input data format. Terminating.")
   exit(1)


#------------------------------------------------------
#--- Define output (target) files:

# Initialize empty list
OUTPUT_FILES = []

# increment by sample:
for sampleLoopi_targets in config["samplelist"]:
  
   # @@@ TODO: implement sample-dependent targets with defaults.
   if ( config["target_out"] == "report" ):
      OUTPUT_FILES.extend(  os.path.join( config["PATHOUT"], config["samplelist"][sampleLoopi_targets]["subdir"], SUBDIR_REPORT, "" + sampleLoopi_targets + "_report.html")  )
   elif ( config["target_out"] == "histlist" ):
      OUTPUT_FILES.extend( 
                     os.path.join( config["PATHOUT"], config["samplelist"][sampleLoopi_targets]["subdir"], SUBDIR_GR, sampleLoopi_targets + "_kmer_histlist.rds")  
                          )
   elif ( config["target_out"] == "flatreads_GRL" ):
      OUTPUT_FILES.extend( 
                          [ os.path.join( config["PATHOUT"], config["samplelist"][sampleLoopi_targets]["subdir"], SUBDIR_GR,  sampleLoopi_targets + "_reads_flat_GRL.rds") ] 
                         )
   elif ( config["target_out"] == "reads_GRL" ):
      OUTPUT_FILES.extend( 
                          [ os.path.join( config["PATHOUT"], config["samplelist"][sampleLoopi_targets]["subdir"], SUBDIR_GR, sampleLoopi_targets + "_reads_GRL.rds") ]
                          )
   elif ( config["target_out"] == "mergedbam" ):
      OUTPUT_FILES.extend( 
                          [ os.path.join( config["PATHOUT"], config["samplelist"][sampleLoopi_targets]["subdir"], SUBDIR_SORTED_MINIMAPPED, sampleLoopi_targets + ".sorted.bam") ]
                         )
   elif ( config["target_out"] == "aligned_chunks"):
      OUTPUT_FILES.extend(  
                          get_chunkfiles( sampleLoopi_targets, os.path.join( config["PATHOUT"], config["samplelist"][sampleLoopi_targets]["subdir"], SUBDIR_ALIGNED_MINIMAP) ,  "", ".sam", 0 )
                         )
   else:
      print("Unrecognized target output file format: ", config["target_out"], " ... Terminating.")
      exit(1)

  
#---  DEBUGGING:
#------------------------------------------------------
# print("input_data_type = " + config["input_data_type"])
# print("target out = " + config["target_out"])
# IPython.embed()
print("---- last check before rules: ------ ")
print ( "len(OUTPUT_FILES)=")
print ( len( OUTPUT_FILES) )
print(" OUTPUT_FILES=")
for x in OUTPUT_FILES:
  print(x)
print("\n finished outputting output files \n\n ")

#
# ========================================================================
#
#   BEGIN RULES
#
# ========================================================================

rule all:
    input:
        OUTPUT_FILES

#------------------------------------------------------
rule make_report:
    # build the final output report in html format
    input:
        aligned_reads_bam = os.path.join( config["PATHOUT"], "{wcreport_sampleDir}", SUBDIR_SORTED_MINIMAPPED, "run_{wcreport_samplename}.sorted.bam"),
        transcriptome     = RefTranscriptome,
        GRobj             = os.path.join( config["PATHOUT"], "{wcreport_sampleDir}", SUBDIR_GR, "{wcreport_samplename}_reads_GRL.rds")
    output:
        os.path.join( config["PATHOUT"], "{wcreport_sampleDir}", SUBDIR_REPORT, "{wcreport_samplename}_report.html")
    params:
        " readcov_THRESH = 10;   ",
        " yplotmax = 10000; "
    log:
        logfile = os.path.join( config["PATHOUT"], "{wcreport_sampleDir}", SUBDIR_REPORT, "final_report_{wcreport_samplename}.log")
    message:
        fmt("Producing final report")
    shell:
        " Rscript -e  '{params} "
        " fin_readalignment_bam = \"{input.aligned_reads_bam}\"; "
        " fin_Transcript        = \"{input.transcriptome}\";"
        " fin_GRobj             = \"{input.GRobj}\";"
        " Genome_version        =\"{GENOME_VERSION}\"; "
        " rmarkdown::render(\"{RmdReportScript}\", output_file = \"{output}\" ) ' "

#------------------------------------------------------
rule np_event_align:
    # Align the events to the reference genome.
    # The wildcard "chunk" can simply be "full", in cases
    # where there are no chunks
    input:
        sortedbam             = os.path.join( config["PATHOUT"], "{wcEvalign_sampleDir}", SUBDIR_SORTED_MINIMAPPED, "bam_chunks", "{wcEvalign_samplename}_{wcEvalign_chunk}.sorted.bam"),
        NOTCALLED_indexedbam  = os.path.join( config["PATHOUT"], "{wcEvalign_sampleDir}", SUBDIR_SORTED_MINIMAPPED, "bam_chunks", "run_{wcEvalign_samplename}_{wcEvalign_chunk}.sorted.bam.bai"),
        fastq_file            = os.path.join( config["PATHOUT"], "{wcEvalign_sampleDir}", SUBDIR_SYMLINKS,  "{wcEvalign_samplename}_{wcEvalign_chunk}" +config["fastq_suffix"]),
        fastq_npi             = os.path.join( config["PATHOUT"], "{wcEvalign_sampleDir}", SUBDIR_SYMLINKS,  "{wcEvalign_samplename}_{wcEvalign_chunk}" + config["fastq_suffix"] + ".index"),
        refgenome_fasta       = os.path.join( DIR_REFGENOME, config['ref']['Genome_version']+ ".fa" ),
        NOTCALLED_bwt         = os.path.join( DIR_REFGENOME, config['ref']['Genome_version']+ ".fa.bwt"),
        NOTCALLED_pac         = os.path.join( DIR_REFGENOME, config['ref']['Genome_version']+ ".fa.pac")
    output:
        Evaligned         = os.path.join( config["PATHOUT"], "{wcEvalign_sampleDir}", SUBDIR_EVENTALIGN, "csv_chunks", 'Ealign_{wcEvalign_samplename}_{wcEvalign_chunk}.csv' )
    log:
        logfile  = os.path.join( config["PATHOUT"], "{wcEvalign_sampleDir}", SUBDIR_EVENTALIGN, "csv_chunks", 'Ealign_{wcEvalign_samplename}_{wcEvalign_chunk}.log')
    message: """---- Align events from sample {wildcards.wcEvalign_samplename}, chunk {wildcards.wcEvalign_chunk} to the genome ----"""
    shell:
        " {nanopolish} eventalign --reads {input.fastq_file} --bam {input.sortedbam} --genome {input.refgenome_fasta} --scale-events  > {output}  2> {log.logfile} "

#------------------------------------------------------
# rule quickcheck: (TODO)
#------------------------------------------------------

rule index_sortedbam:
    # Index the sorted bam file with samtools
    input:
        sortedbam  = os.path.join( config["PATHOUT"], "{wcindexbam_sampleDir}", SUBDIR_SORTED_MINIMAPPED, "bam_chunks", "run_{wcindexbam_samplename}_{wcindexbam_chunk}.sorted.bam")
    output:
        indexedbam = os.path.join( config["PATHOUT"], "{wcindexbam_sampleDir}", SUBDIR_SORTED_MINIMAPPED, "bam_chunks", "run_{wcindexbam_samplename}_{wcindexbam_chunk}.sorted.bam.bai")
    log:
        logfile    = os.path.join( config["PATHOUT"], "{wcindexbam_sampleDir}", SUBDIR_SORTED_MINIMAPPED, "bam_chunks", 'run_{wcindexbam_samplename}_{wcindexbam_chunk}_samtoolsindex.log')
    message: """---- index the bam files for {wildcards.sample} chunk {wildcards.chunk} ----"""
    shell:
        " {SAMTOOLS} index  {input.sortedbam}  2> {log.logfile} "

#------------------------------------------------------

