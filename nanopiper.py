#!/home/bosberg/projects/nanopiper/dev/guix/.guix-profile/bin/python3
# -*- python -*-
# nanopiper Pipeline.
#
# Copyright © 2017 Bren Osberg <brendan.osberg@mdc-berlin.de>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

# ==========================================================
# --- Import necessary modules -----------------------------

import argparse
import os, sys, json, csv, yaml
from os import path
from snakemake.utils import update_config
import shutil, filecmp
from glob import glob
import subprocess
import re # regular expressions.

# --- TO BE GENERALIZED: ----------------------------------


PATH_EXEC="/home/bosberg/projects/nanopiper/"
# ^ execution path from which pipeline executables (+repo) are stored.

snakefilename = PATH_EXEC+"Snakefile"
GUIX_PROFILE  = PATH_EXEC+"dev/guix/"
VERSION       = "0.1.0"

sys.path.append(PATH_EXEC+"scripts/00_SM/")
from func_defs import * 

# shouldn't need this; check back.
# path_rules_chunks="/home/bosberg/projects/nanopiper/scripts/00_SM/rules_chunks.py"
# include: path_rules_chunks

# --- DOCUMENTATION/HELP: ---------------------------------

description = """\
Nanopiper is a data processing pipeline, developed by
B. Osberg at MDC in Berlin in 2017-2018 for nanopore read data. 
It produces current and coverage information and can be used to 
produce information on modification and perhaps cell type.
"""

epilog = """\

Copyright 2019, 2198 Bren Osbeg
License GPLv3+: GNU GPL version 3 or later 
<http://gnu.org/licenses/gpl.html>.

This is free software: you are free to change and redistribute it.
There is NO WARRANTY, to the extent permitted by law.
"""
# --- BUILD ARGUMENT LIST: --------------------------------

def formatter(prog):
    return argparse.RawTextHelpFormatter(prog, max_help_position=80)

# Handling for additional arguments:
parser = argparse.ArgumentParser( description=description,
                                  epilog=epilog,
                                  formatter_class=formatter )

parser.add_argument('-v', '--version', action='version',
                    version=VERSION)
                     
parser.add_argument('-config_defaults', nargs='?', default=PATH_EXEC+'dev/config_defaults.json',
                    help="""\
The config file of default values --to be overwritten as needed.
""")

parser.add_argument('-config_userin', nargs='?', default=PATH_EXEC+'./config.json',
                    help="""\
The config file supplied by the user --to overwrite defaults as needed.
""")

parser.add_argument('config_npSM', nargs='?', default=PATH_EXEC+'config_npSM.json',
                    help="""\
The config file produced by this scripts and supplied to SnakeMake.
""")

parser.add_argument('--target', dest='target', action='append',
                    help="""\
Stop when the named target is completed instead of running the whole
pipeline.  The default target is "final-report".  Pass "--target=help"
to describe all available targets.""")

parser.add_argument('-n', '--dry-run', dest='dry_run', action='store_true',
                    help="""\
Only show what work would be performed.  Do not actually run the
pipeline.""")

parser.add_argument('--clustersub', dest="clustersub", default=False,
                    help="""\
Should we submit this Snakejob to an SGE cluster?
""")

parser.add_argument('--graph', dest='graph',
                    help="""\
Output a graph in PDF format showing the relations between rules of
this pipeline.  You must specify a graph file name such as
"graph.pdf".""")

parser.add_argument('--force', dest='force', action='store_true',
                    help="""\
Force the execution of rules, even though the outputs are considered
fresh.""")

parser.add_argument('--reason', dest='reason', action='store_true',
                    help="""\
Print the reason why a rule is executed.""")

parser.add_argument('--unlock', dest='unlock', action='store_true',
                    help="""\
Recover after a snakemake crash.""")

parser.add_argument('--verbose', dest='verbose', action='store_true',
                    help="""\
Print supplementary info on job execution.""")

parser.add_argument('--printshellcmds', dest='printshellcmds', action='store_true',
                    help="""\
Print commands being executed by snakemake.""")

args = parser.parse_args()

# --- PREPARE config: --------------------------------
   
config = prep_configfile( args.config_defaults,
                          args.config_userin,
                          args.config_npSM )

# @@@ DEBUGGING: REMOVE THIS:
# quit()

# --- Validate that pipeline can be executed:

# check for write access to refgenome dir
DIR_REFGENOME             = config['ref']['Genome_DIR']

if ( not os.access(DIR_REFGENOME, os.W_OK) ):
   print("Write access to refgenome folder is denied. Checking if necessary indexing files already exist: ... ")

   if( not os.path.isfile( os.path.join( DIR_REFGENOME , config['ref']['Genome_version']+ ".mmi")) ):
      bail("minimap index files not found, and cannot be created. Aborting")

   else:
      print("Refgenome index files are present. Continuing... ")

# --- Create symbolic links to PATHIN --- should be performed by rule at this point:
# 
# if ( input_data_type == "raw_minION"):
#    # the snakemake rules defined in the following included script assume
#    # that the minION output has been "chunked" into sequential files, and
#    # will have to be reassembled.
#    include   : os.path.join( config["scripts"]["script_folder"], config["scripts"]["rules_chunks"] )
# 
#    for sample in config["samplelist"]:
# 
#        # get subdir for this sample:
#        samplePATH = os.path.join(config["PATHIN"], config["samplelist"][sample]["subdir"] )
# 
#        # linke to each chunk directly:
#        for linkindex in config["samplelist"][sample]["chunkdirlist"]:
#           linkname = sample + "_" + str(linkindex) + config['samplelist'][sample]["fastq_suffix"]
# 
# 
#           source   = getPathCase( samplePATH, "fastq", "pass",config["samplelist"][sample]["fastq_prefix"] + str(linkindex) +  config["samplelist"][sample]["fastq_suffix"], "raw_minION")
# 
#           makelink(  source, os.path.join( config["PATHOUT"], config["samplelist"][sample]["subdir"], SUBDIR_SYMLINKS, linkname) )
# 
# 
# elif( input_data_type == "fastq"):
#    # Do some other stuff
#    include   : os.path.join( config["scripts"]["script_folder"], config["scripts"]["rules_wholefastq"] )
# 
#    for sample in config["samplelist"]:
# 
#        # get subdir for this sample:
#        samplePATH = os.path.join(config["PATHIN"], config["samplelist"][sample]["subdir"] )
# 
#        linkname = sample + config['samplelist'][sample]["fastq_suffix"]
#        makelink( os.path.join(samplePATH, config["samplelist"][sample]["fastq_prefix"] + config["samplelist"][sample]["fastq_suffix"] ), os.path.join( config["PATHOUT"], config["samplelist"][sample]["subdir"], SUBDIR_SYMLINKS, linkname))
# 
# else:
#    print("Unrecognized input data format. Terminating.")
#    exit(1)
# 

# --- Define command to run Snakemake with ---------


command = [
    os.path.join(GUIX_PROFILE, ".guix-profile", "bin", "snakemake"),
    "--snakefile={}".format(snakefilename),
    "--configfile={}".format(args.config_npSM),
    "--jobs={}".format( str(config['execution']['jobs']) ),
]

if ( args.clustersub ):
# Jobs to be submitted to an SGE cluster:

    print("Commencing snakemake run submission to cluster", flush=True, file=sys.stderr)

    # cluster_config generation is handled in func_defs.
    # generate_cluster_configuration( config )

    # Check if a particular queue has been specified: 
    if ( ('queue' in config['execution']['cluster'] )  or ( 'queue' in config['execution']['rules']['__default__'] )):
        queue_selection_string = " -q {cluster.queue} "
	# User has defined queue name(s) --either generally, or rule-specific.
	# in either case, generate_cluster_configuration() (defined above) has 
        # already printed the apropriate queue name to each rule.
        queue_selection_string = " -q {cluster.queue} "
    else :
        # User has supplied no q value -> let SGE pick the the default.
        queue_selection_string = ""
    # --- done checking if ( queue name(s) supplied by user)

    # check if a contact email has been specified:
    if config['execution']['cluster']['contact-email'].lower() == 'none':
        contact_email_string = ""
    else:
        contact_email_string = "-m a -M %s" % config['execution']['cluster']['contact-email']


    # check if machine allows for queue submission at all
    try:
        subprocess.call(["qsub", "-help"], stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    except OSError as e:
        if e.errno == os.errno.ENOENT:
            print("Error: Your system does not seem to support cluster submission.\nReason: Could not find qsub executable.", file=sys.stderr)
            exit(1)
        else:
            raise

    # create path for Snakemake e/o logfiles: 
    ClusterLogsDir = os.path.join(config['PATHOUT'], 'cluster_log_files')
    os.makedirs( ClusterLogsDir , exist_ok=True)


    # define qsub command:
    # The following were removed (and might be necessary): 
    # " -v R_LIBS_USER "
    # "--jobscript={}/qsub-template.sh".format(config['locations']['pkglibexecdir']),
    qsub = "qsub -e " + ClusterLogsDir + " -o " + ClusterLogsDir + " -v PATH -v GUIX_LOCPATH  %s -l h_stack={cluster.h_stack} -l h_vmem={cluster.MEM} %s -b y -pe smp {cluster.nthreads} -cwd" % ( queue_selection_string, contact_email_string)
    if config['execution']['cluster']['args']:
        qsub += " " + config['execution']['cluster']['args']
    command += [
        "--cluster-config={}".format(cluster_config_file),
        "--cluster={}".format(qsub),
        "--latency-wait={}".format(config['execution']['cluster']['missing-file-timeout'])
    ]
else:
# Jobs to be submitted on local machine:
    print("Commencing snakemake run submission locally", flush=True, file=sys.stderr)

command.append("--rerun-incomplete")
if args.graph:
    # Only output dag graph:
    command.append("--dag | dot -Tpdf > " + config["PATHOUT"] + "dag.pdf")

else:
    # Check for additional arguments/flags:
    prepare_links()
    display_logo()
    if args.force:
        command.append("--forceall")
    if args.dry_run:
        command.append("--dryrun")
    if args.reason:
        command.append("--reason")
    if args.unlock:
        command.append("--unlock")
    if args.verbose:
        command.append("--verbose")
    if args.printshellcmds:
        command.append("--printshellcmds")
    if args.target and 'help' in args.target:
        command.append("help")

    # Now actually execute the pipeline:
    subprocess.run(command)
