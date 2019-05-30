import os, sys, json, csv, yaml
import re # regular expressions.
from snakemake.utils import update_config

# --------------------------------------------------------------
class Rule_wc_struct(object):
    def __init__(self, SampleName_in, SampleDirList_in):
        self.wc_sampleName = SampleName_in
        self.wc_sampleDirs = SampleDirList_in

# --------------------------------------------------------------

def bail( msg ):
    # Terminate, with helpful error message:
    """Print the error message to stderr and exit."""
    print("ERROR: " + msg + "... exiting.", file=sys.stderr)
    exit(1)

# --------------------------------------------------------------

def dumpjson( yamldatin, fout ):
    with open(fout, 'w') as outfile:
        dumps = json.dumps(yamldatin,
                           indent=4, sort_keys=True,
                           separators=(",",": "), ensure_ascii=True)
        outfile.write(dumps)

# --------------------------------------------------------------

def nice(cmd, args, log=None):
    # executable = tool(cmd)
    # line = ["nice -" + str(config['execution']['nice']), executable] + [toolArgs(cmd)] + args
    line = ["nice -" + str(config['execution']['nice']), cmd] + args
    if log:
        line.append("> {} 2>&1".format(log))
    return " ".join(line)

# --------------------------------------------------------------

def fmt(message):
    """Format the MESSAGE string."""
    return "----------  " + message + "  ----------"

# --------------------------------------------------------------

def getPathCase( mainpath, subd1, subd2, filename, input_data_type ):
    # Returns a path dependent on the input data type:
    if ( input_data_type == "raw_minION" ):
       result = os.path.join( mainpath, subd1, subd2, filename )
    elif( input_data_type == "fastq" ):
        result = os.path.join( mainpath, filename ),
    else:
        print("ERROR: Unrecognized input_data_type: " + input_data_type + " in getPathCase")
        exit(1)
    return( result )

# --------------------------------------------------------------
def get_chunkfiles( wcstruct, DIR, subDir, prefix_string, suffix_string, quoted):

    FILES_list = []

    for sampleDirName_i in range( len( config["samplelist"][wcstruct.wc_sampleName]["sampleDirNames"])):
        sampleDirName = config["samplelist"][wcstruct.wc_sampleName]["sampleDirNames"][sampleDirName_i]

        Sample_indices_str     = config["samplelist"][ wcstruct.wc_sampleName ]["sampleDirs"][sampleDirName]["chunkdirlist"]

        FILES_list.extend( [ os.path.join( DIR, sampleDirName, subDir, prefix_string + wcstruct.wc_sampleName + "_" + chunk + suffix_string ) for chunk in Sample_indices_str ] )
    # now we've collected the list of files across each subdirectory and all the corresponding "chunks" in each of them.

    if ( quoted ):
       return( ",".join( FILES_list )  )
    else:
       return( FILES_list )

# --------------------------------------------------------------

def prep_configfile( args ):
    # Create the SM config file, from default, and user-defined inputs:

    config_defaults = args.config_defaults
    config_userin   = args.config_userin
    config_npSM     = args.config_npSM

    if ( not os.path.isfile(config_defaults) ):
        bail("Cannot find default configfile: " + config_defaults + "  " )
    if ( not os.path.isfile(config_userin)  ):
        bail("Cannot find configfile: " + config_userin + "  " )

    config = yaml.safe_load(open(config_defaults, 'r'))
    update_config( config, yaml.safe_load( open( config_userin, 'r')))

    config["execution"]["MM2_align_option"] = config["options"]["minimap2"][config["ReadType"]]

    if args.clustersub:
        # This is false by default, so over-writing will only happen
        # if user has set this to true
        config["execution"]["clustersub"] = args.clustersub

    if ( int(args.jobs) > 1 ):
        # Again, this will be 1 by default. Only over-written with active CL input
        config["execution"]["jobs"] = args.jobs

    if not args.target == "report":
        config["execution"]["target_out"] = args.target

    if len(set(config["samplelist"])) != len(config["samplelist"]):
        raise Exception("sampleIDs are not unique.")

    #TODO: trim possible "/" at the end of the sampleDir names
    for sample in config["samplelist"]:
       config["samplelist"][sample]["sampleDirs"] = {}

       for sampleDirName_i in range(0, len( config["samplelist"][sample]["sampleDirNames"])):
       #  --------- Populate the "chunk" list of reads for each sample : ---------
          sampleDirName = config["samplelist"][sample]["sampleDirNames"][sampleDirName_i]

          DATPATH = os.path.join( config["PATHIN"], sampleDirName, config["samplelist"][sample].get( "fast5dir", config["fast5dir_default"] ) )

          # i.e. the list of subdirectories within this path:
          unsorted_stringlist =  [entry.name for entry in os.scandir( DATPATH ) if entry.is_dir()]
          intlist = list( map(int, unsorted_stringlist) )
          intlist.sort()
          # we now have a list of the "chunks" of reads:

          # create the necessary substructure in config:
          config["samplelist"][sample]["sampleDirs"][ sampleDirName] = {}

          # and create the chunkdir list for each subdir
          config["samplelist"][sample]["sampleDirs"][ sampleDirName]["chunkdirlist"] = list( map(str, intlist))

          #  --------- Now get the hashID, fastq_prefix (if any) : ------------
          # locate the fastq files:
          fqdir  = os.path.join( config["PATHIN"], sampleDirName, config["samplelist"][sample].get( "fastqdir", config["fastqdir_default"] ) )

          # check if the fastq directory exists
          if ( os.path.isdir( fqdir )  ):
             # determine the common prefix that preceeds _{chunk}.fastq in each of these files:
             fastqprefix = list( set ( [ re.sub("_\d+" + config["fastq_suffix"], '_', entry.name ) for entry in  os.scandir( fqdir ) if entry.is_file() ] ) )
             if ( len( fastqprefix ) == 1 ):
                 config["samplelist"][sample][ "sampleDirs" ][sampleDirName]["fastq_prefix"] = fastqprefix[0]
             else:
                 bail("Number of fastq prefixes for sample" + sample +" in directory " + sampleDirName + " != 1 ")
    # FINISHED cycling through samples and prepping configs for each

    # If we are currently runing a SGE cluster submission, then
    # prepare a cluster-config file (specifying mem/time/etc. for each job)
    if( config["execution"]["clustersub"] ):
       generate_cluster_configuration( config )
       del config["execution"]["cluster"]["rules"]
    else:
       del config["execution"]["cluster"]
    # In either case, SMconfig file doesn't need cluster configuration details.

    #  --------- Now dump the config dictionary to the output path : ------------
    with open(config_npSM, 'w') as outfile:
        dumps = json.dumps( config,
                            indent=4, sort_keys=False,
                            separators=(",",": "),
                            ensure_ascii=True )
        outfile.write(dumps)

    return config
# --------------------------------------------------------------

def generate_cluster_configuration( config ):
    rules = config["execution"]["cluster"]["rules"]

    cluster_conf = {}
    for rule in rules:
        # Gather necessary resource requirements for each rule in the pipeline
        cluster_conf[rule] ={
                "nthreads": rules[rule].get( "threads", rules["__default__"]["threads"] ),
                "MEM":      rules[rule].get( "memory",  rules["__default__"]["memory"]  ),
                "h_stack":  rules[rule].get( "h_stack", rules["__default__"]["h_stack"] ),
                "queue":    rules[rule].get( "queue",   rules["__default__"].get("queue", "all") )
                }

    cluster_config_file = config["execution"]["cluster"]["cluster_config_file"]

    # and write them to an outfile:
    with open(cluster_config_file, 'w') as outfile:
        dumps = json.dumps(cluster_conf,
                           indent=4, sort_keys=True,
                           separators=(",",": "), ensure_ascii=True)
        outfile.write(dumps)
    return cluster_config_file

# --------------------------------------------------------------
def display_logo( splash_location ):
    if not os.getenv('NP_UGLY'):
        # Ensure that this is printed without errors, even if the user's
        # locale is not a UTF-8 locale.
        p = open(splash_location, 'r', encoding='utf-8').read()
        sys.stdout.buffer.write(p.encode('utf-8'))

# --------------------------------------------------------------

