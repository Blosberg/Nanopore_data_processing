import os, sys, json, csv, yaml
import re # regular expressions.
from snakemake.utils import update_config
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

# def makelink(src, target):
#     if not os.path.isfile(src):
#         bail("ERROR: Refusing to link non-existent file %s" % src)
#     elif not os.path.isdir(os.path.dirname(target)):
#         bail("%s or subdirectory does not exist for linking  " % target)
#     else:
#         try:
#             # os.symlink( "point_to_here",  "from_here" )
#             os.symlink(src, target)
#         except FileExistsError:
#             pass
# 
# --------------------------------------------------------------

def get_chunkfiles( samplename, DIR, prefix_string, suffix_string, quoted):
    Sample_indices_str     = config["samplelist"][samplename]["chunkdirlist"]
    FILES_list = [ os.path.join( DIR, prefix_string + samplename + "_" + chunk + suffix_string ) for chunk in Sample_indices_str ]
    if ( quoted ):
       return( ",".join( FILES_list )  )
    else:
       return( FILES_list )

# --------------------------------------------------------------

def prep_configfile( config_defaults, config_userin, config_npSM_out):
    # Create the SM config file, from default, and user-defined inputs:

    config = yaml.safe_load(open(config_defaults, 'r'))
    update_config( config, yaml.safe_load( open( config_userin, 'r')))

    if len(set(config["samplelist"])) != len(config["samplelist"]):
        raise Exception("sampleIDs are not unique.")

    for sample in config["samplelist"]:

       #  --------- Populate the "chunk" list of reads for each sample : ---------
       DATPATH = os.path.join( config["PATHIN"], config["samplelist"][sample]["sampledir"], config["samplelist"][sample].get( "fast5dir", config["fast5dir_default"] ) ) 

       # i.e. the list of subdirectories within this path:
       unsorted_stringlist =  [entry.name for entry in os.scandir( DATPATH ) if entry.is_dir()]
       intlist = list( map(int, unsorted_stringlist) )
       intlist.sort()

       # we now have a list of the "chunks" of reads:
       config["samplelist"][sample]["chunkdirlist"] = list( map(str, intlist))

       #  --------- Now get the hashID, fastq_prefix (if any) : ------------
       # locate the files:
       fqdir  = os.path.join( config["PATHIN"], config["samplelist"][sample]["sampledir"], config["samplelist"][sample].get( "fastqdir", config["fastqdir_default"] ) )
       if ( os.path.isdir( fqdir )  ):
          # determine the common prefix that preceeds _{chunk}.fastq in each of these files:
          fastqprefix = list( set (  [ re.sub("_\d+." + config["fastq_suffix"], '_', entry.name ) for entry in  os.scandir( fqdir ) if entry.is_file() ] ) )
          if ( len( fastqprefix ) != 1 ):
              bail("Number of fastq prefixes for a given sample is different from 1 ")
          else:
              config["samplelist"][sample]["fastq_prefix"] = fastqprefix[0]

    # --------- Now dump the config dictionary to the output path : ------------
    if( config["execution"]["submit-to-cluster"] ):
       generate_cluster_configuration( config_defaults )

    #  --------- Now dump the config dictionary to the output path : ------------

    # SMconfig file doesn't need these.
    del config["execution"]["cluster"]
    del config["execution"]["rules"]

    with open(config_npSM_out, 'w') as outfile:
        dumps = json.dumps( config,
                            indent=4, sort_keys=False,
                            separators=(",",": "),
                            ensure_ascii=True )
        outfile.write(dumps)

    return 0

# --------------------------------------------------------------

def generate_cluster_configuration(config_defaults):
    rules = config['execution']['rules']

    cluster_conf = {}
    for rule in rules:
        if 'queue' in config['execution']['cluster']:
            # --- User has supplied general queue name for all rules---
            if 'queue' in rules[rule]:
              bail("ERROR: 'queue' multiply defined in settings file.") #AND per rule ->error
            else:
              cluster_conf[rule] = {
              'nthreads': rules[rule]['threads'],
              'MEM':      rules[rule]['memory'],
              'queue':    config['execution']['cluster']['queue'],
              'h_stack':  config['execution']['cluster']['stack']
              }
        elif ( 'queue' in rules[rule] ):
            # --- User has supplied a queue for this specific rule.
            if not 'queue' in rules['__default__']:
              bail("ERROR: submission queue specified per rule with no default.")
            cluster_conf[rule] = {
              'nthreads': rules[rule]['threads'],
              'MEM':      rules[rule]['memory'],
              'queue':    rules[rule]['queue'],
              'h_stack':  config['execution']['cluster']['stack']
              }
        else:
              # --- User has provided no information on queue for this rule -> default.
              cluster_conf[rule] = {
              'nthreads': rules[rule]['threads'],
              'MEM':      rules[rule]['memory'],
              'h_stack':  config['execution']['cluster']['stack']
              }

    cluster_config_file = "cluster_conf.json"
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

