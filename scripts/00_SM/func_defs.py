def bail( msg ):
    """Print the error message to stderr and exit."""
    print("ERROR: " + msg + "... exiting.", file=sys.stderr)
    exit(1)

# Generate a command line string that can be passed to snakemake's
# "shell".  The string is prefixed with an invocation of "nice".
def tool(name):
    return config['tools'][name]['executable']

def nice(cmd, args, log=None):
    # executable = tool(cmd)
    # line = ["nice -" + str(config['execution']['nice']), executable] + [toolArgs(cmd)] + args
    line = ["nice -" + str(config['execution']['nice']), cmd] + args
    if log:
        line.append("> {} 2>&1".format(log))
    return " ".join(line)

def fmt(message):
    """Format the MESSAGE string."""
    return "----------  " + message + "  ----------"

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

def makelink(src, target):
    if not os.path.isfile(src):
        bail("ERROR: Refusing to link non-existent file %s" % src)
    elif not os.path.isdir(os.path.dirname(target)):
        bail("%s or subdirectory does not exist for linking  " % target)
    else:
        try:
            os.symlink(src, target)
        except FileExistsError:
            pass


def get_chunkfiles( sample, DIR, prefix_string, suffix_string, quoted):
 
    Sample_indices_str     = config["samplelist"][sample]["chunkdirlist"] 

    FILES_list = list( chain( *[ expand ( os.path.join( DIR, prefix_string + "_" + sample + "_" + chunk + suffix_string ), ) for chunk in Sample_indices_str ] ) )

    if ( quoted ):
       return( ",".join( FILES_list )  )
    else:
       return( FILES_list ) 
