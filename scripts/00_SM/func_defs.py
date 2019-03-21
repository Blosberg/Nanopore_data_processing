def bail( msg ):
    """Print the error message to stderr and exit."""
    print("ERROR: " + msg + "... exiting.", file=sys.stderr)
    exit(1)

def dumpjson( yamldatin, fout ):
    with open(fout, 'w') as outfile:
        dumps = json.dumps(yamldatin,
                           indent=4, sort_keys=True,
                           separators=(",",": "), ensure_ascii=True)
        outfile.write(dumps)


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
            # os.symlink( "point_to_here",  "from_here" )
            os.symlink(src, target)
        except FileExistsError:
            pass

def get_chunkfiles( samplename, DIR, prefix_string, suffix_string, quoted):
    Sample_indices_str     = config["samplelist"][samplename]["chunkdirlist"]
    FILES_list = [ os.path.join( DIR, prefix_string + samplename + "_" + chunk + suffix_string ) for chunk in Sample_indices_str ]
    if ( quoted ):
       return( ",".join( FILES_list )  )
    else:
       return( FILES_list )
