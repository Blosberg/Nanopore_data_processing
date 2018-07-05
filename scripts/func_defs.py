# Generate a command line string that can be passed to snakemake's
# "shell".  The string is prefixed with an invocation of "nice".
def nice(cmd, args, log=None):
    executable = tool(cmd)
    line = ["nice -" + str(config['execution']['nice']), executable] + [toolArgs(cmd)] + args
    if log:
        line.append("> {} 2>&1".format(log))
    return " ".join(line)

def fmt(message):
    """Format the MESSAGE string."""
    return "----------  " + message + "  ----------"


