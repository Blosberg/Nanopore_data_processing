# 
rule create_currentGR_obj:
# Dummy rule to produce empty GRanges file. make_report will check the size
# of this file, and if it sees that it is empty, the section will be skipped. 
    input:
        sortedbam   = os.path.join( DIR_SORTED_MINIMAPPED, "{sample}.sorted.bam")
    output:
        GRobj       = os.path.join( DIR_GR, "{sample}_GR.RData")
    log:
        os.path.join( DIR_GR, "{sample}_GR_conversion.log")
    message: fmt("Creating empty (dummy) GRanges file for wholefastq stream")
    shell:
        "touch {output} && echo \"Created empty GRobj file for wholefastq stream\" > {log} "  

# -----------------------------------------------------

rule convert_sort_minimap:
# convert from sam to bam format and sort by position
    input:
        aligned     = os.path.join( DIR_FILTERED_MINIMAP, "{sample}.0filtered.sam")
    output:
        sortedbam   = os.path.join( DIR_SORTED_MINIMAPPED, "{sample}.sorted.bam")
    params:
        options = "-ax splice "
    log:
        logfile = os.path.join( DIR_SORTED_MINIMAPPED, "read_chunks", "{sample}.sortbam.log")
    message: 
        """ --- converting, sorting, and indexing bam file. --- """
    shell:
        'samtools view  -Sb  {input} | samtools sort > {output} && samtools index {output} 2> {log.logfile}'

#------------------------------------------------------

rule filter_nonaligned_minimap:
# Check for alignment filter in sam file: if != 4 then remove this read
    input:
        aligned  = os.path.join( DIR_ALIGNED_MINIMAP, "{sample}.sam" )
    output:
        aligned  = os.path.join( DIR_FILTERED_MINIMAP, "{sample}.0filtered.sam" )
    log:
        log      = os.path.join( DIR_FILTERED_MINIMAP, "{sample}.0filtering.log" )
    message: 
        """--- filtering unaligned reads from alignment data ---"""
    shell:
        " cat {input} | perl -lane 'print if $F[1] ne 4'  >  {output}   2> {log}"

#------------------------------------------------------

rule align_minimap:
# use minimap2 to align the fastq reads to the reference genome
    input:
        mmiref   = os.path.join( config['ref']['Genome_DIR'] , config['ref']['Genome_version']+ ".mmi" ),
        sample   = getPathCase( config['PATHIN'], "fastq", "pass", "{sample}.fq.gz", intype )
    output:
        aligned  = os.path.join( DIR_ALIGNED_MINIMAP, "{sample}.sam" )
    params:
        options  = " -ax splice "
    log:
        log      = os.path.join( DIR_ALIGNED_MINIMAP, "{sample}_alignment.log")
    message: 
        """--- aligning fastq reads to indexed reference"""
    shell:
        "{MM2} {params} {input.mmiref} {input.sample} > {output}  2> {log}"

#------------------------------------------------------

rule minimizer:
# Create indexed version of reference genome for fast alignment with minimap2 later:
    input:
        refgenome_fasta  = os.path.join(config['ref']['Genome_DIR'] , config['ref']['Genome_version']+ ".fa" )
    output:
        refgenome_mmiref = os.path.join(config['ref']['Genome_DIR'] , config['ref']['Genome_version']+ ".mmi")
    params:
        options = " -d  "
    log:
        os.path.join( config['ref']['Genome_DIR'], config['ref']['Genome_version'], "_mmi2_minimizer_creation.log")
    message: 
        """--- creating minimizer index of reference genome for minimap2."""
    shell:
        "{MM2} {params.options}  {output} {input} 2> {log}"

