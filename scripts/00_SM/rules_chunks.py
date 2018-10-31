# -----------------------------------------------------
# BEGIN RULES
# -----------------------------------------------------
# table2GR_script = os.path.join( [config[ "scripts"]["scripts_folder"], [config[ "scripts"]["Rfuncs_tableGRconv_file"] )

# produce GRanges object of current data in .RData 
# format from the .csv files produced by eventalign
rule create_currentGR_obj:
    input:
        csvfile      = lambda wc: get_chunkfiles( wc.sample, os.path.join( DIR_EVENTALIGN, "csv_chunks" ), "Ealign", ".csv", False )
    output:
        GRobj        = os.path.join( DIR_GR, "{sample}_GR.RData")
    params:
        Rfuncs_tableGRconv_file  = os.path.join( config[ "scripts"]["script_folder"], config[ "scripts"]["Rfuncs_tableGRconv_file"] ), 
        output       = os.path.join( DIR_GR, "{sample}_GR.RData"),
        Ealign_files = lambda wc: get_chunkfiles( wc.sample, os.path.join( DIR_EVENTALIGN, "csv_chunks" ), "Ealign", ".csv", True ), 
        samplename   = "{sample}"
    log:
        os.path.join( DIR_GR, "{sample}_GR_conversion.log")
    message:
        fmt("Convert aligned NP reads to GRanges object")
    shell:
        nice('Rscript', [ tables2GR_main,
                         "--Rfuncs_tableGRconv_file={params.Rfuncs_tableGRconv_file}",
                         "--output={params.output}",
                         "--logFile={log}",
                         "--samplename={params.samplename}",
                         "--Ealign_files={params.Ealign_files}"] )


# -----------------------------------------------------

# combine the ~4000 reads from each minimap2 
# alignment into a single bam file
rule merge_bam_files:
    input:
        bami      = lambda wc: get_chunkfiles( wc.sample, os.path.join(DIR_SORTED_MINIMAPPED, "bam_chunks"), "run" , ".sorted.bam", False )
    output:
        sortedbam = os.path.join( DIR_SORTED_MINIMAPPED, "run_{sample}.sorted.bam")
    log:
        logfile   = os.path.join( DIR_SORTED_MINIMAPPED, "{sample}.sortbam.log")
    message: 
        """ --- combining bam files from post-mapping fastq data. --- """
    shell:
        '{SAMTOOLS} merge {output} {input} '

#------------------------------------------------------

# convert minimap2 alignments from sam to bam format 
# and sort by position
rule convert_sort_minimap:
    input:
        aligned     = os.path.join( DIR_FILTERED_MINIMAP, "run_{sample}_{chunk}.0filtered.sam")
    output:
        sortedbam   = os.path.join( DIR_SORTED_MINIMAPPED, "bam_chunks", "run_{sample}_{chunk}.sorted.bam")
    params:
        options = "-ax splice "
    log:
        logfile = os.path.join( DIR_SORTED_MINIMAPPED, "bam_chunks", "run_{sample}_{chunk}.sortbam.log")
    message: 
        """ --- converting, sorting, and indexing bam file. --- """
    shell:
        'samtools view  -Sb  {input} | samtools sort > {output} && samtools index {output} 2> {log.logfile}'


#------------------------------------------------------

# Check for alignment filter in minimap2s sam file:
#  if != 4 then remove this read
rule filter_nonaligned_minimap:
    input:
        aligned  = os.path.join( DIR_ALIGNED_MINIMAP, "run_{sample}_{chunk}.sam" )
    output:
        aligned  = os.path.join( DIR_FILTERED_MINIMAP, "run_{sample}_{chunk}.0filtered.sam" )
    log:
        log      = os.path.join( DIR_FILTERED_MINIMAP, "run_{sample}_{chunk}.0filtering.log" )
    message: 
        """--- filtering unaligned reads from alignment data ---"""
    shell:
        " cat {input} | perl -lane 'print if $F[1] ne 4'  >  {output}   2> {log}"


#------------------------------------------------------

# use minimap2 to align the fastq reads to the reference 
# genome
rule align_minimap:
    input:
        mmiref   = os.path.join( DIR_REFGENOME , config['ref']['Genome_version']+ ".mmi" ),
        sample   = lambda wc: os.path.join( DIR_SYMLINKS,   config['samplelist'][sample]["RUN_ID"] + "_" + str(wc.chunk) + config['samplelist'][sample]["fastq_suffix"] )
    output:
        aligned  = os.path.join( DIR_ALIGNED_MINIMAP, "run_{sample}_{chunk}.sam" )
    params:
        options  = " -ax map-ont "
    log:
        log      = os.path.join( DIR_ALIGNED_MINIMAP, "run_{sample}_{chunk}_alignment.log")
    message: 
        """--- aligning fastq reads to indexed reference"""
    shell:
        "{MM2} {params} {input.mmiref} {input.sample} > {output}  2> {log}"

#------------------------------------------------------

# Create indexed version of reference genome for fast 
# alignment with minimap2 later:
rule minimizer:
    input:
        refgenome_fasta  = os.path.join(DIR_REFGENOME , config['ref']['Genome_version']+ ".fa" )
    output:
        refgenome_mmiref = os.path.join(DIR_REFGENOME , config['ref']['Genome_version']+ ".mmi")
    params:
        options = " -d  "
    log:
        os.path.join( DIR_REFGENOME, config['ref']['Genome_version'], "_mmi2_minimizer_creation.log")
    message: 
        """--- creating minimizer index of reference genome for minimap2."""
    shell:
        "{MM2} {params.options}  {output} {input} 2> {log}"

