# --- enumerate "chunk" files list  
# blocks that the minion divides the ouput into 
# (usually 4000 reads per chunk) 
# This is done for both event alignments and bam files:

Sample_indices_int     = range(config["MAXSAMPLEi"] +1) 
Sample_indices_str     = [ str(item) for item in Sample_indices_int  ]

Ealign_FILES_list = list( chain( *[ expand ( os.path.join( DIR_EVENTALIGN, 'Ealign_'+chunk+'.cvs' ), ) for chunk in Sample_indices_str ] ) )

Ealign_FILES_quoted   = ",".join( Ealign_FILES_list )  
 
bami_FILES_list = list( chain( *[ expand ( os.path.join( DIR_SORTED_MINIMAPPED, "read_chunks", 'run_'+config["RUN_ID"]+ '_'+chunk+'.sorted.bam' ), ) for chunk in Sample_indices_str ] ) )


# -----------------------------------------------------
# BEGIN RULES
# -----------------------------------------------------


rule create_currentGR_obj:
# produce GRanges object of current data from reads using .RData format 
    input:
        table_files  = Ealign_FILES_list
    output:
        GRobj        = os.path.join( DIR_GR, "{sample}_GR.RData")
    params:
        Rfuncs_file  = os.path.join( config[ "scripts"]["script_folder"], config[ "scripts"]["Rfuncs_file"] ), 
        output       = os.path.join( DIR_GR, "{sample}_GR.RData"),
        Ealign_files = Ealign_FILES_quoted 
    log:
        os.path.join( DIR_GR, "{sample}_GR_conversion.log")
    message: fmt("Convert aligned NP reads to GRanges object")
    shell:
        nice('Rscript', ["./scripts/npreads_tables2GR.R",
                         "--Rfuncs_file={params.Rfuncs_file}",
                         "--output={params.output}",
                         "--logFile={log}",
                         "--Ealign_files={params.Ealign_files}"]
)

# -----------------------------------------------------

rule merge_bam_files:
# combine the ~4000 reads from each iteration of the NP data into a single bam file
    input:
        bami      = bami_FILES_list
    output:
        sortedbam = os.path.join( DIR_SORTED_MINIMAPPED, "run_{sample}.sorted.bam")
    log:
        logfile   = os.path.join( DIR_SORTED_MINIMAPPED, "{sample}.sortbam.log")
    message: 
        """ --- combining bam files from post-mapping fastq data. --- """
    shell:
        '{SAMTOOLS} merge {output} {input} '

#------------------------------------------------------

# THIS SHOULD BE THE LEAF NODE WHEN TARGET=.BAM
rule convert_sort_minimap:
# convert from sam to bam format and sort by position
    input:
        aligned     = os.path.join( DIR_FILTERED_MINIMAP, "run_{sample}_{index}.0filtered.sam")
    output:
        sortedbam   = os.path.join( DIR_SORTED_MINIMAPPED, "read_chunks", "run_{sample}_{index}.sorted.bam")
    params:
        options = "-ax splice "
    log:
        logfile = os.path.join( DIR_SORTED_MINIMAPPED, "read_chunks", "run_{sample}_{index}.sortbam.log")
    message: 
        """ --- converting, sorting, and indexing bam file. --- """
    shell:
        'samtools view  -Sb  {input} | samtools sort > {output} && samtools index {output} 2> {log.logfile}'

#------------------------------------------------------

rule filter_nonaligned_minimap:
# Check for alignment filter in sam file: if != 4 then remove this read
    input:
        aligned  = os.path.join( DIR_ALIGNED_MINIMAP, "run_{sample}_{index}.sam" )
    output:
        aligned  = os.path.join( DIR_FILTERED_MINIMAP, "run_{sample}_{index}.0filtered.sam" )
    log:
        log      = os.path.join( DIR_FILTERED_MINIMAP, "run_{sample}_{index}.0filtering.log" )
    message: 
        """--- filtering unaligned reads from alignment data ---"""
    shell:
        " cat {input} | perl -lane 'print if $F[1] ne 4'  >  {output}   2> {log}"

#------------------------------------------------------

rule align_minimap:
# use minimap2 to align the fastq reads to the reference genome
    input:
        mmiref   = os.path.join( DIR_REFGEMONE , config['ref']['Genome_version']+ ".mmi" ),
        sample   = os.path.join( config['PATHIN'], "fastq", "pass", "fastq_runid_{sample}_{index}.fastq" )
    output:
        aligned  = os.path.join( DIR_ALIGNED_MINIMAP, "run_{sample}_{index}.sam" )
    params:
        options  = " -ax splice "
    log:
        log      = os.path.join( DIR_ALIGNED_MINIMAP, "run_{sample}_{index}_alignment.log")
    message: 
        """--- aligning fastq reads to indexed reference"""
    shell:
        "{MM2} {params} {input.mmiref} {input.sample} > {output}  2> {log}"

#------------------------------------------------------

rule minimizer:
# Create indexed version of reference genome for fast alignment with minimap2 later:
    input:
        refgenome_fasta  = os.path.join(DIR_REFGEMONE , config['ref']['Genome_version']+ ".fa" )
    output:
        refgenome_mmiref = os.path.join(DIR_REFGEMONE , config['ref']['Genome_version']+ ".mmi")
    params:
        options = " -d  "
    log:
        os.path.join( DIR_REFGEMONE, config['ref']['Genome_version'], "_mmi2_minimizer_creation.log")
    message: 
        """--- creating minimizer index of reference genome for minimap2."""
    shell:
        "{MM2} {params.options}  {output} {input} 2> {log}"

