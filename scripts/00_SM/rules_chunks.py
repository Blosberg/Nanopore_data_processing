# -----------------------------------------------------
# BEGIN RULES
# -----------------------------------------------------


# -----------------------------------------------------
# Take a read-separated list of events, and assemble a histogram of
# current values for each unique kmer observed.
rule create_kmer_histlist:
    input:
        RDS_GRLreads       = os.path.join( DIR_GR, "{sample}_reads_GRL.rds")
    output:
        RDS_histlist       = os.path.join( DIR_GR, "{sample}_kmer_histlist.rds")
    params:
        current_histmin=50,
        current_histmax=150,
        current_histres=0.5,
        RDS_GRLreads_in=os.path.join( DIR_GR, "{sample}_reads_GRL.rds"),
        RDS_histlist_out=os.path.join( DIR_GR, "{sample}_kmer_histlist.rds"),
        k=5
    log:
        logfile=os.path.join( DIR_GR, "{sample}_histlist.log")
    message:
        fmt("Build list of histograms for unique kmers in dataset.")
    shell:
        nice('Rscript', [ R_build_histlist_main,
                          "--current_histmin={params.current_histmin}",
                          "--current_histmax={params.current_histmax}",
                          "--current_histres={params.current_histres}",
                          "--rds_fin_readdat={params.RDS_GRLreads_in}",
                          "--rds_fout_histlist={params.RDS_histlist_out}",
                          "--k={params.k}",
                          "--logfile={log.logfile}",] )

# -----------------------------------------------------
# flatten reads from GRangesList in .RDS format
# by "flatten", we mean collaps multiple 'events' assigned to the same segment.

rule flatten_reads:
    input:
        GRobj_in     = os.path.join( DIR_GR, "{sample}_reads_GRL.rds")
    output:
        GRflat_out   = os.path.join( DIR_GR, "{sample}_reads_flat_GRL.rds")
    params:
        Rfuncs_flatten_reads = R_flattenreads_funcs,
        readsGRL_in          = os.path.join( DIR_GR, "{sample}_reads_GRL.rds"),
        flatreadsGRL_out    = os.path.join( DIR_GR, "{sample}_reads_flat_GRL.rds"),
        samplename           = "{sample}"
    log:
        os.path.join( DIR_GR, "{sample}_flattenreads_GRL.log")
    message:
        fmt("Flatten intra-read events with coincident alignments")
    shell:
        nice('Rscript', [ R_flattenreads_main,
                         "--Rfuncs_flattenreads={params.Rfuncs_flatten_reads}",
                         "--readsGRL_in={params.readsGRL_in}",
                         "--flatreadsGRL_out={params.flatreadsGRL_out}",
                         "--logFile={log}",
                         "--samplename={params.samplename}"] )

# -----------------------------------------------------
# produce GRangesList object of current data in .RDS format
# from the .csv files produced by eventalign
# each entry of the list is a read:

rule create_readcurrent_GRL_obj:
    input:
        csvfile      = lambda wc: get_chunkfiles( wc.sample, os.path.join( DIR_EVENTALIGN, "csv_chunks" ), "Ealign", ".csv", False )
    output:
        GRobj        = os.path.join( DIR_GR, "{sample}_reads_GRL.rds"),
        poremodel    = os.path.join( DIR_GR, "{sample}_poremodel.tsv")
    params:
        Rfuncs_table2GRconv = R_tables2GR_funcs,
        output_reads_GRL    = os.path.join( DIR_GR, "{sample}_reads_GRL.rds"),
        output_poremodel    = os.path.join( DIR_GR, "{sample}_poremodel.tsv"),
        Ealign_files        = lambda wc: get_chunkfiles( wc.sample, os.path.join( DIR_EVENTALIGN, "csv_chunks" ), "Ealign", ".csv", True ),
        samplename          = "{sample}"
    log:
        os.path.join( DIR_GR, "{sample}_reads_GRL_conversion.log")
    message:
        fmt("Convert aligned NP reads to GRanges object")
    shell:
        nice('Rscript', [ R_tables2GR_main,
                         "--Rfuncs_table2GRconv={params.Rfuncs_table2GRconv}",
                         "--output_reads_GRL={params.output_reads_GRL}",
                         "--output_poremodel={params.output_poremodel}",
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
        sample   = lambda wc: os.path.join( DIR_SYMLINKS,   sample + "_" + str(wc.chunk) + config['samplelist'][sample]["fastq_suffix"] )
    output:
        aligned  = os.path.join( DIR_ALIGNED_MINIMAP, "run_{sample}_{chunk}.sam" )
    params:
        options  = " -ax splice -uf -k14"
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

