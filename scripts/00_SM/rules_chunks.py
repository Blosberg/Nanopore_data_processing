# -----------------------------------------------------
# BEGIN RULES
# -----------------------------------------------------


# -----------------------------------------------------
# Take a read-separated list of events, and assemble a histogram of
# current values for each unique kmer observed.
rule create_kmer_histlist:
    input:
        RDS_GRLreads       = os.path.join( config["PATHOUT"], "{wckmerhist_sampleDir}", SUBDIR_GR, "{wckmerhist_samplename}_reads_GRL.rds")
    output:
        RDS_histlist       = os.path.join( config["PATHOUT"], "{wckmerhist_sampleDir}", SUBDIR_GR, "{wckmerhist_samplename}_kmer_histlist.rds")
    params:
        current_histmin=50,
        current_histmax=150,
        current_histres=0.5,
        RDS_GRLreads_in=os.path.join( config["PATHOUT"], "{wckmerhist_sampleDir}", SUBDIR_GR, "{wckmerhist_samplename}_reads_GRL.rds"),
        RDS_histlist_out=os.path.join( config["PATHOUT"], "{wckmerhist_sampleDir}", SUBDIR_GR, "{wckmerhist_samplename}_kmer_histlist.rds"),
        k=5
    log:
        logfile=os.path.join( config["PATHOUT"], "{wckmerhist_sampleDir}", SUBDIR_GR, "{wckmerhist_samplename}_histlist.log")
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
        GRobj_in     = os.path.join( config["PATHOUT"], "{wcflatten_sampleDir}", SUBDIR_GR, "{wcflatten_samplename}_reads_GRL.rds")
    output:
        GRflat_out   = os.path.join( config["PATHOUT"], "{wcflatten_sampleDir}", SUBDIR_GR, "{wcflatten_samplename}_reads_flat_GRL.rds")
    params:
        Rfuncs_flatten_reads = R_flattenreads_funcs,
        readsGRL_in          = os.path.join( config["PATHOUT"], "{wcflatten_sampleDir}", SUBDIR_GR, "{wcflatten_samplename}_reads_GRL.rds"),
        flatreadsGRL_out     = os.path.join( config["PATHOUT"], "{wcflatten_sampleDir}", SUBDIR_GR, "{wcflatten_samplename}_reads_flat_GRL.rds"),
        samplename           = "{wcflatten_samplename}"
    log:
        os.path.join( config["PATHOUT"], "{wcflatten_sampleDir}", SUBDIR_GR, "{wcflatten_samplename}_flattenreads_GRL.log")
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
# from the .tsv files produced by eventalign
# each entry of the list is a read:

rule create_readcurrent_GRL_obj:
    input:
        tsvfile      = lambda wc: get_chunkfiles( wc.wcreadGRL_samplename, os.path.join( config["PATHOUT"], wc.wcreadGRL_sampleDir, SUBDIR_EVENTALIGN, "tsv_chunks" ), "Ealign_", ".tsv", False )
    output:
        GRobj        = os.path.join( config["PATHOUT"], "{wcreadGRL_sampleDir}", SUBDIR_GR, "{wcreadGRL_samplename}_reads_GRL.rds"),
        poremodel    = os.path.join( config["PATHOUT"], "{wcreadGRL_sampleDir}", SUBDIR_GR, "{wcreadGRL_samplename}_poremodel.tsv")
    params:
        Rfuncs_tsv2GRconv = R_tables2GR_funcs,
        output_reads_GRL  = os.path.join( config["PATHOUT"], "{wcreadGRL_sampleDir}", SUBDIR_GR, "{wcreadGRL_samplename}_reads_GRL.rds"),
        output_poremodel  = os.path.join( config["PATHOUT"], "{wcreadGRL_sampleDir}", SUBDIR_GR, "{wcreadGRL_samplename}_poremodel.tsv"),
        Flatten           = config["Flatten"],
        Ealign_files      = lambda wc: get_chunkfiles( wc.wcreadGRL_samplename, os.path.join( config["PATHOUT"], wc.wcreadGRL_sampleDir, SUBDIR_EVENTALIGN, "tsv_chunks" ), "Ealign_", ".tsv", True ),
        samplename          = "{wcreadGRL_samplename}"
    log:
        os.path.join( config["PATHOUT"], "{wcreadGRL_sampleDir}", SUBDIR_GR, "{wcreadGRL_samplename}_reads_GRL_conversion.log")
    message:
        fmt("Convert aligned NP reads to GRanges object")
    shell:
        nice('Rscript', [ R_tables2GR_main,
                         "--Rfuncs_tsv2GRconv={params.Rfuncs_tsv2GRconv}",
                         "--output_reads_GRL={params.output_reads_GRL}",
                         "--output_poremodel={params.output_poremodel}",
                         "--samplename={params.samplename}",
                         "--Flatten={params.Flatten}",
                         "--logFile={log}",
                         "--Ealign_files={params.Ealign_files}"] )

# -----------------------------------------------------

rule np_index:
    # Index the reads and the fast5 files for nanopolish
    input:
        fast5_folder = os.path.join( config["PATHIN"], "{wcnpindex_sampleDir}/fast5/pass", "{wcnpindex_chunk}"),
#        fast5_folder = lambda wc: getPathCase( os.path.join( config["PATHIN"], "{wcnpindex_sampleDir}"), 'fast5', 'pass', wc.wcnpindex_chunk, input_data_type ),
        fastq_file   = os.path.join( config["PATHOUT"], "{wcnpindex_sampleDir}", SUBDIR_SYMLINKS, "{wcnpindex_samplename}_{wcnpindex_chunk}." + config["fastq_suffix"])
#        fastq_file   = lambda wc: os.path.join( config["PATHOUT"], "{wcnpindex_sampleDir}", SUBDIR_SYMLINKS, wc.wcnpindex_samplename + "_{wcnpindex_chunk}." + config["fastq_suffix"] )
    output:
        npi    = os.path.join( config["PATHOUT"], "{wcnpindex_sampleDir}", SUBDIR_SYMLINKS, "{wcnpindex_samplename}_{wcnpindex_chunk}."+config["fastq_suffix"]+".index"  ),
        fai    = os.path.join( config["PATHOUT"], "{wcnpindex_sampleDir}", SUBDIR_SYMLINKS, "{wcnpindex_samplename}_{wcnpindex_chunk}."+config["fastq_suffix"]+".index.fai" ),
        gzi    = os.path.join( config["PATHOUT"], "{wcnpindex_sampleDir}", SUBDIR_SYMLINKS, "{wcnpindex_samplename}_{wcnpindex_chunk}."+config["fastq_suffix"]+".index.gzi" ),
        readdb = os.path.join( config["PATHOUT"], "{wcnpindex_sampleDir}", SUBDIR_SYMLINKS, "{wcnpindex_samplename}_{wcnpindex_chunk}."+config["fastq_suffix"]+".index.readdb")
    params:
        options    = " index -d "
    log:
        logfile  = os.path.join( config["PATHOUT"], "{wcnpindex_sampleDir}", SUBDIR_SYMLINKS,  "{wcnpindex_samplename}_{wcnpindex_chunk}_npi.log" )
    message: """---- Index the reads from chunk {wildcards.wcnpindex_chunk} against the fast5 files from the same. ----"""
    shell:
        " nice -19 {nanopolish} {params.options} {input.fast5_folder} {input.fastq_file} 2> {log.logfile} "


# -----------------------------------------------------
rule merge_bam_files:
    # combine the various ~4000 read chunks into a single bam file
    input:
        bami      = lambda wc: get_chunkfiles( wc.wcmerge_samplename, os.path.join( config["PATHOUT"], wc.wcmerge_sampleDir, SUBDIR_SORTED_MINIMAPPED, "bam_chunks"), "" , ".sorted.bam", False )
    output:
        sortedbam = os.path.join( config["PATHOUT"], "{wcmerge_sampleDir}", SUBDIR_SORTED_MINIMAPPED, "{wcmerge_samplename}.sorted.bam")
    log:
        logfile   = os.path.join( config["PATHOUT"], "{wcmerge_sampleDir}", SUBDIR_SORTED_MINIMAPPED, "{wcmerge_samplename}.mergingbam.log")
    message:
        """ --- combining bam files from post-mapping fastq data. --- """
    shell:
        '{SAMTOOLS} merge {output} {input} '

#------------------------------------------------------

# convert minimap2 alignments from sam to bam format
# and sort by position
rule convert_sort_minimap:
    input:
        aligned     = os.path.join( config["PATHOUT"], "{wcsort_sampleDir}", SUBDIR_FILTERED_MINIMAP, "{wcsort_samplename}_{wcsort_chunk}.0filtered.sam")
    output:
        sortedbam   = os.path.join( config["PATHOUT"], "{wcsort_sampleDir}", SUBDIR_SORTED_MINIMAPPED, "bam_chunks", "{wcsort_samplename}_{wcsort_chunk}.sorted.bam")
    params:
        options = "-ax splice "
    log:
        logfile = os.path.join( config["PATHOUT"], "{wcsort_sampleDir}", SUBDIR_SORTED_MINIMAPPED, "bam_chunks", "{wcsort_samplename}_{wcsort_chunk}.sortingbam.log")
    message:
        """ --- converting, sorting, and indexing bam file. --- """
    shell:
        'samtools view  -Sb  {input} | samtools sort > {output} && samtools index {output} 2> {log.logfile}'


#------------------------------------------------------

# Check for alignment filter in minimap2s sam file:
#  if != 4 then remove this read
rule filter_nonaligned_minimap:
    input:
        aligned  = os.path.join( config["PATHOUT"], "{wcfilter_sampleDir}", SUBDIR_ALIGNED_MINIMAP, "{wcfilter_samplename}_{wcfilter_chunk}.sam" )
    output:
        filtered = os.path.join( config["PATHOUT"], "{wcfilter_sampleDir}", SUBDIR_FILTERED_MINIMAP, "{wcfilter_samplename}_{wcfilter_chunk}.0filtered.sam" )
    log:
        log      = os.path.join( config["PATHOUT"], "{wcfilter_sampleDir}", SUBDIR_FILTERED_MINIMAP, "{wcfilter_samplename}_{wcfilter_chunk}.0filtering.log" )
    message:
        """--- filtering unaligned reads from alignment data ---"""
    shell:
        " cat {input} | perl -lane 'print if $F[1] ne 4'  >  {output}   2> {log}"

#------------------------------------------------------

rule align_minimap:
    # use minimap2 to align the fastq reads to the reference
    # genome
    input:
        mmiref   = os.path.join( DIR_REFGENOME , config['ref']['Genome_version']+ ".mmi" ),
        fqlink   = os.path.join( config["PATHOUT"], "{wcalign_sampleDir}", SUBDIR_SYMLINKS, "{wcalign_samplename}_{wcalign_chunk}." + config["fastq_suffix"] )
    output:
        aligned  = os.path.join( config["PATHOUT"], "{wcalign_sampleDir}", SUBDIR_ALIGNED_MINIMAP,  "{wcalign_samplename}_{wcalign_chunk}.sam" )
    params:
        options  = " -ax splice -uf -k14"
    log:
        log      = os.path.join( config["PATHOUT"], "{wcalign_sampleDir}", SUBDIR_ALIGNED_MINIMAP,  "{wcalign_samplename}_{wcalign_chunk}_alignment.log" )
    message:
        """--- aligning fastq reads from {input.fqlink} to {input.mmiref}; outputting to {output} --- """
    shell:
        "{MM2} {params} {input.mmiref} {input.fqlink} > {output}  2> {log}"

#------------------------------------------------------
rule create_fqsymlink:
    # create symlink to data:
    input:
        chunk_src     = lambda wc: os.path.join( config["PATHIN"], wc.wcfqlink_sampleDir, "fastq", "pass", config["samplelist"][wc.wcfqlink_samplename]["fastq_prefix"] + wc.wcfqlink_chunk + "." + config["fastq_suffix"] )
    output:
        chunk_linkloc = os.path.join( config["PATHOUT"], "{wcfqlink_sampleDir}", SUBDIR_SYMLINKS, "{wcfqlink_samplename}_{wcfqlink_chunk}." + config["fastq_suffix"] )
    params:
        options  = " -s "
    message:
        """--- creating symbolic link to chunk data ---"""
    shell:
        " ln {params} {input} {output}"

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

