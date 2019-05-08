# -----------------------------------------------------
# BEGIN RULES
# -----------------------------------------------------


rule np_event_align:
    # Align the events to the reference genome.
    # The wildcard "chunk" can simply be "full", in cases
    # where there are no chunks
    input:
        sortedbam             = os.path.join( config["PATHOUT"], "{wcEvalign_sampleDir}", SUBDIR_SORTED_MINIMAPPED, "bam_chunks", "{wcEvalign_sampleName}_{wcEvalign_chunk}.sorted.bam"),
        NOTCALLED_indexedbam  = os.path.join( config["PATHOUT"], "{wcEvalign_sampleDir}", SUBDIR_SORTED_MINIMAPPED, "bam_chunks", "{wcEvalign_sampleName}_{wcEvalign_chunk}.sorted.bam.bai"),
        fastq_file            = os.path.join( config["PATHOUT"], "{wcEvalign_sampleDir}", SUBDIR_SYMLINKS,  "{wcEvalign_sampleName}_{wcEvalign_chunk}" +config["fastq_suffix"]),
        fastq_npi             = os.path.join( config["PATHOUT"], "{wcEvalign_sampleDir}", SUBDIR_SYMLINKS,  "{wcEvalign_sampleName}_{wcEvalign_chunk}" + config["fastq_suffix"] + ".index"),
        refgenome_fasta       = os.path.join( DIR_REFGENOME, config['ref']['Genome_version']+ ".fa" ),
        NOTCALLED_bwt         = os.path.join( DIR_REFGENOME, config['ref']['Genome_version']+ ".fa.bwt"),
        NOTCALLED_pac         = os.path.join( DIR_REFGENOME, config['ref']['Genome_version']+ ".fa.pac")
    output:
        Evaligned         = os.path.join( config["PATHOUT"], "{wcEvalign_sampleDir}", SUBDIR_EVENTALIGN, "tsv_chunks", 'Ealign_{wcEvalign_sampleName}_{wcEvalign_chunk}.tsv' )
    log:
        logfile  = os.path.join( config["PATHOUT"], "{wcEvalign_sampleDir}", SUBDIR_EVENTALIGN, "tsv_chunks", 'Ealign_{wcEvalign_sampleName}_{wcEvalign_chunk}.log')
    message: """---- Align events from sample {wildcards.wcEvalign_sampleName}, chunk {wildcards.wcEvalign_chunk} to the genome ----"""
    shell:
        " {nanopolish} eventalign --reads {input.fastq_file} --bam {input.sortedbam} --genome {input.refgenome_fasta} --scale-events  > {output}  2> {log.logfile} "

#------------------------------------------------------
# rule quickcheck: (TODO)
#------------------------------------------------------

rule index_sortedbam:
    # Index the sorted bam file with samtools
    input:
        sortedbam  = os.path.join( config["PATHOUT"], "{wcindexbam_sampleDir}", SUBDIR_SORTED_MINIMAPPED, "bam_chunks", "{wcindexbam_sampleName}_{wcindexbam_chunk}.sorted.bam")
    output:
        indexedbam = os.path.join( config["PATHOUT"], "{wcindexbam_sampleDir}", SUBDIR_SORTED_MINIMAPPED, "bam_chunks", "{wcindexbam_sampleName}_{wcindexbam_chunk}.sorted.bam.bai")
    log:
        logfile    = os.path.join( config["PATHOUT"], "{wcindexbam_sampleDir}", SUBDIR_SORTED_MINIMAPPED, "bam_chunks", '{wcindexbam_sampleName}_{wcindexbam_chunk}_samtoolsindex.log')
    message: """---- index the bam files for {wildcards.wcindexbam_sampleName} chunk {wildcards.wcindexbam_chunk} ----"""
    shell:
        " {SAMTOOLS} index  {input.sortedbam}  2> {log.logfile} "


# -----------------------------------------------------

rule create_readcurrent_GRL_obj:
    # produce GRangesList object of current data in .RDS format
    # from the .tsv files produced by eventalign
    # each entry of the list is a read:
    input:
        tsvfile      = lambda wc: get_chunkfiles( wc.wcreadGRL_samplename, os.path.join( config["PATHOUT"], wc.wcreadGRL_sampleDir, SUBDIR_EVENTALIGN, "tsv_chunks" ), "Ealign_", ".tsv", False )
    output:
        GRLreads     = os.path.join( config["PATHOUT"], "{wcreadGRL_sampleDir}", SUBDIR_GR, "{wcreadGRL_samplename}_reads_GRL.rds"),
        poremodel    = os.path.join( config["PATHOUT"], "{wcreadGRL_sampleDir}", SUBDIR_GR, "{wcreadGRL_samplename}_poremodel.tsv")
    params:
        Rfuncs_tsv2GRconv = Rfuncs_tsv2GRL,
        output_reads_GRL  = os.path.join( config["PATHOUT"], "{wcreadGRL_sampleDir}", SUBDIR_GR, "{wcreadGRL_samplename}_reads_GRL.rds"),
        output_poremodel  = os.path.join( config["PATHOUT"], "{wcreadGRL_sampleDir}", SUBDIR_GR, "{wcreadGRL_samplename}_poremodel.tsv"),
        Flatten           = config["execution"]["FlattenReads"],
        Ealign_files      = lambda wc: get_chunkfiles( wc.wcreadGRL_samplename, os.path.join( config["PATHOUT"], wc.wcreadGRL_sampleDir, SUBDIR_EVENTALIGN, "tsv_chunks" ), "Ealign_", ".tsv", True ),
        samplename          = "{wcreadGRL_samplename}"
    log:
        os.path.join( config["PATHOUT"], "{wcreadGRL_sampleDir}", SUBDIR_GR, "{wcreadGRL_samplename}_reads_GRL_conversion.log")
    message:
        fmt("Convert aligned NP reads from {input} to GRangesList object")
    shell:
        nice('Rscript', [ Rmain_tsv2GRL,
                         "--Rfuncs_tsv2GRconv={params.Rfuncs_tsv2GRconv}",
                         "--output_reads_GRL={params.output_reads_GRL}",
                         "--output_poremodel={params.output_poremodel}",
                         "--samplename={params.samplename}",
                         "--Flatten_reads={params.Flatten}",
                         "--logFile={log}",
                         "--Ealign_files={params.Ealign_files}"] )

# -----------------------------------------------------

rule np_index:
    # Index the reads and the fast5 files for nanopolish
    input:
        fast5_folder = os.path.join( config["PATHIN"], "{wcnpindex_sampleDir}/fast5/pass", "{wcnpindex_chunk}"),
#        fast5_folder = lambda wc: getPathCase( os.path.join( config["PATHIN"], "{wcnpindex_sampleDir}"), 'fast5', 'pass', wc.wcnpindex_chunk, input_data_type ),
        fastq_file   = os.path.join( config["PATHOUT"], "{wcnpindex_sampleDir}", SUBDIR_SYMLINKS, "{wcnpindex_samplename}_{wcnpindex_chunk}" + config["fastq_suffix"])
#        fastq_file   = lambda wc: os.path.join( config["PATHOUT"], "{wcnpindex_sampleDir}", SUBDIR_SYMLINKS, wc.wcnpindex_samplename + "_{wcnpindex_chunk}." + config["fastq_suffix"] )
    output:
        npi    = os.path.join( config["PATHOUT"], "{wcnpindex_sampleDir}", SUBDIR_SYMLINKS, "{wcnpindex_samplename}_{wcnpindex_chunk}"+config["fastq_suffix"]+".index"  ),
        fai    = os.path.join( config["PATHOUT"], "{wcnpindex_sampleDir}", SUBDIR_SYMLINKS, "{wcnpindex_samplename}_{wcnpindex_chunk}"+config["fastq_suffix"]+".index.fai" ),
        gzi    = os.path.join( config["PATHOUT"], "{wcnpindex_sampleDir}", SUBDIR_SYMLINKS, "{wcnpindex_samplename}_{wcnpindex_chunk}"+config["fastq_suffix"]+".index.gzi" ),
        readdb = os.path.join( config["PATHOUT"], "{wcnpindex_sampleDir}", SUBDIR_SYMLINKS, "{wcnpindex_samplename}_{wcnpindex_chunk}"+config["fastq_suffix"]+".index.readdb")
    params:
        options    = " index -d "
    log:
        logfile  = os.path.join( config["PATHOUT"], "{wcnpindex_sampleDir}", SUBDIR_SYMLINKS,  "{wcnpindex_samplename}_{wcnpindex_chunk}_npi.log" )
    message:
        fmt("Index the reads from chunk {wildcards.wcnpindex_chunk} against the fast5 files from the same.")
    shell:
        " nice -19 {nanopolish} {params.options} {input.fast5_folder} {input.fastq_file} 2> {log.logfile} "

# -----------------------------------------------------

rule merge_bam_files:
    # combine the ~4000 reads from each minimap2
    # alignment into a single bam file
    input:
        bami      = lambda wc: get_chunkfiles( wc.wcmerge_samplename, os.path.join( config["PATHOUT"], wc.wcmerge_sampleDir, SUBDIR_SORTED_MINIMAPPED, "bam_chunks"), "" , ".sorted.bam", False )
    output:
        sortedbam = os.path.join( config["PATHOUT"], "{wcmerge_sampleDir}", SUBDIR_SORTED_MINIMAPPED, "{wcmerge_samplename}.sorted.bam")
    log:
        logfile   = os.path.join( config["PATHOUT"], "{wcmerge_sampleDir}", SUBDIR_SORTED_MINIMAPPED, "{wcmerge_samplename}.mergingbam.log")
    message:
        fmt("Combining bam files from post-mapping fastq data.")
    shell:
        '{SAMTOOLS} merge {output} {input} '

#------------------------------------------------------
rule convert_sort_minimap:
    # convert minimap2 alignments from sam to bam format
    # and sort by position
    input:
        aligned     = os.path.join( config["PATHOUT"], "{wcsort_sampleDir}", SUBDIR_FILTERED_MINIMAP, "{wcsort_samplename}_{wcsort_chunk}.0filtered.sam")
    output:
        sortedbam   = os.path.join( config["PATHOUT"], "{wcsort_sampleDir}", SUBDIR_SORTED_MINIMAPPED, "bam_chunks", "{wcsort_samplename}_{wcsort_chunk}.sorted.bam")
#    params:
#
    log:
        logfile = os.path.join( config["PATHOUT"], "{wcsort_sampleDir}", SUBDIR_SORTED_MINIMAPPED, "bam_chunks", "{wcsort_samplename}_{wcsort_chunk}.sortingbam.log")
    message:
        fmt("Converting, sorting, and indexing input file: {input}")
    shell:
        'samtools view  -Sb  {input} | samtools sort > {output} && samtools index {output} 2> {log.logfile}'


#------------------------------------------------------
rule filter_nonaligned_minimap:
    # Check for alignment filter in minimap2s sam file:
    #  if != 4 then remove this read
    input:
        aligned  = os.path.join( config["PATHOUT"], "{wcfilter_sampleDir}", SUBDIR_ALIGNED_MINIMAP, "{wcfilter_samplename}_{wcfilter_chunk}.sam" )
    output:
        filtered = os.path.join( config["PATHOUT"], "{wcfilter_sampleDir}", SUBDIR_FILTERED_MINIMAP, "{wcfilter_samplename}_{wcfilter_chunk}.0filtered.sam" )
    log:
        log      = os.path.join( config["PATHOUT"], "{wcfilter_sampleDir}", SUBDIR_FILTERED_MINIMAP, "{wcfilter_samplename}_{wcfilter_chunk}.0filtering.log" )
    message:
        fmt("Filtering unaligned reads from alignment data")
    shell:
        " cat {input} | perl -lane 'print if $F[1] ne 4'  >  {output}   2> {log}"

#------------------------------------------------------
rule align_minimap:
    # use minimap2 to align the fastq reads to the reference
    # genome
    input:
        mmiref   = os.path.join( DIR_REFGENOME , config['ref']['Genome_version']+ ".mmi" ),
        fqlink   = os.path.join( config["PATHOUT"], "{wcalign_sampleDir}", SUBDIR_SYMLINKS, "{wcalign_samplename}_{wcalign_chunk}" + config["fastq_suffix"] )
    output:
        aligned  = os.path.join( config["PATHOUT"], "{wcalign_sampleDir}", SUBDIR_ALIGNED_MINIMAP,  "{wcalign_samplename}_{wcalign_chunk}.sam" )
    params:
        options  = config["execution"]["MM2_align_option"]
    log:
        log      = os.path.join( config["PATHOUT"], "{wcalign_sampleDir}", SUBDIR_ALIGNED_MINIMAP,  "{wcalign_samplename}_{wcalign_chunk}_alignment.log" )
    message:
        fmt("Aligning fastq reads from {input.fqlink} to {input.mmiref}; outputting to {output}")
    shell:
        "{MM2} {params} {input.mmiref} {input.fqlink} > {output}  2> {log}"

#------------------------------------------------------
rule create_fqsymlink:
    # create symlink to data:
    input:
        chunk_src     = lambda wc: os.path.join( config["PATHIN"], wc.wcfqlink_sampleDir, "fastq", "pass", config["samplelist"][wc.wcfqlink_samplename]["fastq_prefix"] + wc.wcfqlink_chunk + config["fastq_suffix"] )
    output:
        chunk_linkloc = os.path.join( config["PATHOUT"], "{wcfqlink_sampleDir}", SUBDIR_SYMLINKS, "{wcfqlink_samplename}_{wcfqlink_chunk}" + config["fastq_suffix"] )
    params:
        options  = " -s "
    message:
        fmt("Creating symbolic link to chunk data")
    shell:
        " ln {params} {input} {output}"


#------------------------------------------------------
# rule basecall:
#     # If we have fast5 data, but no fastq, then perform
#     # basecalling ourselves:
#     # docs: https://medium.com/@kepler_00/nanopore-gpu-basecalling-using-guppy-on-ubuntu-18-04-and-nvidia-docker-v2-with-a-rtx-2080-d875945e5c8d
#     # https://community.nanoporetech.com/protocols/Guppy-protocol-preRev/v/gpb_2003_v1_revj_14dec2018
#     input:
#         fast5_folder =
#     output:
#         fastq_folder =
#     params:
#         ???
#     log:
#         os.path.join( DIR_???, ??? )
#      message:
#         fmt("Performing base-calling of fast5 data from {input}")
#     shell:
#         guppy_basecaller -i <fast5_dir> -o <output_folder> -c dna_r9.4.1_450bps -x "cuda:0"
#
#------------------------------------------------------

