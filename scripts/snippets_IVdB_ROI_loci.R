
ROIrds_template=readRDS("/fast/AG_Akalin/bosberg/nanopore/ref/Regions_of_interest/m6A_putlocs_Linder.rds")
IV_out = readRDS("/scratch/AG_Akalin/bosberg/nanopore/pipeline_output/20190702_1647_IVdBm6A_FAK57816_38717116/07_GRprocessing/m6a_0_10_100_read_ROIolap_barcodes.rds") 
ref_Genome_IVdb    <- readDNAStringSet( "/scratch/AG_Akalin/refGenomes/IV_DeBRuijn_barcoded/IV_dB_Bcoded.fa" )


ROI_BC_GR = GRanges(seqnames = paste0("BC",as.character(seq(1:7))),
                    ranges = IRanges( start=1538, end = 1549 ),
                    strand = "+"
                    )
IvdB_ROI_barcodes = list( 
                         "refdatName"  = "IVdb_barcodes",
                         "assembly"    = "IVdb" ,
                         "Region_groups" = list( "barcodes" = ROI_BC_GR)
    )
saveRDS( object = IvdB_ROI_barcodes,
         file = "/fast/AG_Akalin/bosberg/nanopore/ref/Regions_of_interest/IV_dB/IVdB_ROI_barcodes.rds" )


ROI_homopolymer_GR = GRanges(  seqnames = rep( paste0("BC",as.character(seq(1:7))), 4),
                               ranges   = IRanges( start = rep(c(512,767,1021,1277), each = 7), 
                                                   end   = rep(c(512,767,1021,1277), each = 7) ),
                               strand = "+"
                             )
IvdB_ROI_homoP_rds = list(
    "refdatName" = "IVdb_homoP",
    "assembly"   = "IVdb" ,
    "Region_groups" = list( "homoP" = ROI_homopolymer_GR)
    )
saveRDS( object = IvdB_ROI_homoP_rds,
         file = "/fast/AG_Akalin/bosberg/nanopore/ref/Regions_of_interest/IV_dB/IVdB_ROI_homoP.rds" )

         IvdB_ROI_homoP_rds =

  #=====================================
         
IV_dB_homoP_olaps   = readRDS("/scratch/AG_Akalin/bosberg/nanopore/pipeline_output/20190702_1647_IVdBm6A_FAK57816_38717116/07_GRprocessing/m6a_0_10_100_read_ROIolap_homoP.rds")
IV_dB_barcode_olaps = readRDS("/scratch/AG_Akalin/bosberg/nanopore/pipeline_output/20190702_1647_IVdBm6A_FAK57816_38717116/07_GRprocessing/m6a_0_10_100_read_ROIolap_barcodes.rds")
         
         