#===============================PRODUCE MODEL VARIABLE FILES FROM SAMPLE
EXTRA_build_model_variables_from_sample <- function(model,package_output){
#---------------------------------Input the clonal populations in sample
    package_clonal_evolution            <- package_output[[1]]

    package_sample                      <- package_output[[2]]

    genotype_list_ploidy_chrom          <- package_clonal_evolution[[5]]
    genotype_list_ploidy_block          <- package_clonal_evolution[[6]]
    genotype_list_ploidy_allele         <- package_clonal_evolution[[7]]
    genotype_list_driver_count          <- package_clonal_evolution[[8]]
    genotype_list_driver_map            <- package_clonal_evolution[[9]]

    sample_clone_ID                     <- package_sample[[3]]
#----------------------Output CN-driver profiles as model variable files
#---Find all unique clones in sample
    all_clones_ID                       <- unique(sample_clone_ID)
    all_clones_population               <- rep(0,length(all_clones_ID))
    for (i in 1:length(all_clones_ID)){
        all_clones_population[i]        <- length(which(sample_clone_ID==all_clones_ID[i]))
    }
    N_all_clones                        <- length(all_clones_ID)
#---Create and output model variables - copy number state
    TABLE_CLONAL_COPY_NUMBER_PROFILES   <- data.frame()
    HEADER_CLONAL_COPY_NUMBER_PROFILES  <- c()
#   Initialize table of CN profiles with chromosome and bin positions
    vec_chrom                           <- c()
    vec_bin                             <- c()
    for (chrom in 1:N_chromosomes){
        bin_count                       <- vec_CN_block_no[chrom]
        vec_chrom                       <- c(vec_chrom,rep(chrom,bin_count))
        vec_bin                         <- c(vec_bin (1:bin_count))
    }
    TABLE_CLONAL_COPY_NUMBER_PROFILES   <- data.frame(vec_chrom,vec_bin)
    colnames(TABLE_CLONAL_COPY_NUMBER_PROFILES) <- c('Chromosome','Bin')





print(TABLE_CLONAL_COPY_NUMBER_PROFILES)


}
