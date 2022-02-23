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


print(genotype_list_driver_count)


}
