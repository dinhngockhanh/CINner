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


print(all_clones_ID)
print(all_clones_population)


}
