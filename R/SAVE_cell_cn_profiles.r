#====================================SAVE THE SAMPLED CELLS' CN PROFILES
SAVE_cell_cn_profiles <- function(package_simulation,
                                    filename=''){
    package_sample                <- package_simulation[[2]]
    sample_genotype_profiles      <- package_sample[[1]]

    filename    <- paste(filename,'.csv',sep='')
    write.csv(sample_genotype_profiles,filename,row.names=FALSE)
}
