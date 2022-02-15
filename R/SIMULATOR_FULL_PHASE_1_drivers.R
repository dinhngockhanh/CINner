#=================================================SIMULATE DRIVER EVENTS
SIMULATOR_FULL_PHASE_1_drivers <- function(genotype_to_react,genotype_daughter_1,genotype_daughter_2) {
#------------------------------------Find the new CN and driver profiles
#   Find the mother cell's CN and driver profiles
    ploidy_chrom            <- genotype_list_ploidy_chrom[[genotype_to_react]]
    ploidy_allele           <- genotype_list_ploidy_allele[[genotype_to_react]]
    ploidy_block            <- genotype_list_ploidy_block[[genotype_to_react]]
    driver_count            <- genotype_list_driver_count[[genotype_to_react]]
    driver_map              <- genotype_list_driver_map[[genotype_to_react]]
#   Find the daughter cells' current CN and driver profiles
    ploidy_chrom_1          <- genotype_list_ploidy_chrom[[genotype_daughter_1]]
    ploidy_allele_1         <- genotype_list_ploidy_allele[[genotype_daughter_1]]
    ploidy_block_1          <- genotype_list_ploidy_block[[genotype_daughter_1]]
    driver_count_1          <- genotype_list_driver_count[[genotype_daughter_1]]
    driver_map_1            <- genotype_list_driver_map[[genotype_daughter_1]]

    ploidy_chrom_2          <- genotype_list_ploidy_chrom[[genotype_daughter_2]]
    ploidy_allele_2         <- genotype_list_ploidy_allele[[genotype_daughter_2]]
    ploidy_block_2          <- genotype_list_ploidy_block[[genotype_daughter_2]]
    driver_count_2          <- genotype_list_driver_count[[genotype_daughter_2]]
    driver_map_2            <- genotype_list_driver_map[[genotype_daughter_2]]
#---------------------------Find the eligible driver genes to be mutated
#   Eliminate already mutated genes from list of eligible genes to mutate
    driver_library_eligible                             <- driver_library
    if (driver_count>0){
        vec_loc                                         <- driver_map[,1]
        driver_library_eligible                         <- driver_library_eligible[-vec_loc,]
    }
#   If no more genes to mutate then no new drivers

print(driver_library_eligible)

    if (all(is.na(driver_library_eligible))){
        return()
    }
#---------------------------Find copy count of each eligible driver gene
#   Find copy number of each eligible driver gene
    for (driver in 1:nrow(driver_library_eligible)){
        chrom                                           <- driver_library_eligible$Chromosome[driver]
        block                                           <- driver_library_eligible$Bin[driver]
        no_strands                                      <- ploidy_chrom[chrom]
        driver_copy                                     <- 0
        for (strand in 1:no_strands){
            driver_copy                                 <- driver_copy+ploidy_block[[chrom]][[strand]][block]
        }
        driver_library_eligible$Copy_count[driver]      <- driver_copy
    }
#   Delete driver genes with zero copy
    vec_delete                                          <- which(driver_library_eligible$Copy_count==0)
    driver_library_eligible                             <- driver_library_eligible[-vec_delete,]
#   If no more genes to mutate then no new drivers
    if (all(is.na(driver_library_eligible))){
        return()
    }
#----------------Place the new driver in the daughter cells' driver maps
#---Choose one driver from the eligible list to mutate
    driver                                              <- sample.int(nrow(driver_library_eligible),size=1)
    driver_library_eligible                             <- driver_library_eligible[driver,]
#---Find the new driver's address
#   Find the new driver's ID
    driver_ID                                           <- which(driver_library$Gene_ID==driver_library_eligible$Gene_ID)
#   Find the new driver's chromosome
    chrom                                               <- driver_library_eligible$Chromosome
    while (1){
#       Find the new driver's strand
        no_strands                                      <- ploidy_chrom[chrom]
        strand                                          <- sample.int(no_strands,size=1)
#       Find the new driver's block
        block                                           <- driver_library_eligible$Bin
        no_units                                        <- ploidy_block[[chrom]][[strand]][block]
        if (no_units<=0){
            next
        }
#       Find the new driver's unit
        unit                                            <- sample.int(no_units,size=1)
        break
    }
#---Place the new driver in the daughter cells' driver maps
    driver_count_1                                      <- driver_count_1+1
    if (driver_count_1==1){
        driver_map_1                                    <- matrix(c(driver_ID,chrom,strand,block,unit),nrow=1)
    }
    else{
        driver_map_1                                    <- rbind(driver_map_1,c(driver_ID,chrom,strand,block,unit))
    }
    driver_count_2                                      <- driver_count_2+1
    if (driver_count_2==1){
        driver_map_2                                    <- matrix(c(driver_ID,chrom,strand,block,unit),nrow=1)
    }
    else{
        driver_map_2                                    <- rbind(driver_map_2,c(driver_ID,chrom,strand,block,unit))
    }

print(driver_count)
print(driver_map)
print('---------------')
print(driver_count)
print(driver_map)

#-----------------------------------------------Output the new genotypes
    genotype_list_driver_count[genotype_daughter_1]     <- driver_count_1
    genotype_list_driver_map[[genotype_daughter_1]]     <- driver_map_1
    evolution_genotype_changes[[genotype_daughter_1]]   <- c(evolution_genotype_changes[[genotype_daughter_1]],list('new-driver',driver_ID))

    genotype_list_driver_count(genotype_daughter_2)     <- driver_count_2
    genotype_list_driver_map[[genotype_daughter_2]]     <- driver_map_2
    evolution_genotype_changes[[genotype_daughter_2]]   <- c(evolution_genotype_changes[[genotype_daughter_2]],list('new-driver',driver_ID))
}
