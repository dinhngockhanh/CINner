#=================================================SIMULATE DRIVER EVENTS
SIMULATOR_FULL_PHASE_1_drivers <- function(genotype_to_react,genotype_daughter_1,genotype_daughter_2) {
#-----------------------------------------Find the number of new drivers
    DNA_length              <- genotype_list_DNA_length[[genotype_to_react]]
    no_new_drivers          <- 0
    while (no_new_drivers<=0) {
        no_new_drivers      <- 1
        # no_new_drivers  = poissrnd(rate_driver*DNA_length);
    }
#------------------------------------Place the new drivers on the genome
#   Find the mother cell's CN and driver profiles
    ploidy_chrom            <- genotype_list_ploidy_chrom[[genotype_to_react]]
    ploidy_block            <- genotype_list_ploidy_block[[genotype_to_react]]
    driver_count            <- genotype_list_driver_count[genotype_to_react]
    driver_map              <- genotype_list_driver_map[[genotype_to_react]]
    if (is.null(driver_map)) {
        vec_old_driver_ID   <- c()
        }
    else { if (nrow(driver_map)==1) {
        vec_old_driver_ID   <- driver_map[1]
        }
    else {
        vec_old_driver_ID   <- t(driver_map[,1])
    }}
    vec_new_driver_ID       <- c()
#   Find the daughter cells' current CN and driver profiles
    ploidy_chrom_1          <- genotype_list_ploidy_chrom[[genotype_daughter_1]]
    ploidy_block_1          <- genotype_list_ploidy_block[[genotype_daughter_1]]
    driver_count_1          <- genotype_list_driver_count[genotype_daughter_1]
    driver_map_1            <- genotype_list_driver_map[[genotype_daughter_1]]

    ploidy_chrom_2          <- genotype_list_ploidy_chrom[[genotype_daughter_2]]
    ploidy_block_2          <- genotype_list_ploidy_block[[genotype_daughter_2]]
    driver_count_2          <- genotype_list_driver_count[genotype_daughter_2]
    driver_map_2            <- genotype_list_driver_map[[genotype_daughter_2]]
#   Find total number of drivers in library
#    no_drivers_in_library   = size(driver_library,1); The first dimension of driver_library?
    no_drivers_in_library   <- dim(driver_library)[1]
#   Distribute the new driver mutations
    for (driver in 1:no_new_drivers) {
#       Choose a random gene and its location
        while (1==1) {
#           If cell has acquired all possible drivers, then no new driver mutations
            if (length(vec_old_driver_ID)>=no_drivers_in_library) {
                break
            }
#           Choose a random gene in the driver library
#            driver_ID               = randi(no_drivers_in_library);
            driver_ID                   <- sample.int(no_drivers_in_library,size=1)
#           Make sure the gene is not already mutated
            if (is.element(driver_ID,vec_old_driver_ID)) {
                next # It is continue of Matlab
            }
#           Find its chromosome
            chrom                       <- driver_library[driver_ID,3] # driver_library is a data.frame
            no_strands                  <- ploidy_chrom[chrom]
            if (no_strands<=0) {
                next
            }
#           Find its strand
            strand                      <- sample.int(no_strands,size=1)
#           Find its block
            block                       <- driver_library[driver_ID,4]
            no_units                    <- ploidy_block[[chrom]][[strand]][block]
            if (no_units<=0) {
                next
            }
#           Find its unit
            unit                        <- sample.int(no_units,size=1)
            break
        }
#       Place the driver on the site for daughter cells
#        vec_new_driver_ID(end+1)        = driver_ID; append one more element
        vec_new_driver_ID               <- c(vec_new_driver_ID, driver_ID)

        driver_count_1                  <- driver_count_1+1
        driver_map_1[driver_count_1,]   <- c(driver_ID, chrom, strand, block, unit)

        driver_count_2                  <- driver_count_2+1
        driver_map_2[driver_count_2,]   <- c(driver_ID, chrom, strand, block, unit)
    }
#-------------Delete new drivers if they don't follow the mutation order
    SIMULATOR_FULL_PHASE_1_driver_order(driver_count_1,driver_map_1)
    SIMULATOR_FULL_PHASE_1_driver_order(driver_count_2,driver_map_2)
#-----------------------------------------------Output the new genotypes
    genotype_list_driver_count[genotype_daughter_1]         <- driver_count_1
    genotype_list_driver_map[[genotype_daughter_1]]         <- driver_map_1
#    evolution_genotype_changes{genotype_daughter_1}{end+1}  = {'new-driver',vec_new_driver_ID};
    evolution_genotype_changes[[genotype_daughter_1]]         <- c(evolution_genotype_changes[[genotype_daughter_1]],list('new-driver',vec_new_driver_ID))

    genotype_list_driver_count(genotype_daughter_2)         <- driver_count_2
    genotype_list_driver_map[[genotype_daughter_2]]           <- driver_map_2
    evolution_genotype_changes[[genotype_daughter_2]]         <- c(evolution_genotype_changes[[genotype_daughter_2]],list('new-driver',vec_new_driver_ID))
}
