#================================================SIMULATE FOCAL DELETION
SIMULATOR_FULL_PHASE_1_CN_focal_deletion <- function(genotype_to_react,genotype_daughter) {
#------------------------------------Find the new CN and driver profiles
#   Initialize the daughter's CN and driver profiles
    ploidy_chrom        <- genotype_list_ploidy_chrom[[genotype_daughter]]
    ploidy_allele       <- genotype_list_ploidy_allele[[genotype_daughter]]
    ploidy_block        <- genotype_list_ploidy_block[[genotype_daughter]]
    driver_count        <- genotype_list_driver_count[genotype_daughter]
    driver_map          <- genotype_list_driver_map[[genotype_daughter]]
#   Find information about the focal deletion
    while (1) {
#       Choose the chromosome to be focally deleted
        chrom           <- sample.int(N_chromosomes,size=1)
        no_strands      = ploidy_chrom[chrom]
        if (no_strands<=0) {
            next
        }
#       Choose the strand to be focally deleted
        strand          <- sample.int(no_strands,size=1)
#       Find the chromosome's centromere location and length
        centromere      <- vec_centromere_location[chrom]
        chrom_length    <- vec_CN_block_no[chrom]
#       Choose the chromosome arm to be focally deleted
        chrom_arm       <- sample.int(2,size=1)
        if (chrom_arm==1) {
            max_length  <- centromere
            }
        else { if (chrom_arm==2) {
            max_length  <- chrom_length-centromere
            }
        }
#       Choose the length of the focal deletion
        focal_length    <- max_length+1;
        while (focal_length>max_length) {
            focal_length <- 1+rgeom(n=1,prob=prob_CN_focal_deletion_length)
        }
#       Choose the region to be focally deleted
        block_start     <- (chrom_arm-1)*centromere + sample.int(max_length-focal_length+1,size=1)
        block_end       <- block_start+focal_length-1
        break
    }
#   Find all drivers located on this region
    if ((driver_count==0)||(length(intersect(which(driver_map[,4]>=block_start),which(driver_map[,4]<=block_end)))==0)) {
        pos_drivers_to_delete   <- c()
        }
    else {
        pos_drivers_to_delete   <- intersect(intersect(which(driver_map[,2]==chrom),which(driver_map[,3]==strand)),intersect(which(driver_map[,4]>=block_start),which(driver_map[,4]<=block_end)))
    }
    N_drivers_to_delete         <- length(pos_drivers_to_delete)
#   Update the chromosome strand allele identity
    ploidy_allele[[chrom]][[strand]][,block_start:block_end]    <- 0
#   Change the local CN on the deleted region
    ploidy_block[[chrom]][[strand]][block_start:block_end]      <- 0
#   Change the driver count
    driver_count                <- driver_count-N_drivers_to_delete
#   Delete the drivers
    if (N_drivers_to_delete>0) {
        driver_map              <- driver_map[-pos_drivers_to_delete,]
    }
#------------------------------------------------Output the new genotype
    genotype_list_ploidy_chrom[[genotype_daughter]]           <<- ploidy_chrom
    genotype_list_ploidy_allele[[genotype_daughter]]          <<- ploidy_allele
    genotype_list_ploidy_block[[genotype_daughter]]           <<- ploidy_block
    genotype_list_driver_count[genotype_daughter]             <<- driver_count
    genotype_list_driver_map[[genotype_daughter]]             <<- driver_map
    end                                                       <- length(evolution_genotype_changes[[genotype_daughter]])
    evolution_genotype_changes[[genotype_daughter]][[end+1]]  <<- list('focal-deletion',chrom,strand,block_start,block_end)
}
