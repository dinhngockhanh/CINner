#===========================================SIMULATE FOCAL AMPLIFICATION
SIMULATOR_FULL_PHASE_1_CN_focal_amplification <- function(genotype_to_react,genotype_daughter) {
#------------------------------------Find the new CN and driver profiles
#   Initialize the daughter's CN and driver profiles
    ploidy_chrom        <- genotype_list_ploidy_chrom[[genotype_daughter]]
    ploidy_allele       <- genotype_list_ploidy_allele[[genotype_daughter]]
    ploidy_block        <- genotype_list_ploidy_block[[genotype_daughter]]
    driver_count        <- genotype_list_driver_count[genotype_daughter]
    driver_map          <- genotype_list_driver_map[[genotype_daughter]]
#   Find information about the focal amplification
    while (1) {
#       Choose the chromosome to be focally amplified
        chrom           <- sample.int(N_chromosomes,size=1)
        no_strands      <- ploidy_chrom[chrom]
        if (no_strands<=0) {
            next
        }
#       Choose the strand to be focally amplified
        strand          <- sample.int(no_strands,size=1)
#       Find the chromosome's centromere location and length
        centromere      <- vec_centromere_location[chrom]
        chrom_length    <- vec_CN_block_no[chrom]
#       Choose the chromosome arm to be focally amplified
        chrom_arm       <- sample.int(2,size=1)
        if (chrom_arm==1) {
            max_length  <- centromere
            }
        else { if (chrom_arm==2) {
            max_length  <- chrom_length-centromere
            }
        }
#       Choose the length of the focal amplification
        focal_length    <- max_length+1
        while (focal_length>max_length) {
            focal_length<- 1+rgeom(n=1,prob_CN_focal_amplified_length)
        }
#       Choose the region to be focally amplified
        block_start     <- (chrom_arm-1)*centromere + sample.int(max_length-focal_length+1,size=1)
        block_end       <- block_start+focal_length-1
        break
    }
#   Find all drivers located on this region
    if ((driver_count==0)||(length(intersect(which(driver_map[,4]>=block_start),which(driver_map[,4]<=block_end)))==0)) {
        pos_drivers_to_amplify   <- c()
        }
    else {
        pos_drivers_to_amplify   <- intersect(intersect(which(driver_map[,2]==chrom),which(driver_map[,3]==strand)),intersect(which(driver_map[,4]>=block_start),which(driver_map[,4]<=block_end)))
    }
    N_drivers_to_amplify         <- length(pos_drivers_to_amplify)
#   Update the chromosome strand allele identity
    for(block in block_start:block_end){
        block_CN                <- ploidy_block[[chrom]][[strand]][block]
        if(2*block_CN > nrow(ploidy_allele[[chrom]][[strand]])){
            ploidy_allele[[chrom]][[strand]]                            <- rbind(ploidy_allele[[chrom]][[strand]],matrix(0,2*block_CN-nrow(ploidy_allele[[chrom]][[strand]]),ncol(ploidy_allele[[chrom]][[strand]])))
        }

        # ploidy_allele[[chrom]][[strand]][block_CN+1:2*block_CN,block]   <- ploidy_allele[[chrom]][[strand]][1:block_CN,block]

        for(i_unit in block_CN+1:2*block_CN){
print('>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>')
print(block_CN+1)
print(2*block_CN)
print(i_unit)
print(nrow(ploidy_allele[[chrom]][[strand]]))
print(block)
print(ncol(ploidy_allele[[chrom]][[strand]]))

            ploidy_allele[[chrom]][[strand]][i_unit,block]              <- ploidy_allele[[chrom]][[strand]][i_unit-block_CN,block]
        }
    }
#   Change the local CN on the amplified region
    ploidy_block[[chrom]][[strand]][block_start:block_end] <- 2*ploidy_block[[chrom]][[strand]][block_start:block_end]
#   Change the driver count
    driver_count                 <- driver_count+N_drivers_to_amplify
#   Amplify the drivers
    if (N_drivers_to_amplify>0) {
#        driver_map(end+1:end+N_drivers_to_amplify,:)        = driver_map(pos_drivers_to_amplify,:);
        driver_map               <- rbind(driver_map,driver_map[pos_drivers_to_amplify,])
        for (driver in (nrow(driver_map)+1-N_drivers_to_amplify):nrow(driver_map)) {
            block                <- driver_map[driver,4]
            ploidy_block_driver  <- ploidy_block[[chrom]][[strand]][block]/2
            driver_map[driver,5] <- driver_map[driver,5]+ploidy_block_driver
        }
    }
#------------------------------------------------Output the new genotype
    genotype_list_ploidy_chrom[[genotype_daughter]]           <<- ploidy_chrom
    genotype_list_ploidy_allele[[genotype_daughter]]          <<- ploidy_allele
    genotype_list_ploidy_block[[genotype_daughter]]           <<- ploidy_block
    genotype_list_driver_count[genotype_daughter]             <<- driver_count
    genotype_list_driver_map[[genotype_daughter]]             <<- driver_map
    evolution_genotype_changes[[genotype_daughter]]           <<- c(evolution_genotype_changes[[genotype_daughter]],list('focal-amplification',chrom,strand,block_start,block_end))
}
