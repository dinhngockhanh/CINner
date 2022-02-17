#=================================SIMULATE CHROMOSOME ARM MISSEGREGATION
SIMULATOR_FULL_PHASE_1_CN_chrom_arm_missegregation <- function(genotype_to_react,genotype_daughter_1,genotype_daughter_2) {
#------------------------------------Find the new CN and driver profiles
#   Find the daughter cells' current CN and driver profiles
    ploidy_chrom_1      <- genotype_list_ploidy_chrom[[genotype_daughter_1]]
    ploidy_allele_1     <- genotype_list_ploidy_allele[[genotype_daughter_1]]
    ploidy_block_1      <- genotype_list_ploidy_block[[genotype_daughter_1]]
    driver_count_1      <- genotype_list_driver_count[genotype_daughter_1]
    driver_map_1        <- genotype_list_driver_map[[genotype_daughter_1]]

    ploidy_chrom_2      <- genotype_list_ploidy_chrom[[genotype_daughter_2]]
    ploidy_allele_2     <- genotype_list_ploidy_allele[[genotype_daughter_2]]
    ploidy_block_2      <- genotype_list_ploidy_block[[genotype_daughter_2]]
    driver_count_2      <- genotype_list_driver_count[genotype_daughter_2]
    driver_map_2        <- genotype_list_driver_map[[genotype_daughter_2]]
#   Find information about the chromosome-arm missegregation
    while (1) {
#       Choose which cell to gain/lose the chromosome arm
        i_gain          <- sample.int(2,size=1)
#       Choose the chromosome to be mis-segregated
        chrom           <- sample.int(N_chromosomes,size=1)
        if (i_gain==1) {
            no_strands  <- ploidy_chrom_2[chrom]
            }
        else { if (i_gain==2) {
            no_strands  <- ploidy_chrom_1[chrom]
            }
        }
        if (no_strands<=0) {
            next
        }
#       Choose the strand to be mis-segregated
        strand          <- sample.int(no_strands,size=1)
#       Choose the chromosome arm to be mis-segregated
        centromere      <- vec_centromere_location[chrom]
        chrom_length    <- vec_CN_block_no[chrom]
        chrom_arm       <- sample.int(2,size=1)
#       Find the region to be mis-segregated
        if (chrom_arm==1) {
            block_start <- 1
            block_end   <- centromere
            }
        else { if (chrom_arm==2) {
            block_start <- centromere+1
            block_end   <- chrom_length
            }
        }
        break
    }
#   Find all drivers located on this strand arm in the losing cell
    if (i_gain==1) {
        driver_count_lose   <- driver_count_2
        driver_map_lose     <- driver_map_2
        }
    else { if (i_gain==2) {
        driver_count_lose   <- driver_count_1
        driver_map_lose     <- driver_map_1
        }
    }
#    if ((driver_count_lose==0)||(isempty(intersect(find(driver_map_lose(:,4)>=block_start),find(driver_map_lose(:,4)<=block_end))))) {
# check if A is empty by length(A) == 0
    if ((driver_count_lose==0)||(length(intersect(which(driver_map_lose[,4]>=block_start),which(driver_map_lose[,4]<=block_end)))==0)) {
        pos_drivers_to_move <- c()
        }
    else {
        pos_drivers_to_move <- intersect(intersect(which(driver_map_lose[,2]==chrom),which(driver_map_lose[,3]==strand)),intersect(which(driver_map_lose[,4]>=block_start),which(driver_map_lose[,4]<=block_end)))
    }
    N_drivers_to_move       <- length(pos_drivers_to_move)
#   Change the chromosome ploidy of daughter cell gaining the arm
    if (i_gain==1) {
        ploidy_chrom_1[chrom]   <- ploidy_chrom_1[chrom]+1
        }
    else { if (i_gain==2) {
        ploidy_chrom_2[chrom]   <- ploidy_chrom_2[chrom]+1
        }
    }
#   Update the chromosome strand allele identities of daughter cells
    if (i_gain==1) {
        chrom_ploidy                                                    <- ploidy_chrom_1[chrom]
        ploidy_allele_1[[chrom]][[chrom_ploidy]]                        <- matrix(0,nrow(ploidy_allele_2[[chrom]][[strand]]),ncol(ploidy_allele_2[[chrom]][[strand]]))
        ploidy_allele_1[[chrom]][[chrom_ploidy]][,block_start:block_end]<- ploidy_allele_2[[chrom]][[strand]][,block_start:block_end]
        ploidy_allele_2[[chrom]][[strand]][,block_start:block_end]      <- 0
        }
    else { if (i_gain==2) {
        chrom_ploidy                                                    <- ploidy_chrom_2[chrom]
        ploidy_allele_2[[chrom]][[chrom_ploidy]]                        <- matrix(0,nrow(ploidy_allele_1[[chrom]][[strand]]),ncol(ploidy_allele_1[[chrom]][[strand]]))
        ploidy_allele_2[[chrom]][[chrom_ploidy]][,block_start:block_end]<- ploidy_allele_1[[chrom]][[strand]][,block_start:block_end]
        ploidy_allele_1[[chrom]][[strand]][,block_start:block_end]      <- 0
        }
    }
#   Move the chromosome arm from losing cell to winning cell
    if (i_gain==1) {
        chrom_ploidy                                                    <- ploidy_chrom_1[chrom]
        ploidy_block_1[[chrom]][[chrom_ploidy]]                         <- rep(0,length(ploidy_block_2[[chrom]][[strand]]))
        ploidy_block_1[[chrom]][[chrom_ploidy]][block_start:block_end]  <- ploidy_block_2[[chrom]][[strand]][block_start:block_end]
        ploidy_block_2[[chrom]][[strand]][block_start:block_end]        <- 0
        }
    else { if (i_gain==2) {
        chrom_ploidy                                                    <- ploidy_chrom_2[chrom]
        ploidy_block_2[[chrom]][[chrom_ploidy]]                         <- rep(0,length(ploidy_block_1[[chrom]][[strand]]))
        ploidy_block_2[[chrom]][[chrom_ploidy]][block_start:block_end]  <- ploidy_block_1[[chrom]][[strand]][block_start:block_end]
        ploidy_block_1[[chrom]][[strand]][block_start:block_end]        <- 0
        }
    }
#   Move the drivers from losing cell to winning cell
    if ((i_gain==1)&&(N_drivers_to_move>0)){
        driver_map_new_1                                    <- driver_map_2[pos_drivers_to_move,]
        if (!is.matrix(driver_map_new_1)){
            driver_map_new_1                                <- matrix(driver_map_new_1,nrow=1)
        }
        driver_map_new_1[,3]                                <- ploidy_chrom_1[chrom]
        driver_map_1                                        <- rbind(driver_map_1,driver_map_new_1)
        driver_map_2                                        <- driver_map_2[-pos_drivers_to_move,]
        if (!is.matrix(driver_map_2)){
            driver_map_2                                    <- matrix(driver_map_2,nrow=1)
        }
    }
    else{if ((i_gain==2)&&(N_drivers_to_move>0)){
        driver_map_new_2                                    <- driver_map_1[pos_drivers_to_move,]
        if (!is.matrix(driver_map_new_2)){
            driver_map_new_2                                <- matrix(driver_map_new_2,nrow=1)
        }
        driver_map_new_2[,3]                                <- ploidy_chrom_2[chrom]
        driver_map_2                                        <- rbind(driver_map_2,driver_map_new_2)
        driver_map_1                                        <- driver_map_1[-pos_drivers_to_move,]
        if (!is.matrix(driver_map_1)){
            driver_map_1                                    <- matrix(driver_map_1,nrow=1)
        }
    }}
#   Change the driver count in each daughter cell
    driver_unique_1                     <- unique(driver_map_1[,1])
    driver_unique_1                     <- driver_unique_1[driver_unique_1!=0]
    driver_count_1                      <- length(driver_unique_1)

    driver_unique_2                     <- unique(driver_map_2[,1])
    driver_unique_2                     <- driver_unique_2[driver_unique_2!=0]
    driver_count_2                      <- length(driver_unique_2)
#-----------------------------------------------Output the new genotypes
    genotype_list_ploidy_chrom[[genotype_daughter_1]]       <<- ploidy_chrom_1
    genotype_list_ploidy_allele[[genotype_daughter_1]]      <<- ploidy_allele_1
    genotype_list_ploidy_block[[genotype_daughter_1]]       <<- ploidy_block_1
    genotype_list_driver_count[genotype_daughter_1]         <<- driver_count_1
    genotype_list_driver_map[[genotype_daughter_1]]         <<- driver_map_1
    if (i_gain==1) {
        loc_end                                                         <- length(evolution_genotype_changes[[genotype_daughter_1]])
        evolution_genotype_changes[[genotype_daughter_1]][[loc_end+1]]  <<- c('chromosome-arm-missegregation',chrom,strand,chrom_arm,1)
        }
    else { if (i_gain==2) {
        loc_end                                                         <- length(evolution_genotype_changes[[genotype_daughter_1]])
        evolution_genotype_changes[[genotype_daughter_1]][[loc_end+1]]  <<- c('chromosome-arm-missegregation',chrom,strand,chrom_arm,-1)
        }
    }

    genotype_list_ploidy_chrom[[genotype_daughter_2]]       <<- ploidy_chrom_2
    genotype_list_ploidy_allele[[genotype_daughter_2]]      <<- ploidy_allele_2
    genotype_list_ploidy_block[[genotype_daughter_2]]       <<- ploidy_block_2
    genotype_list_driver_count[genotype_daughter_2]         <<- driver_count_2
    genotype_list_driver_map[[genotype_daughter_2]]         <<- driver_map_2
    if (i_gain==1) {
        loc_end                                                         <- length(evolution_genotype_changes[[genotype_daughter_2]])
        evolution_genotype_changes[[genotype_daughter_2]][[loc_end+1]]  <<- c('chromosome-arm-missegregation',chrom,strand,chrom_arm,-1)
        }
    else { if (i_gain==2) {
        loc_end                                                         <- length(evolution_genotype_changes[[genotype_daughter_2]])
        evolution_genotype_changes[[genotype_daughter_2]][[loc_end+1]]  <<- c('chromosome-arm-missegregation',chrom,strand,chrom_arm,1)        
        }
    }
}
