#===================SIMULATE TERMINAL COPY-NEUTRAL LOSS OF HETEROGENEITY
SIMULATOR_FULL_PHASE_1_CN_cnloh_terminal <- function(genotype_to_react,genotype_daughter) {
#------------------------------------Find the new CN and driver profiles
#   Initialize the daughter's CN and driver profiles
    ploidy_chrom        <- genotype_list_ploidy_chrom[[genotype_daughter]]
    ploidy_allele       <- genotype_list_ploidy_allele[[genotype_daughter]]
    ploidy_block        <- genotype_list_ploidy_block[[genotype_daughter]]
    driver_count        <- genotype_list_driver_count[genotype_daughter]
    driver_map          <- genotype_list_driver_map[[genotype_daughter]]
#   Find information about the interstitial CN-LOH
    while (1) {
#       Choose the chromosome to harbor the interstitial CN-LOH
        chrom           <- sample.int(N_chromosomes,size=1)
        no_strands      <- ploidy_chrom[chrom]
        if (no_strands<=1) {
            next
        }
#       Choose the strands to donate and receive the DNA
        vec_strands     <- sample.int(no_strands,size=2)
        strand_give     <- vec_strands[1]
        strand_take     <- vec_strands[2]
#       Find the chromosome's centromere location and length
        centromere      <- vec_centromere_location[chrom]
        chrom_length    <- vec_CN_block_no[chrom]
#       Choose the chromosome arm to harbor the terminal CN-LOH
        chrom_arm       <- sample.int(2,size=1)
        if (chrom_arm==1){
            max_length  <- centromere
            }
        else {if (chrom_arm==2){
            max_length  <- chrom_length-centromere
            }
        }
#       Choose the length of the terminal CN-LOH
        cnloh_length    <- max_length+1
        while (cnloh_length>max_length) {
            cnloh_length<- 1+rgeom(n=1,prob_CN_cnloh_terminal_length)
        }
#       Choose the region to harbor the terminal CN-LOH
        if (chrom_arm==1){
            block_start <- 1
            block_end   <- cnloh_length
        }
        else {if (chrom_arm==2){
            block_start <- chrom_length-cnloh_length+1
            block_end   <- chrom_length
        }}
        break
    }
#   Find all drivers to lose in the strand to receive the DNA
    if ((driver_count==0)||(length(intersect(which(driver_map[,4]>=block_start),which(driver_map[,4]<=block_end)))==0)) {
        pos_drivers_to_delete    <- c()
        }
    else {
        pos_drivers_to_delete    <- intersect(intersect(which(driver_map[,2]==chrom),which(driver_map[,3]==strand_take)),intersect(which(driver_map[,4]>=block_start),which(driver_map[,4]<=block_end)))
    }
    N_drivers_to_delete          <- length(pos_drivers_to_delete)
#   Find all drivers to gain in the strand to donate the DNA
    if ((driver_count==0)||(length(intersect(which(driver_map[,4]>=block_start),which(driver_map[,4]<=block_end)))==0)) {
        pos_drivers_to_gain     <- c()
        }
    else {
        pos_drivers_to_gain     <- intersect(intersect(which(driver_map[,2]==chrom),which(driver_map[,3]==strand_give)),intersect(which(driver_map[,4]>=block_start),which(driver_map[,4]<=block_end)))
    }
    N_drivers_to_gain           <- length(pos_drivers_to_gain)
#   Update the chromosome strand allele identity of the strand to receive the DNA
    ploidy_allele_strand_give   <- ploidy_allele[[chrom]][[strand_give]]
    if (nrow(ploidy_allele_strand_give)>nrow(ploidy_allele[[chrom]][[strand_take]])){
        ploidy_allele[[chrom]][[strand_take]]                                                                   <- rbind(ploidy_allele[[chrom]][[strand_take]],matrix(0,nrow(ploidy_allele_strand_give)-nrow(ploidy_allele[[chrom]][[strand_take]]),ncol(ploidy_allele_strand_give)))
    }
    ploidy_allele[[chrom]][[strand_take]][1:nrow(ploidy_allele[[chrom]][[strand_take]]),block_start:block_end]  <- 0;
    ploidy_allele[[chrom]][[strand_take]][1:nrow(ploidy_allele_strand_give),block_start:block_end]              <- ploidy_allele_strand_give[nrow(ploidy_allele_strand_give),block_start:block_end]
#   Change the local CN on the strand to receive the DNA
    ploidy_block[[chrom]][[strand_take]][block_start:block_end]                                                 <- ploidy_block[[chrom]][[strand_give]][block_start:block_end]
#   Copy the drivers gained during terminal CN-LOH
    if (N_drivers_to_gain>0){
        driver_map_new                      <- driver_map[pos_drivers_to_gain,]
        if (!is.matrix(driver_map_new)){
            driver_map_new                  <- matrix(driver_map_new,nrow=1)
        }
        driver_map_new[,3]                  <- strand_take
        driver_map                          <- rbind(driver_map,driver_map_new)
    }
#   Remove the drivers lost during terminal CN-LOH
    if (N_drivers_to_delete>0){
        driver_map                          <- driver_map[-pos_drivers_to_delete,]
        if (!is.matrix(driver_map)){
            driver_map                      <- matrix(driver_map,nrow=1)
        }
    }
#   Change the driver count
    driver_unique                           <- unique(driver_map[,1])
    driver_unique                           <- driver_unique[driver_unique!=0]
    driver_count                            <- length(driver_unique)
#------------------------------------------------Output the new genotype
    genotype_list_ploidy_chrom[[genotype_daughter]]                 <<- ploidy_chrom
    genotype_list_ploidy_allele[[genotype_daughter]]                <<- ploidy_allele
    genotype_list_ploidy_block[[genotype_daughter]]                 <<- ploidy_block
    genotype_list_driver_count[genotype_daughter]                   <<- driver_count
    genotype_list_driver_map[[genotype_daughter]]                   <<- driver_map
    loc_end                                                         <- length(evolution_genotype_changes[[genotype_daughter]])
    evolution_genotype_changes[[genotype_daughter]][[loc_end+1]]    <<- c('cnloh-terminal',chrom,strand_give,strand_take,block_start,block_end)
}
