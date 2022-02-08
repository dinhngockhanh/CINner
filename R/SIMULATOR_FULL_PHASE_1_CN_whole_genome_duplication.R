#======================================SIMULATE WHOLE GENOME DUPLICATION
SIMULATOR_FULL_PHASE_1_CN_whole_genome_duplication <- function(genotype_to_react,genotype_daughter_1,genotype_daughter_2) {
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
#   Change the chromosome ploidy of daughter cells
    ploidy_chrom_1      <- 2*ploidy_chrom_1
    ploidy_chrom_2      <- 2*ploidy_chrom_2

#   Update the chromosome strand allele identities of daughter cells
    for (chrom in 1:N_chromosomes) {
        chrom_ploidy    <- ploidy_chrom_1[chrom]
        for (strand in 1:chrom_ploidy/2) {

print('+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++')
print(chrom_ploidy/2+strand)
print(ploidy_allele_1[[chrom]][[chrom_ploidy/2+strand]])
print('+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++')
print(ploidy_allele_1[[chrom]][[strand]])

            ploidy_allele_1[[chrom]][[chrom_ploidy/2+strand]]   <- ploidy_allele_1[[chrom]][[strand]]
        }
        chrom_ploidy    <- ploidy_chrom_2[chrom]
        for (strand in 1:chrom_ploidy/2) {
            ploidy_allele_2[[chrom]][[chrom_ploidy/2+strand]]   <- ploidy_allele_2[[chrom]][[strand]]
        }
    }

#   Multiply the chromosome strands in each daughter cell
    for (chrom in 1:N_chromosomes) {
        chrom_ploidy    <- ploidy_chrom_1[chrom]
        for (strand in 1:chrom_ploidy/2) {
            ploidy_block_1[[chrom]][[chrom_ploidy/2+strand]]    <- ploidy_block_1[[chrom]][[strand]]
        }
        chrom_ploidy    <- ploidy_chrom_2[chrom]
        for (strand in 1:chrom_ploidy/2) {
            ploidy_block_2[[chrom]][[chrom_ploidy/2+strand]]    <- ploidy_block_2[[chrom]][[strand]]
        }
    }
#   Change the driver count in each daughter cell
    driver_count_1      <- 2*driver_count_1
    driver_count_2      <- 2*driver_count_2
#   Multiply the drivers in each daughter cell
    if (driver_count_1>0){
        driver_map_1[driver_count_1/2+1:driver_count_1,]   <- driver_map_1[1:driver_count_1/2,]
        for (driver in length(driver_map_1)/2+1:length(driver_map_1)) {
            chrom                           <- driver_map_1[driver,2]
            chrom_ploidy                    <- ploidy_chrom_1[chrom]
            driver_map_1[driver,3]          <- driver_map_1[driver,3]+chrom_ploidy
        }
    }
    if (driver_count_2>0){
        driver_map_2[driver_count_2/2+1:driver_count_2,]   <- driver_map_2[1:driver_count_2/2,]
        for (driver in length(driver_map_2)/2+1:length(driver_map_2)) {
            chrom                           <- driver_map_2[driver,2]
            chrom_ploidy                    <- ploidy_chrom_2[chrom]
            driver_map_2[driver,3]          <- driver_map_2[driver,3]+chrom_ploidy
        }
    }
#-----------------------------------------------Output the new genotypes
    genotype_list_ploidy_chrom[[genotype_daughter_1]]         <<- ploidy_chrom_1
    genotype_list_ploidy_allele[[genotype_daughter_1]]        <<- ploidy_allele_1
    genotype_list_ploidy_block[[genotype_daughter_1]]         <<- ploidy_block_1
    genotype_list_driver_count[genotype_daughter_1]           <<- driver_count_1
    genotype_list_driver_map[[genotype_daughter_1]]           <<- driver_map_1
    evolution_genotype_changes[[genotype_daughter_1]]         <<- c(evolution_genotype_changes[[genotype_daughter_1]],list('whole-genome-duplication'))

    genotype_list_ploidy_chrom[[genotype_daughter_2]]         <<- ploidy_chrom_2
    genotype_list_ploidy_allele[[genotype_daughter_2]]        <<- ploidy_allele_2
    genotype_list_ploidy_block[[genotype_daughter_2]]         <<- ploidy_block_2
    genotype_list_driver_count[genotype_daughter_2]           <<- driver_count_2
    genotype_list_driver_map[[genotype_daughter_2]]           <<- driver_map_2
    evolution_genotype_changes[[genotype_daughter_2]]         <<- c(evolution_genotype_changes[[genotype_daughter_2]],list('whole-genome-duplication'))
}
