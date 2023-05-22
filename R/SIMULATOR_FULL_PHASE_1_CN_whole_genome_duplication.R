# ======================================SIMULATE WHOLE GENOME DUPLICATION
#' @export
SIMULATOR_FULL_PHASE_1_CN_whole_genome_duplication <- function(genotype_to_react,
                                                               genotype_daughter_1,
                                                               genotype_daughter_2 = NULL) {
    #------------------------------------Find the new CN and driver profiles
    #   Find the daughter cells' current CN and driver profiles
    ploidy_chrom_1 <- genotype_list_ploidy_chrom[[genotype_daughter_1]]
    ploidy_allele_1 <- genotype_list_ploidy_allele[[genotype_daughter_1]]
    ploidy_block_1 <- genotype_list_ploidy_block[[genotype_daughter_1]]
    WGD_count_1 <- genotype_list_WGD_count[genotype_daughter_1]
    driver_count_1 <- genotype_list_driver_count[genotype_daughter_1]
    driver_map_1 <- genotype_list_driver_map[[genotype_daughter_1]]
    if (!is.null(genotype_daughter_2)) {
        ploidy_chrom_2 <- genotype_list_ploidy_chrom[[genotype_daughter_2]]
        ploidy_allele_2 <- genotype_list_ploidy_allele[[genotype_daughter_2]]
        ploidy_block_2 <- genotype_list_ploidy_block[[genotype_daughter_2]]
        WGD_count_2 <- genotype_list_WGD_count[genotype_daughter_2]
        driver_count_2 <- genotype_list_driver_count[genotype_daughter_2]
        driver_map_2 <- genotype_list_driver_map[[genotype_daughter_2]]
    }
    #   Change the WGD count of daughter cells
    WGD_count_1 <- WGD_count_1 + 1
    if (!is.null(genotype_daughter_2)) {
        WGD_count_2 <- WGD_count_2 + 1
    }
    #   Change the chromosome ploidy of daughter cells
    ploidy_chrom_1 <- 2 * ploidy_chrom_1
    if (!is.null(genotype_daughter_2)) {
        ploidy_chrom_2 <- 2 * ploidy_chrom_2
    }
    #   Update the chromosome strand allele identities of daughter cells
    for (chrom in 1:N_chromosomes) {
        chrom_ploidy <- ploidy_chrom_1[chrom]
        if (chrom_ploidy >= 1) {
            for (strand in 1:(chrom_ploidy / 2)) {
                ploidy_allele_1[[chrom]][[chrom_ploidy / 2 + strand]] <- ploidy_allele_1[[chrom]][[strand]]
            }
        }
        if (!is.null(genotype_daughter_2)) {
            chrom_ploidy <- ploidy_chrom_2[chrom]
            if (chrom_ploidy >= 1) {
                for (strand in 1:(chrom_ploidy / 2)) {
                    ploidy_allele_2[[chrom]][[chrom_ploidy / 2 + strand]] <- ploidy_allele_2[[chrom]][[strand]]
                }
            }
        }
    }
    #   Multiply the chromosome strands in each daughter cell
    for (chrom in 1:N_chromosomes) {
        chrom_ploidy <- ploidy_chrom_1[chrom]
        if (chrom_ploidy >= 1) {
            for (strand in 1:(chrom_ploidy / 2)) {
                ploidy_block_1[[chrom]][[chrom_ploidy / 2 + strand]] <- ploidy_block_1[[chrom]][[strand]]
            }
        }
        if (!is.null(genotype_daughter_2)) {
            chrom_ploidy <- ploidy_chrom_2[chrom]
            if (chrom_ploidy >= 1) {
                for (strand in 1:(chrom_ploidy / 2)) {
                    ploidy_block_2[[chrom]][[chrom_ploidy / 2 + strand]] <- ploidy_block_2[[chrom]][[strand]]
                }
            }
        }
    }
    #   Multiply the drivers in each daughter cell
    if (driver_count_1 > 0) {
        driver_map_new_1 <- driver_map_1
        for (driver in 1:nrow(driver_map_new_1)) {
            chrom <- driver_map_1[driver, 2]
            chrom_ploidy <- ploidy_chrom_1[chrom]
            driver_map_new_1[driver, 3] <- driver_map_new_1[driver, 3] + chrom_ploidy / 2
        }
        driver_map_1 <- rbind(driver_map_1, driver_map_new_1)
    }
    if (!is.null(genotype_daughter_2)) {
        if (driver_count_2 > 0) {
            driver_map_new_2 <- driver_map_2
            for (driver in 1:nrow(driver_map_new_2)) {
                chrom <- driver_map_2[driver, 2]
                chrom_ploidy <- ploidy_chrom_2[chrom]
                driver_map_new_2[driver, 3] <- driver_map_new_2[driver, 3] + chrom_ploidy / 2
            }
            driver_map_2 <- rbind(driver_map_2, driver_map_new_2)
        }
    }
    #   Change the driver count in each daughter cell
    driver_unique_1 <- unique(driver_map_1[, 1])
    driver_unique_1 <- driver_unique_1[driver_unique_1 != 0]
    driver_count_1 <- length(driver_unique_1)
    if (!is.null(genotype_daughter_2)) {
        driver_unique_2 <- unique(driver_map_2[, 1])
        driver_unique_2 <- driver_unique_2[driver_unique_2 != 0]
        driver_count_2 <- length(driver_unique_2)
    }
    #-----------------------------------------------Output the new genotypes
    genotype_list_ploidy_chrom[[genotype_daughter_1]] <<- ploidy_chrom_1
    genotype_list_ploidy_allele[[genotype_daughter_1]] <<- ploidy_allele_1
    genotype_list_ploidy_block[[genotype_daughter_1]] <<- ploidy_block_1
    genotype_list_WGD_count[genotype_daughter_1] <<- WGD_count_1
    genotype_list_driver_count[genotype_daughter_1] <<- driver_count_1
    genotype_list_driver_map[[genotype_daughter_1]] <<- driver_map_1
    loc_end <- length(evolution_genotype_changes[[genotype_daughter_1]])
    evolution_genotype_changes[[genotype_daughter_1]][[loc_end + 1]] <<- c("whole-genome-duplication")
    if (!is.null(genotype_daughter_2)) {
        genotype_list_ploidy_chrom[[genotype_daughter_2]] <<- ploidy_chrom_2
        genotype_list_ploidy_allele[[genotype_daughter_2]] <<- ploidy_allele_2
        genotype_list_ploidy_block[[genotype_daughter_2]] <<- ploidy_block_2
        genotype_list_WGD_count[genotype_daughter_2] <<- WGD_count_2
        genotype_list_driver_count[genotype_daughter_2] <<- driver_count_2
        genotype_list_driver_map[[genotype_daughter_2]] <<- driver_map_2
        loc_end <- length(evolution_genotype_changes[[genotype_daughter_2]])
        evolution_genotype_changes[[genotype_daughter_2]][[loc_end + 1]] <<- c("whole-genome-duplication")
    }
}
