# ================================================SIMULATE MISSEGREGATION
#' @export
SIMULATOR_FULL_PHASE_1_CN_missegregation <- function(genotype_to_react,
                                                     genotype_daughter_1,
                                                     genotype_daughter_2 = NULL,
                                                     chromosomes_excluded = NULL,
                                                     event = NULL) {
    #------------------------------------Find the new CN and driver profiles
    #   Find the daughter cells' current CN and driver profiles
    ploidy_chrom_1 <- genotype_list_ploidy_chrom[[genotype_daughter_1]]
    ploidy_allele_1 <- genotype_list_ploidy_allele[[genotype_daughter_1]]
    ploidy_block_1 <- genotype_list_ploidy_block[[genotype_daughter_1]]
    driver_count_1 <- genotype_list_driver_count[genotype_daughter_1]
    driver_map_1 <- genotype_list_driver_map[[genotype_daughter_1]]
    if (!is.null(genotype_daughter_2)) {
        ploidy_chrom_2 <- genotype_list_ploidy_chrom[[genotype_daughter_2]]
        ploidy_allele_2 <- genotype_list_ploidy_allele[[genotype_daughter_2]]
        ploidy_block_2 <- genotype_list_ploidy_block[[genotype_daughter_2]]
        driver_count_2 <- genotype_list_driver_count[genotype_daughter_2]
        driver_map_2 <- genotype_list_driver_map[[genotype_daughter_2]]
    }
    #   Find information about the missegregation
    if (!is.null(event)) {
        chrom <- strtoi(event[2])
        strand <- strtoi(event[3])
        if (strtoi(event[4]) == 1) {
            i_gain <- 1
        } else if (strtoi(event[4]) == -1) {
            i_gain <- 2
        }
    } else if (length(unique(chromosomes_excluded)) >= N_chromosomes) {
        return()
    } else {
        while (1) {
            #   Choose which cell to gain/lose the strand
            i_gain <- sample.int(2, size = 1)
            #   Choose the chromosome to be mis-segregated
            chrom <- sample.int(N_chromosomes, size = 1)
            # chrom <- sample.int(N_chromosomes, size = 1, prob = ploidy_chrom_1 / sum(ploidy_chrom_1))
            if (!is.null(chromosomes_excluded) & (chrom %in% chromosomes_excluded)) {
                next
            }
            if (i_gain == 1 & !is.null(genotype_daughter_2)) {
                no_strands <- ploidy_chrom_2[chrom]
            } else if (i_gain == 1 & is.null(genotype_daughter_2)) {
                no_strands <- ploidy_chrom_1[chrom]
            } else if (i_gain == 2) {
                no_strands <- ploidy_chrom_1[chrom]
            }
            if (no_strands <= 0) {
                next
            }
            #   Choose the strand to be mis-segregated
            strand <- sample.int(no_strands, size = 1)
            break
        }
    }
    #   Find all drivers located on this strand in the losing cell
    if ((i_gain == 1) & (!(is.null(genotype_daughter_2)))) {
        if (driver_count_2 > 0) {
            pos_drivers_to_move <- intersect(which((driver_map_2[, 2] == chrom)), which((driver_map_2[, 3] == strand)))
        } else {
            pos_drivers_to_move <- c()
        }
    } else if ((i_gain == 1) & ((is.null(genotype_daughter_2)) & (driver_count_1 > 0))) {
        pos_drivers_to_move <- intersect(which((driver_map_1[, 2] == chrom)), which((driver_map_1[, 3] == strand)))
    } else if ((i_gain == 2) & (driver_count_1 > 0)) {
        pos_drivers_to_move <- intersect(which((driver_map_1[, 2] == chrom)), which((driver_map_1[, 3] == strand)))
    } else {
        pos_drivers_to_move <- c()
    }
    N_drivers_to_move <- length(pos_drivers_to_move)
    #   Change the chromosome ploidy of daughter cells
    if (i_gain == 1) {
        ploidy_chrom_1[chrom] <- ploidy_chrom_1[chrom] + 1
        if (!(is.null(genotype_daughter_2))) {
            ploidy_chrom_2[chrom] <- ploidy_chrom_2[chrom] - 1
        }
    } else if (i_gain == 2) {
        if (!(is.null(genotype_daughter_2))) {
            ploidy_chrom_2[chrom] <- ploidy_chrom_2[chrom] + 1
        }
        ploidy_chrom_1[chrom] <- ploidy_chrom_1[chrom] - 1
    }
    #   Update the chromosome strand allele identities of daughter cells
    if ((i_gain == 1) & (!is.null(genotype_daughter_2))) {
        chrom_ploidy <- ploidy_chrom_1[chrom]
        ploidy_allele_1[[chrom]][[chrom_ploidy]] <- ploidy_allele_2[[chrom]][[strand]]
        chrom_ploidy <- ploidy_chrom_2[chrom]
        if (strand <= chrom_ploidy) {
            for (i_strand in strand:chrom_ploidy) {
                ploidy_allele_2[[chrom]][[i_strand]] <- ploidy_allele_2[[chrom]][[i_strand + 1]]
            }
        }
        ploidy_allele_2[[chrom]] <- ploidy_allele_2[[chrom]][-(chrom_ploidy + 1)]
    } else if ((i_gain == 1) & (is.null(genotype_daughter_2))) {
        chrom_ploidy <- ploidy_chrom_1[chrom]
        ploidy_allele_1[[chrom]][[chrom_ploidy]] <- ploidy_allele_1[[chrom]][[strand]]
    } else if ((i_gain == 2) & (!is.null(genotype_daughter_2))) {
        chrom_ploidy <- ploidy_chrom_2[chrom]
        ploidy_allele_2[[chrom]][[chrom_ploidy]] <- ploidy_allele_1[[chrom]][[strand]]
        chrom_ploidy <- ploidy_chrom_1[chrom]
        if (strand <= chrom_ploidy) {
            for (i_strand in strand:chrom_ploidy) {
                ploidy_allele_1[[chrom]][[i_strand]] <- ploidy_allele_1[[chrom]][[i_strand + 1]]
            }
        }
        ploidy_allele_1[[chrom]] <- ploidy_allele_1[[chrom]][-(chrom_ploidy + 1)]
    } else if ((i_gain == 2) & (is.null(genotype_daughter_2))) {
        chrom_ploidy <- ploidy_chrom_1[chrom]
        if (strand <= chrom_ploidy) {
            for (i_strand in strand:chrom_ploidy) {
                ploidy_allele_1[[chrom]][[i_strand]] <- ploidy_allele_1[[chrom]][[i_strand + 1]]
            }
        }
        ploidy_allele_1[[chrom]] <- ploidy_allele_1[[chrom]][-(chrom_ploidy + 1)]
    }
    #   Move the chromosome strand from losing cell to winning cell
    if ((i_gain == 1) & (!is.null(genotype_daughter_2))) {
        chrom_ploidy <- ploidy_chrom_1[chrom]
        ploidy_block_1[[chrom]][[chrom_ploidy]] <- ploidy_block_2[[chrom]][[strand]]
        chrom_ploidy <- ploidy_chrom_2[chrom]
        if (strand <= chrom_ploidy) {
            for (i_strand in strand:chrom_ploidy) {
                ploidy_block_2[[chrom]][[i_strand]] <- ploidy_block_2[[chrom]][[i_strand + 1]]
            }
        }
        ploidy_block_2[[chrom]] <- ploidy_block_2[[chrom]][-(chrom_ploidy + 1)]
    } else if ((i_gain == 1) & (is.null(genotype_daughter_2))) {
        chrom_ploidy <- ploidy_chrom_1[chrom]
        ploidy_block_1[[chrom]][[chrom_ploidy]] <- ploidy_block_1[[chrom]][[strand]]
    } else if ((i_gain == 2) & (!is.null(genotype_daughter_2))) {
        chrom_ploidy <- ploidy_chrom_2[chrom]
        ploidy_block_2[[chrom]][[chrom_ploidy]] <- ploidy_block_1[[chrom]][[strand]]
        chrom_ploidy <- ploidy_chrom_1[chrom]
        if (strand <= chrom_ploidy) {
            for (i_strand in strand:chrom_ploidy) {
                ploidy_block_1[[chrom]][[i_strand]] <- ploidy_block_1[[chrom]][[i_strand + 1]]
            }
        }
        ploidy_block_1[[chrom]] <- ploidy_block_1[[chrom]][-(chrom_ploidy + 1)]
    } else if ((i_gain == 2) & (is.null(genotype_daughter_2))) {
        chrom_ploidy <- ploidy_chrom_1[chrom]
        if (strand <= chrom_ploidy) {
            for (i_strand in strand:chrom_ploidy) {
                ploidy_block_1[[chrom]][[i_strand]] <- ploidy_block_1[[chrom]][[i_strand + 1]]
            }
        }
        ploidy_block_1[[chrom]] <- ploidy_block_1[[chrom]][-(chrom_ploidy + 1)]
    }
    #   Move the drivers from losing cell to winning cell
    if ((i_gain == 1) & (N_drivers_to_move > 0) & (!is.null(genotype_daughter_2))) {
        driver_map_new_1 <- driver_map_2[pos_drivers_to_move, ]
        if (!is.matrix(driver_map_new_1)) {
            driver_map_new_1 <- matrix(driver_map_new_1, nrow = 1)
        }
        driver_map_new_1[, 3] <- ploidy_chrom_1[chrom]
        driver_map_1 <- rbind(driver_map_1, driver_map_new_1)
        driver_map_2 <- driver_map_2[-pos_drivers_to_move, ]
        if (!is.matrix(driver_map_2)) {
            driver_map_2 <- matrix(driver_map_2, nrow = 1)
        }
    } else if ((i_gain == 1) & (N_drivers_to_move > 0) & (is.null(genotype_daughter_2))) {
        driver_map_new_1 <- driver_map_1[pos_drivers_to_move, ]
        if (!is.matrix(driver_map_new_1)) {
            driver_map_new_1 <- matrix(driver_map_new_1, nrow = 1)
        }
        driver_map_new_1[, 3] <- ploidy_chrom_1[chrom]
        driver_map_1 <- rbind(driver_map_1, driver_map_new_1)
    } else if ((i_gain == 2) & (N_drivers_to_move > 0) & (!is.null(genotype_daughter_2))) {
        driver_map_new_2 <- driver_map_1[pos_drivers_to_move, ]
        if (!is.matrix(driver_map_new_2)) {
            driver_map_new_2 <- matrix(driver_map_new_2, nrow = 1)
        }
        driver_map_new_2[, 3] <- ploidy_chrom_2[chrom]
        driver_map_2 <- rbind(driver_map_2, driver_map_new_2)
        driver_map_1 <- driver_map_1[-pos_drivers_to_move, ]
        if (!is.matrix(driver_map_1)) {
            driver_map_1 <- matrix(driver_map_1, nrow = 1)
        }
    } else if ((i_gain == 2) & (N_drivers_to_move > 0) & (is.null(genotype_daughter_2))) {
        driver_map_1 <- driver_map_1[-pos_drivers_to_move, ]
        if (!is.matrix(driver_map_1)) {
            driver_map_1 <- matrix(driver_map_1, nrow = 1)
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
    genotype_list_driver_count[genotype_daughter_1] <<- driver_count_1
    genotype_list_driver_map[[genotype_daughter_1]] <<- driver_map_1
    if (i_gain == 1) {
        loc_end <- length(evolution_genotype_changes[[genotype_daughter_1]])
        evolution_genotype_changes[[genotype_daughter_1]][[loc_end + 1]] <<- c("missegregation", chrom, strand, 1)
    } else {
        if (i_gain == 2) {
            loc_end <- length(evolution_genotype_changes[[genotype_daughter_1]])
            evolution_genotype_changes[[genotype_daughter_1]][[loc_end + 1]] <<- c("missegregation", chrom, strand, -1)
        }
    }
    if (!is.null(genotype_daughter_2)) {
        genotype_list_ploidy_chrom[[genotype_daughter_2]] <<- ploidy_chrom_2
        genotype_list_ploidy_allele[[genotype_daughter_2]] <<- ploidy_allele_2
        genotype_list_ploidy_block[[genotype_daughter_2]] <<- ploidy_block_2
        genotype_list_driver_count[genotype_daughter_2] <<- driver_count_2
        genotype_list_driver_map[[genotype_daughter_2]] <<- driver_map_2
        if (i_gain == 1) {
            loc_end <- length(evolution_genotype_changes[[genotype_daughter_2]])
            evolution_genotype_changes[[genotype_daughter_2]][[loc_end + 1]] <<- c("missegregation", chrom, strand, -1)
        } else {
            if (i_gain == 2) {
                loc_end <- length(evolution_genotype_changes[[genotype_daughter_2]])
                evolution_genotype_changes[[genotype_daughter_2]][[loc_end + 1]] <<- c("missegregation", chrom, strand, 1)
            }
        }
    }
}
