# ==================================COMPUTE THE SELECTION RATE OF A CLONE
#' @export
SIMULATOR_FULL_PHASE_1_selection_rate <- function(driver_count, driver_map, ploidy_chrom, ploidy_block, ploidy_allele) {
    #---------------------Cell is not viable if losing whole chromosomes
    #--------------------------------or exceeding maximum average ploidy
    vec_CN_all <- c()
    for (chrom in 1:N_chromosomes) {
        vec_CN <- rep(0, vec_CN_block_no[chrom])
        no_strands <- ploidy_chrom[chrom]
        if (no_strands > 0) {
            for (strand in 1:no_strands) {
                vec_CN <- vec_CN + ploidy_block[[chrom]][[strand]]
            }
        }
        vec_CN_all <- c(vec_CN_all, vec_CN)
    }
    ploidy <- mean(vec_CN_all)
    if (mean(vec_CN_all) > bound_average_ploidy) {
        clone_selection_rate <- 0
        return(clone_selection_rate)
    }
    if (max(vec_CN_all) > bound_maximum_CN) {
        clone_selection_rate <- 0
        return(clone_selection_rate)
    }
    if (max(vec_CN_all / ploidy) > bound_maximum_CN_normalized) {
        clone_selection_rate <- 0
        return(clone_selection_rate)
    }
    if (length(which(vec_CN_all == 0)) > bound_homozygosity) {
        clone_selection_rate <- 0
        return(clone_selection_rate)
    }
    # #-----Cell is not viable if exceeding maximum length of homozygosity
    # L_homozygosity <- 0
    # for (chrom in 1:N_chromosomes) {
    #     no_strands <- ploidy_chrom[chrom]
    #     vec_CN <- rep(0, vec_CN_block_no[chrom])
    #     for (strand in 1:no_strands) {
    #         vec_CN <- vec_CN + ploidy_block[[chrom]][[strand]]
    #     }
    #     L_homozygosity <- L_homozygosity + length(which(vec_CN == 0))
    #     if (L_homozygosity > bound_homozygosity) {
    #         clone_selection_rate <- 0
    #         return(clone_selection_rate)
    #     }
    # }
    #---------------Cell is not viable if exceeding maximum driver count
    if (driver_count > 0) {
        driver_count_unique <- unique(driver_map[, 1])
        if (length(driver_count_unique) > bound_driver) {
            clone_selection_rate <- 0
            return(clone_selection_rate)
        }
    }
    #--------------------------Compute selection rates for viable clones
    if (selection_model == "chrom-arm-selection-old") {
        #--------------------------------------------Find average ploidy
        ploidy <- round(mean(vec_CN_all))
        #-----------------------------------------Compute selection rate
        #   Find average CN per chromosome arm
        chrom_arm_library_copy <- chrom_arm_library
        chrom_arm_library_copy$cn <- 0
        for (i_arm in 1:nrow(chrom_arm_library_copy)) {
            chrom <- which(vec_chromosome_id == chrom_arm_library_copy$Chromosome[i_arm])
            start <- chrom_arm_library_copy$Bin_start[i_arm]
            end <- chrom_arm_library_copy$Bin_end[i_arm]
            no_strands <- ploidy_chrom[chrom]
            if (no_strands == 0) {
                cn <- 0
            } else {
                vec_cn <- ploidy_block[[chrom]][[1]]
                if (no_strands > 1) {
                    for (strand in 2:no_strands) {
                        vec_cn <- vec_cn + ploidy_block[[chrom]][[strand]]
                    }
                }
                cn <- round(mean(vec_cn[start:end]))
            }
            chrom_arm_library_copy$cn[i_arm] <- cn
        }
        clone_selection_rate <- prod(chrom_arm_library_copy$s_rate^(chrom_arm_library_copy$cn / ploidy))
    } else if (selection_model == "chrom-arm-selection") {
        #--------------------------------------------Find average ploidy
        ploidy <- mean(vec_CN_all)
        #-----------------------------------------Compute selection rate
        #   Find average CN per chromosome arm
        chrom_arm_library_copy <- chrom_arm_library
        chrom_arm_library_copy$cn <- 0
        for (i_arm in 1:nrow(chrom_arm_library_copy)) {
            chrom <- which(vec_chromosome_id == chrom_arm_library_copy$Chromosome[i_arm])
            start <- chrom_arm_library_copy$Bin_start[i_arm]
            end <- chrom_arm_library_copy$Bin_end[i_arm]
            no_strands <- ploidy_chrom[chrom]
            if (no_strands == 0) {
                cn <- 0
            } else {
                vec_cn <- ploidy_block[[chrom]][[1]]
                if (no_strands > 1) {
                    for (strand in 2:no_strands) {
                        vec_cn <- vec_cn + ploidy_block[[chrom]][[strand]]
                    }
                }
                cn <- round(mean(vec_cn[start:end]))
            }
            chrom_arm_library_copy$cn[i_arm] <- cn
        }
        clone_selection_rate <- prod(chrom_arm_library_copy$s_rate^(chrom_arm_library_copy$cn / ploidy))

        # #   ???
        # #   ???
        # #   ???
        # # ploidy <- 2 * round(mean(chrom_arm_library_copy$cn) / 2)
        # # ploidy <- round(mean(chrom_arm_library_copy$cn))
        # # ploidy <- (mean(chrom_arm_library_copy$cn))
        #
        # #   ???
        # #   ???
        # #   ???
        #
        # clone_selection_rate <- prod(chrom_arm_library_copy$s_rate^(chrom_arm_library_copy$cn / ploidy))
        # if (TMPTMP == 1) {
        #     cat("CN profile:        ", chrom_arm_library_copy$cn, "\n")
        #     cat("Ploidy:            ", ploidy, "\n")
        #     cat("Real ploidy:       ", (mean(chrom_arm_library_copy$cn)), "\n")
        #     # cat("Minimum CN:        ", min(chrom_arm_library_copy$cn), "\n")
        #     # cat("Maximum CN:        ", max(chrom_arm_library_copy$cn), "\n")
        #     cn_unique <- sort(unique(chrom_arm_library_copy$cn))
        #     for (i in 1:length(cn_unique)) {
        #         cat("Arms with CN=", cn_unique[i], ":  ", length(which(chrom_arm_library_copy$cn == cn_unique[i])), "\n")
        #     }
        #     cat("\n")
        #     cat("Fitness:           ", clone_selection_rate, "\n")
        #     cat("Misseg count:      ", sum(abs(chrom_arm_library_copy$cn - ploidy)), "\n")
        #
        #
        #
        #     # if (ploidy > 3) {
        #     #     cn <- chrom_arm_library_copy$cn
        #     #     cn[which(cn == 2)] <- 1
        #     #     ploidy <- (mean(cn))
        #     #     if (max(cn / ploidy) > 2) {
        #     #         tmp_clone_selection_rate <- 0
        #     #     }
        #     #     tmp_clone_selection_rate <- prod(chrom_arm_library_copy$s_rate^(cn / ploidy))
        #     #     cat("\n\n\nWHAT COULD HAVE BEEN:\n")
        #     #     cat("CN profile:        ", cn, "\n")
        #     #     cat("Ploidy:            ", ploidy, "\n")
        #     #     cat("Real ploidy:       ", (mean(cn)), "\n")
        #     #     # cat("Minimum CN:        ", min(cn), "\n")
        #     #     # cat("Maximum CN:        ", max(cn), "\n")
        #     #     cn_unique <- sort(unique(cn))
        #     #     for (i in 1:length(cn_unique)) {
        #     #         cat("Arms with CN=", cn_unique[i], ":  ", length(which(cn == cn_unique[i])), "\n")
        #     #     }
        #     #     cat("\n")
        #     #     cat("Fitness:           ", tmp_clone_selection_rate, "\n")
        #     #     cat("Misseg count:      ", sum(abs(cn - ploidy)), "\n")
        #     # }
        # }
    } else if (selection_model == "driver-gene-selection") {
        #--------------------------------------------Find average ploidy
        ploidy <- max(1, round(mean(vec_CN_all)))
        #--If driver library is empty, then viable cells have sel rate 1
        if (nrow(driver_library) == 0) {
            clone_selection_rate <- 1
            return(clone_selection_rate)
        }
        #-----------------------------------------Compute selection rate
        #   Find WT and MUT allele counts for each driver
        driver_library_copy <- driver_library
        driver_library_copy$Copy_WT <- 0
        driver_library_copy$Copy_MUT <- 0
        for (i_driver in 1:nrow(driver_library_copy)) {
            chrom <- driver_library_copy$Chromosome[i_driver]
            block <- driver_library_copy$Bin[i_driver]
            no_strands <- ploidy_chrom[chrom]
            driver_copy <- 0
            if (no_strands > 0) {
                for (strand in 1:no_strands) {
                    driver_copy <- driver_copy + ploidy_block[[chrom]][[strand]][block]
                }
            }
            driver_library_copy$Copy_WT[i_driver] <- driver_copy
        }
        if (driver_count >= 1) {
            for (i_driver in 1:driver_count) {
                driver_ID <- driver_map[i_driver, 1]
                driver_library_copy$Copy_MUT[driver_ID] <- driver_library_copy$Copy_MUT[driver_ID] + 1
                driver_library_copy$Copy_WT[driver_ID] <- driver_library_copy$Copy_WT[driver_ID] - 1
            }
        }
        #   Compute selection rate
        clone_selection_rate <- prod(driver_library_copy$s_rate_WT^(2 * driver_library_copy$Copy_WT / ploidy)) *
            prod(driver_library_copy$s_rate_MUT^(2 * driver_library_copy$Copy_MUT / ploidy))
    } else if (selection_model == "chrom-arm-and-driver-gene-selection") {
        #--------------------------------------------Find average ploidy
        ploidy <- max(1, round(mean(vec_CN_all)))
        #-----------------------------------------Compute selection rate
        #   Find average CN per chromosome arm
        chrom_arm_library_copy <- chrom_arm_library
        chrom_arm_library_copy$cn <- 0
        for (i_arm in 1:nrow(chrom_arm_library_copy)) {
            chrom <- which(vec_chromosome_id == chrom_arm_library_copy$Chromosome[i_arm])
            start <- chrom_arm_library_copy$Bin_start[i_arm]
            end <- chrom_arm_library_copy$Bin_end[i_arm]
            no_strands <- ploidy_chrom[chrom]
            if (no_strands == 0) {
                cn <- 0
            } else {
                vec_cn <- ploidy_block[[chrom]][[1]]
                if (no_strands > 1) {
                    for (strand in 2:no_strands) {
                        vec_cn <- vec_cn + ploidy_block[[chrom]][[strand]]
                    }
                }
                cn <- round(mean(vec_cn[start:end]))
            }
            chrom_arm_library_copy$cn[i_arm] <- cn
        }
        #   If driver library is empty, then viable cells have sel rate 1
        if (nrow(driver_library) == 0) {
            clone_selection_rate <- prod(chrom_arm_library_copy$s_rate^(chrom_arm_library_copy$cn / ploidy))
            return(clone_selection_rate)
        }
        #   Find WT and MUT allele counts for each driver
        driver_library_copy <- driver_library
        driver_library_copy$Copy_WT <- 0
        driver_library_copy$Copy_MUT <- 0
        for (i_driver in 1:nrow(driver_library_copy)) {
            chrom <- which(vec_chromosome_id == driver_library_copy$Chromosome[i_driver])
            # chrom <- driver_library_copy$Chromosome[i_driver]
            block <- driver_library_copy$Bin[i_driver]
            no_strands <- ploidy_chrom[chrom]
            driver_copy <- 0
            if (no_strands >= 1) {
                for (strand in 1:no_strands) {
                    driver_copy <- driver_copy + ploidy_block[[chrom]][[strand]][block]
                }
            }
            driver_library_copy$Copy_WT[i_driver] <- driver_copy
        }
        if (driver_count >= 1) {
            for (i_driver in 1:driver_count) {
                driver_ID <- driver_map[i_driver, 1]
                driver_library_copy$Copy_MUT[driver_ID] <- driver_library_copy$Copy_MUT[driver_ID] + 1
                driver_library_copy$Copy_WT[driver_ID] <- driver_library_copy$Copy_WT[driver_ID] - 1
            }
        }
        #   Compute selection rate
        clone_selection_rate <- prod(chrom_arm_library_copy$s_rate^(chrom_arm_library_copy$cn / ploidy)) *
            prod(driver_library_copy$s_rate_WT^(2 * driver_library_copy$Copy_WT / ploidy)) *
            prod(driver_library_copy$s_rate_MUT^(2 * driver_library_copy$Copy_MUT / ploidy))
    } else if (selection_model == "chrom-arm-and-driver-gene-selection-diploid-base") {
        #--------------------------------------------Find average ploidy
        ploidy <- 2
        #-----------------------------------------Compute selection rate
        #   Find average CN per chromosome arm
        chrom_arm_library_copy <- chrom_arm_library
        chrom_arm_library_copy$cn <- 0
        for (i_arm in 1:nrow(chrom_arm_library_copy)) {
            chrom <- which(vec_chromosome_id == chrom_arm_library_copy$Chromosome[i_arm])
            start <- chrom_arm_library_copy$Bin_start[i_arm]
            end <- chrom_arm_library_copy$Bin_end[i_arm]
            no_strands <- ploidy_chrom[chrom]
            if (no_strands == 0) {
                cn <- 0
            } else {
                vec_cn <- ploidy_block[[chrom]][[1]]
                if (no_strands > 1) {
                    for (strand in 2:no_strands) {
                        vec_cn <- vec_cn + ploidy_block[[chrom]][[strand]]
                    }
                }
                cn <- round(mean(vec_cn[start:end]))
            }
            chrom_arm_library_copy$cn[i_arm] <- cn
        }
        #   If driver library is empty, then viable cells have sel rate 1
        if (nrow(driver_library) == 0) {
            clone_selection_rate <- prod(chrom_arm_library_copy$s_rate^(chrom_arm_library_copy$cn / ploidy))
            return(clone_selection_rate)
        }
        #   Find WT and MUT allele counts for each driver
        driver_library_copy <- driver_library
        driver_library_copy$Copy_WT <- 0
        driver_library_copy$Copy_MUT <- 0
        for (i_driver in 1:nrow(driver_library_copy)) {
            chrom <- which(vec_chromosome_id == driver_library_copy$Chromosome[i_driver])
            # chrom <- driver_library_copy$Chromosome[i_driver]
            block <- driver_library_copy$Bin[i_driver]
            no_strands <- ploidy_chrom[chrom]
            driver_copy <- 0
            for (strand in 1:no_strands) {
                driver_copy <- driver_copy + ploidy_block[[chrom]][[strand]][block]
            }
            driver_library_copy$Copy_WT[i_driver] <- driver_copy
        }
        if (driver_count >= 1) {
            for (i_driver in 1:driver_count) {
                driver_ID <- driver_map[i_driver, 1]
                driver_library_copy$Copy_MUT[driver_ID] <- driver_library_copy$Copy_MUT[driver_ID] + 1
                driver_library_copy$Copy_WT[driver_ID] <- driver_library_copy$Copy_WT[driver_ID] - 1
            }
        }
        #   Compute selection rate
        clone_selection_rate <- prod(chrom_arm_library_copy$s_rate^(chrom_arm_library_copy$cn / ploidy)) *
            prod(driver_library_copy$s_rate_WT^(2 * driver_library_copy$Copy_WT / ploidy)) *
            prod(driver_library_copy$s_rate_MUT^(2 * driver_library_copy$Copy_MUT / ploidy))
    } else if (selection_model == "ancient") {
        #--If driver library is empty, then viable cells have sel rate 1
        if (nrow(driver_library) == 0) {
            clone_selection_rate <- 1
            return(clone_selection_rate)
        }
        #-----------------------------------------Compute selection rate
        driver_library_copy <- driver_library
        driver_library_copy$Copy_WT <- 0
        driver_library_copy$Copy_MUT <- 0
        #   Find WT and MUT allele counts for each driver
        for (i_driver in 1:nrow(driver_library_copy)) {
            chrom <- driver_library_copy$Chromosome[i_driver]
            block <- driver_library_copy$Bin[i_driver]
            no_strands <- ploidy_chrom[chrom]
            driver_copy <- 0
            for (strand in 1:no_strands) {
                driver_copy <- driver_copy + ploidy_block[[chrom]][[strand]][block]
            }
            driver_library_copy$Copy_WT[i_driver] <- driver_copy
        }
        if (driver_count >= 1) {
            for (i_driver in 1:driver_count) {
                driver_ID <- driver_map[i_driver, 1]
                driver_library_copy$Copy_MUT[driver_ID] <- driver_library_copy$Copy_MUT[driver_ID] + 1
                driver_library_copy$Copy_WT[driver_ID] <- driver_library_copy$Copy_WT[driver_ID] - 1
            }
        }
        #   Compute selection rate
        clone_selection_rate <- (s_normalization^mean(vec_CN_all)) *
            prod(driver_library_copy$s_rate_WT^driver_library_copy$Copy_WT) *
            prod(driver_library_copy$s_rate_MUT^driver_library_copy$Copy_MUT)
    } else if (selection_model == "old") {
        #--If driver library is empty, then viable cells have sel rate 1
        if (nrow(driver_library) == 0) {
            clone_selection_rate <- 1
            return(clone_selection_rate)
        }
        #-----------------------------------------Compute selection rate
        driver_library_copy <- driver_library
        driver_library_copy$Copy_WT <- 0
        driver_library_copy$Copy_MUT <- 0
        #   Find WT and MUT allele counts for each driver
        for (i_driver in 1:nrow(driver_library_copy)) {
            chrom <- driver_library_copy$Chromosome[i_driver]
            block <- driver_library_copy$Bin[i_driver]
            no_strands <- ploidy_chrom[chrom]
            driver_copy <- 0
            for (strand in 1:no_strands) {
                driver_copy <- driver_copy + ploidy_block[[chrom]][[strand]][block]
            }
            driver_library_copy$Copy_WT[i_driver] <- driver_copy
        }
        if (driver_count >= 1) {
            for (i_driver in 1:driver_count) {
                driver_ID <- driver_map[i_driver, 1]
                driver_library_copy$Copy_MUT[driver_ID] <- driver_library_copy$Copy_MUT[driver_ID] + 1
                driver_library_copy$Copy_WT[driver_ID] <- driver_library_copy$Copy_WT[driver_ID] - 1
            }
        }
        #   Compute selection rate
        clone_selection_rate <- prod(driver_library_copy$s_rate_WT^driver_library_copy$Copy_WT) *
            prod(driver_library_copy$s_rate_MUT^driver_library_copy$Copy_MUT)
    }
    return(clone_selection_rate)
}
