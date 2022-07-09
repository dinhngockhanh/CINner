#' @export
GET_clone_from_simulation <- function(model_variables = list(),
                                      cell_count = 0,
                                      old_simulation = "",
                                      clone = "") {
    TABLE_CHROMOSOME_CN_INFO <- model_variables$cn_info
    TABLE_CANCER_GENES <- model_variables$driver_library
    #--------------------------------------Input the previous simulation
    filename <- paste(old_simulation, ".rda", sep = "")
    load(filename)
    table_clone_ID_vs_letter <- simulation$sample$table_clone_ID_vs_letters
    genotype_list_ploidy_chrom <- simulation$clonal_evolution$genotype_list_ploidy_chrom
    genotype_list_ploidy_block <- simulation$clonal_evolution$genotype_list_ploidy_block
    genotype_list_ploidy_allele <- simulation$clonal_evolution$genotype_list_ploidy_allele

    genotype_list_driver_count <- simulation$clonal_evolution$genotype_list_driver_count
    genotype_list_driver_map <- simulation$clonal_evolution$genotype_list_driver_map
    #---------------------------------------Get requested clone's record
    clone_ID <- table_clone_ID_vs_letter$Clone_ID_number[which(table_clone_ID_vs_letter$Clone_ID_letter == clone)]

    clone_ploidy_chrom <- genotype_list_ploidy_chrom[[clone_ID]]
    clone_ploidy_block <- genotype_list_ploidy_block[[clone_ID]]
    clone_ploidy_allele <- genotype_list_ploidy_allele[[clone_ID]]

    clone_driver_count <- genotype_list_driver_count[[clone_ID]]
    clone_driver_map <- genotype_list_driver_map[[clone_ID]]
    #-----------------------------Build requested clone's driver profile
    driver_profile <- list()
    if (clone_driver_count > 0) {
        for (driver in 1:clone_driver_count) {
            driver_ID <- TABLE_CANCER_GENES$Gene_ID[clone_driver_map[driver, 1]]
            strand <- clone_driver_map[driver, 3]
            unit <- clone_driver_map[driver, 5]
            driver_profile[[driver]] <- list(strand = strand, unit = unit, driver = driver_ID)
        }
    }
    #---------------------------------Build requested clone's CN profile
    CN_profile <- data.frame(matrix(ncol = 5, nrow = 0))
    current_row <- 0
    columns <- c("Chromosome", "Strand", "Bin_start", "Bin_end", "Allele")
    colnames(CN_profile) <- columns
    #   Update clone's CN profile by chrom/strand
    for (i_chrom in 1:N_chromosomes) {
        chrom <- TABLE_CHROMOSOME_CN_INFO$Chromosome[i_chrom]
        no_bins <- vec_CN_block_no[i_chrom]
        no_strands <- clone_ploidy_chrom[i_chrom]
        for (strand in 1:no_strands) {
            strand_ploidy_block <- clone_ploidy_block[[i_chrom]][[strand]]
            strand_ploidy_allele <- clone_ploidy_allele[[i_chrom]][[strand]]
            strand_char_allele <- TRANSLATE_allele(strand_ploidy_block, strand_ploidy_allele)

            current_row <- current_row + 1
            current_allele <- strand_char_allele[[1]]

            CN_profile[current_row, 1] <- chrom
            CN_profile[current_row, 2:4] <- c(strand, 1, 1)
            CN_profile[current_row, 5] <- current_allele

            for (bin in 2:no_bins) {
                next_allele <- strand_char_allele[[bin]]
                if (next_allele != current_allele) {
                    current_row <- current_row + 1
                    current_allele <- next_allele

                    CN_profile[current_row, 1] <- chrom
                    CN_profile[current_row, 2:4] <- c(strand, bin, bin)
                    CN_profile[current_row, 5] <- current_allele
                } else {
                    CN_profile$Bin_end[current_row] <- bin
                }
            }
        }
    }
    #-----------------------------Put together requested clone's profile
    clone <- list()
    clone$CN_profile <- CN_profile
    clone$driver_profile <- driver_profile
    return(clone)
}
#' @export
TRANSLATE_allele <- function(ploidy_block, ploidy_allele) {
    ploidy_char_allele <- vector("list", length = length(ploidy_block))
    for (bin in 1:length(ploidy_block)) {
        char_allele <- ""
        no_units <- ploidy_block[bin]
        if (no_units > 0) {
            for (i in 1:no_units) {
                char_allele <- paste(char_allele, LETTERS[ploidy_allele[i, bin]], sep = "")
            }
        }
        ploidy_char_allele[bin] <- char_allele
    }
    return(ploidy_char_allele)
}
