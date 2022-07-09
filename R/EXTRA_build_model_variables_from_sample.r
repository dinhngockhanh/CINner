# ===============================PRODUCE MODEL VARIABLE FILES FROM SAMPLE
#' @export
EXTRA_build_model_variables_from_sample <- function(model, package_output) {
    #---------------------------------Input the clonal populations in sample
    package_clonal_evolution <- package_output[[1]]

    package_sample <- package_output[[2]]

    genotype_list_ploidy_chrom <- package_clonal_evolution[[5]]
    genotype_list_ploidy_block <- package_clonal_evolution[[6]]
    genotype_list_ploidy_allele <- package_clonal_evolution[[7]]
    genotype_list_driver_count <- package_clonal_evolution[[8]]
    genotype_list_driver_map <- package_clonal_evolution[[9]]

    sample_clone_ID <- package_sample[[3]]
    #----------------------Output CN-driver profiles as model variable files
    #---Find all unique clones in sample
    all_clones_ID <- unique(sample_clone_ID)
    all_clones_population <- rep(0, length(all_clones_ID))
    for (i in 1:length(all_clones_ID)) {
        all_clones_population[i] <- length(which(sample_clone_ID == all_clones_ID[i]))
    }
    N_all_clones <- length(all_clones_ID)
    #---Create and output model variables - copy number state
    TABLE_CLONAL_COPY_NUMBER_PROFILES <- data.frame()
    # HEADER_CLONAL_COPY_NUMBER_PROFILES  <- c()
    #   Initialize table of CN profiles with chromosome and bin positions
    vec_chrom <- c()
    vec_bin <- c()
    for (chrom in 1:N_chromosomes) {
        bin_count <- vec_CN_block_no[chrom]
        vec_chrom <- c(vec_chrom, rep(chrom, bin_count))
        vec_bin <- c(vec_bin, (1:bin_count))
    }
    TABLE_CLONAL_COPY_NUMBER_PROFILES <- data.frame(vec_chrom, vec_bin)
    colnames(TABLE_CLONAL_COPY_NUMBER_PROFILES) <- c("Chromosome", "Bin")
    #   Update table of CN profiles with each clone found in sample
    for (clone in 1:N_all_clones) {
        clone_ID <- all_clones_ID[clone]
        #       Get the CN profile of this clone
        ploidy_chrom <- genotype_list_ploidy_chrom[[clone_ID]]
        ploidy_block <- genotype_list_ploidy_block[[clone_ID]]
        ploidy_allele <- genotype_list_ploidy_allele[[clone_ID]]
        #       Find maximum strand count for any chromosome
        max_no_strands <- max(ploidy_chrom)
        #       Translate CN profile of this clone into table
        TABLE_CLONE_CURRENT <- data.frame()
        N_row <- 0
        for (chrom in 1:N_chromosomes) {
            no_strands <- ploidy_chrom[chrom]
            bin_count <- vec_CN_block_no[chrom]
            for (bin in 1:bin_count) {
                N_row <- N_row + 1
                for (strand in 1:no_strands) {
                    unit_count <- ploidy_block[[chrom]][[strand]][bin]
                    allele <- ""
                    for (unit in 1:unit_count) {
                        allele <- paste(allele, intToUtf8(ploidy_allele[[chrom]][[strand]][unit, bin] + 64), sep = "")
                    }
                    if (allele == "") {
                        allele <- "NA"
                    }
                    TABLE_CLONE_CURRENT[N_row, strand] <- allele
                }
                if (strand == max_no_strands) {
                    next
                }
                for (strand in (no_strands + 1):max_no_strands) {
                    TABLE_CLONE_CURRENT[N_row, strand] <- "NA"
                }
            }
        }
        vec_header <- c()
        for (strand in 1:max_no_strands) {
            vec_header <- c(vec_header, paste("Clone_", clone, "_strand_", strand, sep = ""))
        }
        colnames(TABLE_CLONE_CURRENT) <- vec_header
        TABLE_CLONAL_COPY_NUMBER_PROFILES <- cbind(TABLE_CLONAL_COPY_NUMBER_PROFILES, TABLE_CLONE_CURRENT)
    }
    #   Output model variables - copy number state
    filename <- paste(model, "-input-initial-cn-profiles.csv", sep = "")
    write.csv(TABLE_CLONAL_COPY_NUMBER_PROFILES, filename, row.names = FALSE)
    #---Create and output model variables - other information
    vec_clone <- (1:N_all_clones)
    vec_cell_count <- all_clones_population
    vec_drivers <- c()
    for (clone in 1:N_all_clones) {
        DRIVERS <- ""
        driver_map <- genotype_list_driver_map[[clone]]
        for (row in 1:nrow(driver_map)) {
            driver_index <- driver_map[row, 1]
            driver_ID <- driver_library$Gene_ID[driver_index]
            driver_strand <- driver_map[row, 3]
            driver_unit <- driver_map[row, 5]
            if (DRIVERS != "") {
                DRIVERS <- paste(DRIVERS, ";", sep = "")
            }
            DRIVERS <- paste(DRIVERS, driver_ID, "_strand", driver_strand, "_unit", driver_unit, sep = "")
        }
        vec_drivers <- c(vec_drivers, DRIVERS)
    }
    #   Output model variables - other information
    CLONAL_OTHERS <- data.frame(vec_clone, vec_cell_count, vec_drivers)
    colnames(CLONAL_OTHERS) <- c("Clone", "Cell_count", "Drivers")
    filename <- paste(model, "-input-initial-others.csv", sep = "")
    write.csv(CLONAL_OTHERS, filename, row.names = FALSE)
}
