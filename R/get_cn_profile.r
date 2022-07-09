#' @export
get_cn_profile <- function(package_clonal_evolution, clone_ID) {
    genotype_list_ploidy_chrom <- package_clonal_evolution$genotype_list_ploidy_chrom
    genotype_list_ploidy_block <- package_clonal_evolution$genotype_list_ploidy_block
    genotype_list_ploidy_allele <- package_clonal_evolution$genotype_list_ploidy_allele

    ploidy_chrom <- genotype_list_ploidy_chrom[[clone_ID]]
    ploidy_block <- genotype_list_ploidy_block[[clone_ID]]
    ploidy_allele <- genotype_list_ploidy_allele[[clone_ID]]
    #   Build the CN profile in SIGNALS style for the clone
    vec_clone_chr <- c()
    vec_clone_start <- c()
    vec_clone_end <- c()
    vec_clone_copy <- c()
    vec_clone_state <- c()
    vec_clone_Min <- c()
    vec_clone_Maj <- c()
    for (chrom in 1:N_chromosomes) {
        chrom_block_count <- vec_CN_block_no[chrom]
        chrom_ploidy <- ploidy_chrom[chrom]
        #   Find location information of each chromosome block
        vec_chr <- rep(vec_chromosome_id[chrom], 1, chrom_block_count)
        vec_start <- seq(0, size_CN_block_DNA * (chrom_block_count - 1), by = size_CN_block_DNA) + 1
        vec_end <- seq(size_CN_block_DNA, size_CN_block_DNA * chrom_block_count, by = size_CN_block_DNA)
        #   Find CN counts for each allele of each chromosome block
        vec_Allele_1 <- rep(0, 1, chrom_block_count)
        vec_Allele_2 <- rep(0, 1, chrom_block_count)
        if (chrom_ploidy >= 1) {
            for (strand in 1:chrom_ploidy) {
                mat_allele <- ploidy_allele[[chrom]][[strand]]
                for (CN_row in 1:nrow(mat_allele)) {
                    list_1 <- which(mat_allele[CN_row, ] == 1)
                    vec_Allele_1[list_1] <- vec_Allele_1[list_1] + 1
                    list_2 <- which(mat_allele[CN_row, ] == 2)
                    vec_Allele_2[list_2] <- vec_Allele_2[list_2] + 1
                }
            }
        }
        #   Find Major/Minor CN counts of each chromosome block
        if (mean(vec_Allele_1) <= mean(vec_Allele_2)) {
            vec_Min <- vec_Allele_1
            vec_Maj <- vec_Allele_2
        } else {
            vec_Min <- vec_Allele_2
            vec_Maj <- vec_Allele_1
        }
        #   Find total CN count of each chromosome block
        vec_copy <- vec_Min + vec_Maj
        vec_state <- vec_copy
        #   Update the CN information of the clone
        vec_clone_chr <- c(vec_clone_chr, vec_chr)
        vec_clone_start <- c(vec_clone_start, vec_start)
        vec_clone_end <- c(vec_clone_end, vec_end)
        vec_clone_copy <- c(vec_clone_copy, vec_copy)
        vec_clone_state <- c(vec_clone_state, vec_state)
        vec_clone_Min <- c(vec_clone_Min, vec_Min)
        vec_clone_Maj <- c(vec_clone_Maj, vec_Maj)
    }
    #   Store the CN profile for the clone
    CN_profile <- data.frame(vec_clone_chr, vec_clone_start, vec_clone_end, vec_clone_copy, vec_clone_state, vec_clone_Min, vec_clone_Maj)
    names(CN_profile) <- c("chr", "start", "end", "copy", "state", "Min", "Maj")
    return(CN_profile)
}
