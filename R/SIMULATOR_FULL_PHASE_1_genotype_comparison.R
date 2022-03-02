#=================================COMPARE IF TWO GENOTYPES ARE IDENTICAL
SIMULATOR_FULL_PHASE_1_genotype_comparison <- function(genotype_1,genotype_2) {
#------------------------------------------Initialize flag of comparison
    output              <- 1
#---------------------------------------Compare chromosome strand counts
    ploidy_chrom_1      <- genotype_list_ploidy_chrom[[genotype_1]]
    ploidy_chrom_2      <- genotype_list_ploidy_chrom[[genotype_2]]
    if (any(ploidy_chrom_1!=ploidy_chrom_2)) {
        output          <- 0
        return(output)
    }
#-----------------------------Compare block CN on each chromosome strand
    ploidy_block_1      <- genotype_list_ploidy_block[[genotype_1]]
    ploidy_block_2      <- genotype_list_ploidy_block[[genotype_2]]
    for (chrom in 1:N_chromosomes) {
        chrom_ploidy    <- ploidy_chrom_1[chrom]
        for (strand in 1:chrom_ploidy) {

            # if (strand==0){
            #     next
            # }

            if (any(ploidy_block_1[[chrom]][[strand]] != ploidy_block_2[[chrom]][[strand]])) {
                output  <- 0
                return(output)
            }
        }
    }
#--------------------------------------Compare chromosome strand alleles
    ploidy_allele_1     <- genotype_list_ploidy_allele[[genotype_1]]
    ploidy_allele_2     <- genotype_list_ploidy_allele[[genotype_2]]
    for (chrom in 1:N_chromosomes) {
        chrom_ploidy    <- ploidy_chrom_1[chrom]
        for (strand in 1:chrom_ploidy) {
            if ((nrow(ploidy_allele_1[[chrom]][[strand]]) != nrow(ploidy_allele_2[[chrom]][[strand]])) || (ncol(ploidy_allele_1[[chrom]][[strand]]) != ncol(ploidy_allele_2[[chrom]][[strand]]))){
                output  <- 0
                return(output)
            }
            if (any(ploidy_allele_1[[chrom]][[strand]] != ploidy_allele_2[[chrom]][[strand]])) {
                output  <- 0
                return(output)
            }
        }
    }
#--------------------------------------------------Compare driver counts
    driver_count_1      <- genotype_list_driver_count[genotype_1]
    driver_count_2      <- genotype_list_driver_count[genotype_2]
    if (driver_count_1 != driver_count_2) {
        output          <- 0
        return(output)
    }
#----------------------------------------------------Compare driver maps
    driver_map_1        <- genotype_list_driver_map[[genotype_1]]
    driver_map_2        <- genotype_list_driver_map[[genotype_2]]
    if ((nrow(driver_map_1) != nrow(driver_map_2)) || (ncol(driver_map_1) != ncol(driver_map_2))){
        output          <- 0
        return(output)
    }
    if (any(driver_map_1 != driver_map_2)) {
        output          <- 0
        return(output)
    }
    return(output)
}
