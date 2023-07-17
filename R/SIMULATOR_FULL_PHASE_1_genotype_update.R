# =================UPDATE DNA LENGTH AND SELECTION RATES OF NEW GENOTYPES
#' @export
SIMULATOR_FULL_PHASE_1_genotype_update <- function(genotype, time) {
    #-----------------------------Get the CN profile of the new genotype
    ploidy_chrom <- genotype_list_ploidy_chrom[[genotype]]
    ploidy_allele <- genotype_list_ploidy_allele[[genotype]]
    ploidy_block <- genotype_list_ploidy_block[[genotype]]
    WGD_count <- genotype_list_WGD_count[genotype]
    driver_count <- genotype_list_driver_count[genotype]
    driver_map <- genotype_list_driver_map[[genotype]]
    #------------------------Update the DNA length for the new genotypes
    #   Compute the DNA length for genotype 1
    DNA_length <- 0
    for (chrom in 1:N_chromosomes) {
        no_strands <- ploidy_chrom[chrom]
        if (no_strands < 1) {
            next
        }
        for (strand in 1:no_strands) {
            DNA_length <- DNA_length + sum(ploidy_block[[chrom]][[strand]])
        }
    }
    DNA_length <- size_CN_block_DNA * DNA_length
    genotype_list_DNA_length[[genotype]] <<- DNA_length
    #--------------------Update the selection rate for the new genotypes
    genotype_list_selection_rate[genotype] <<- SIMULATOR_FULL_PHASE_1_selection_rate(
        WGD_count = WGD_count,
        driver_count = driver_count,
        driver_map = driver_map,
        ploidy_chrom = ploidy_chrom,
        ploidy_block = ploidy_block,
        ploidy_allele = ploidy_allele
    )
    #--------------------------Update prob(driver) for the new genotypes
    genotype_list_prob_new_drivers[genotype] <<- 1 - dpois(x = 0, lambda = rate_driver * DNA_length)
    #----------------------------Update prob(CNAs) for the new genotypes
    #   Find ingredients to compute probabilities of CNA/driver mutation
    chrom_ploidy <- genotype_list_ploidy_chrom[[genotype]]
    DNA_length <- genotype_list_DNA_length[[genotype]]
    WGD_count <- genotype_list_WGD_count[genotype]
    #   Find probability of new genotype
    if (mode_CN_WGD == "per_division") {
        prob_CN_WGD <- min(1, eval(parse(text = sub(".*:", "", formula_CN_whole_genome_duplication))))
    }
    if (mode_CN_misseg == "per_division") {
        prob_CN_misseg <- min(1, eval(parse(text = sub(".*:", "", formula_CN_missegregation))))
    } else if (mode_CN_misseg == "per_homolog") {
        prob_CN_misseg_homolog <- min(1, eval(parse(text = sub(".*:", "", formula_CN_missegregation))))
        genotype_list_prob_CN_misseg_homolog[genotype] <<- prob_CN_misseg_homolog
        prob_CN_misseg <- 1 - (1 - prob_CN_misseg_homolog)^sum(chrom_ploidy)
    }
    if (mode_CN_arm_misseg == "per_division") {
        prob_CN_arm_misseg <- min(1, eval(parse(text = sub(".*:", "", formula_CN_chrom_arm_missegregation))))
    } else if (mode_CN_arm_misseg == "per_homolog") {
        prob_CN_arm_misseg_homolog <- min(1, eval(parse(text = sub(".*:", "", formula_CN_chrom_arm_missegregation))))
        genotype_list_prob_CN_arm_misseg_homolog[genotype] <<- prob_CN_arm_misseg_homolog
        prob_CN_arm_misseg <- 1 - (1 - prob_CN_arm_misseg_homolog)^sum(chrom_ploidy)
    }
    if (mode_CN_foc_amp == "per_division") {
        prob_CN_foc_amp <- min(1, eval(parse(text = sub(".*:", "", formula_CN_focal_amplification))))
    }
    if (mode_CN_foc_del == "per_division") {
        prob_CN_foc_del <- min(1, eval(parse(text = sub(".*:", "", formula_CN_focal_deletion))))
    }
    if (mode_CN_cnloh_i == "per_division") {
        prob_CN_cnloh_i <- min(1, eval(parse(text = sub(".*:", "", formula_CN_cnloh_interstitial))))
    }
    if (mode_CN_cnloh_t == "per_division") {
        prob_CN_cnloh_t <- min(1, eval(parse(text = sub(".*:", "", formula_CN_cnloh_terminal))))
    }
    #   Update prob(CNAs) for the new genotypes
    prob_CNAs <- c(
        prob_CN_WGD,
        prob_CN_misseg, prob_CN_arm_misseg,
        prob_CN_foc_amp, prob_CN_foc_del,
        prob_CN_cnloh_i, prob_CN_cnloh_t
    )
    genotype_list_prob_CNAs[[genotype]] <<- prob_CNAs
}
