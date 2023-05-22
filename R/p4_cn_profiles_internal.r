#' @export
p4_cn_profiles_internal <- function(simulation) {
    hclust_nodes <- simulation$sample_phylogeny$package_cell_phylogeny_hclust_extra$hclust_nodes
    hclust_internal_genotypes <- simulation$sample_phylogeny$package_cell_phylogeny_hclust_extra$hclust_internal_genotypes
    neuvar_phylogeny_genotype <- simulation$neutral_variations$cell_phylogeny$neuvar_phylogeny_genotype
    clonal_evolution <- simulation$neutral_variations$clonal_evolution
    #-----------------------------------Find the internal node clonal ID
    hclust_internal_genotypes_neuvar <- rep(0, length(hclust_internal_genotypes))
    for (i in 1:length(hclust_internal_genotypes)) {
        if (length(which(hclust_nodes == i)) == 0) {
            #   Artificial nodes for merging unmerged nodes at time 0
            hclust_internal_genotypes_neuvar[i] <- 0
        } else {
            #   Normal nodes with defined clonal CN
            hclust_internal_genotypes_neuvar[i] <- neuvar_phylogeny_genotype[which(hclust_nodes == i)]
        }
    }
    #------------------------Find the CN profiles for each internal node
    sample_genotype_internal_nodes_unique_neuvar <- unique(hclust_internal_genotypes_neuvar)
    sample_genotype_internal_nodes_unique_neuvar <- sample_genotype_internal_nodes_unique_neuvar[which(sample_genotype_internal_nodes_unique_neuvar > 0)]
    sample_genotype_internal_nodes_unique_profile_neuvar <- list()
    for (i_node in 1:length(sample_genotype_internal_nodes_unique_neuvar)) {
        #   Extract CN information for the node from clonal evolution data
        clone_ID <- sample_genotype_internal_nodes_unique_neuvar[i_node]
        genotype_unique_profile <- get_cn_profile(clonal_evolution, clone_ID)
        sample_genotype_internal_nodes_unique_profile_neuvar[[i_node]] <- genotype_unique_profile
    }
    #-----------------------------------------------------Output results
    simulation$sample$hclust_internal_genotypes_neuvar <- hclust_internal_genotypes_neuvar
    simulation$sample$sample_genotype_internal_nodes_unique_neuvar <- sample_genotype_internal_nodes_unique_neuvar
    simulation$sample$sample_genotype_internal_nodes_unique_profile_neuvar <- sample_genotype_internal_nodes_unique_profile_neuvar
    return(simulation)
}
