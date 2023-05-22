#' @export
p3_cn_profiles_internal <- function(simulation) {
    sample_genotype_internal_nodes <- simulation$sample_phylogeny$package_cell_phylogeny_hclust_extra$hclust_internal_genotypes
    clonal_evolution <- simulation$clonal_evolution
    #------------------------Find the CN profiles for each internal node
    sample_genotype_internal_nodes_unique <- unique(sample_genotype_internal_nodes)
    sample_genotype_internal_nodes_unique <- sample_genotype_internal_nodes_unique[which(sample_genotype_internal_nodes_unique > 0)]
    sample_genotype_internal_nodes_unique_profile <- list()
    for (i_node in 1:length(sample_genotype_internal_nodes_unique)) {
        #   Extract CN information for the node from clonal evolution data
        clone_ID <- sample_genotype_internal_nodes_unique[i_node]
        genotype_unique_profile <- get_cn_profile(clonal_evolution, clone_ID)
        sample_genotype_internal_nodes_unique_profile[[i_node]] <- genotype_unique_profile
    }
    #---------------------------------------Find ID's for internal nodes
    sample_internal_node_ID <- paste("Internal-node-", 1:length(sample_genotype_internal_nodes), sep = "")
    #-----------------------------------------------------Output results
    simulation$sample$sample_internal_node_ID <- sample_internal_node_ID
    simulation$sample$sample_genotype_internal_nodes_unique <- sample_genotype_internal_nodes_unique
    simulation$sample$sample_genotype_internal_nodes_unique_profile <- sample_genotype_internal_nodes_unique_profile
    return(simulation)
}
