p3_internal_node_cn_profiles_wide <- function(simulation) {
    cn_profiles_wide <- simulation$sample$cn_profiles_wide

    sample_genotype_internal_nodes <- simulation$sample_phylogeny$package_cell_phylogeny_hclust_extra$hclust_internal_genotypes
    sample_internal_node_ID <- simulation$sample$sample_internal_node_ID

    sample_genotype_internal_nodes_unique <- simulation$sample$sample_genotype_internal_nodes_unique
    sample_genotype_internal_nodes_unique_profile <- simulation$sample$sample_genotype_internal_nodes_unique_profile
    #-------------------Find the CN profiles for each node in the sample
    for (node in 1:length(sample_genotype_internal_nodes)) {
        clone_ID <- sample_genotype_internal_nodes[node]
        i_clone <- which(sample_genotype_internal_nodes_unique == clone_ID)
        #   Find the CN profile for this node
        if (length(i_clone) == 0) {
            cell_genotype_profile <- sample_genotype_internal_nodes_unique_profile[[1]]
            cell_genotype_profile$copy <- NA
            cell_genotype_profile$state <- NA
            cell_genotype_profile$Min <- NA
            cell_genotype_profile$Maj <- NA
        } else {
            cell_genotype_profile <- sample_genotype_internal_nodes_unique_profile[[i_clone]]
        }
        #   Add column for node ID
        cell_id <- sample_internal_node_ID[node]
        cn_profiles_wide[cell_id] <- cell_genotype_profile$state
    }
    simulation$sample$cn_profiles_wide <- cn_profiles_wide
    return(simulation)
}
