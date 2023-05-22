#' @export
p4_internal_node_cn_profiles_wide <- function(simulation) {
    cn_profiles_wide <- simulation$neutral_variations$sample$cn_profiles_wide

    hclust_internal_genotypes_neuvar <- simulation$sample$hclust_internal_genotypes_neuvar
    sample_internal_node_ID <- simulation$sample$sample_internal_node_ID

    sample_genotype_internal_nodes_unique_neuvar <- simulation$sample$sample_genotype_internal_nodes_unique_neuvar
    sample_genotype_internal_nodes_unique_profile_neuvar <- simulation$sample$sample_genotype_internal_nodes_unique_profile_neuvar
    #-------------------Find the CN profiles for each node in the sample
    for (node in 1:length(hclust_internal_genotypes_neuvar)) {
        clone_ID <- hclust_internal_genotypes_neuvar[node]
        i_clone <- which(sample_genotype_internal_nodes_unique_neuvar == clone_ID)
        #   Find the CN profile for this node
        if (length(i_clone) == 0) {
            cell_genotype_profile <- sample_genotype_internal_nodes_unique_profile_neuvar[[1]]
            cell_genotype_profile$copy <- NA
            cell_genotype_profile$state <- NA
            cell_genotype_profile$Min <- NA
            cell_genotype_profile$Maj <- NA
        } else {
            cell_genotype_profile <- sample_genotype_internal_nodes_unique_profile_neuvar[[i_clone]]
        }
        #   Add column for node ID
        cell_id <- sample_internal_node_ID[node]
        cn_profiles_wide[cell_id] <- cell_genotype_profile$state
    }
    simulation$neutral_variations$sample$cn_profiles_wide <- cn_profiles_wide
    return(simulation)
}
