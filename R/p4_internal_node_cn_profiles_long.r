#' @export
p4_internal_node_cn_profiles_long <- function(simulation) {
    cn_profiles_long <- simulation$neutral_variations$sample$cn_profiles_long

    hclust_internal_genotypes_neuvar <- simulation$sample$hclust_internal_genotypes_neuvar
    sample_internal_node_ID <- simulation$sample$sample_internal_node_ID

    sample_genotype_internal_nodes_unique_neuvar <- simulation$sample$sample_genotype_internal_nodes_unique_neuvar
    sample_genotype_internal_nodes_unique_profile_neuvar <- simulation$sample$sample_genotype_internal_nodes_unique_profile_neuvar
    #-------------------Find the CN profiles for each node in the sample
    cn_profiles_long_list <- vector("list", length = (length(hclust_internal_genotypes_neuvar) + 1))
    cn_profiles_long_list[[1]] <- cn_profiles_long
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
        cell_genotype_profile$cell_id <- sample_internal_node_ID[node]
        names(cell_genotype_profile) <- c("chr", "start", "end", "copy", "state", "Min", "Maj", "cell_id")
        #   Add record for this node to list
        cn_profiles_long_list[[node + 1]] <- cell_genotype_profile
    }
    #   Bind all cells' CN profiles together
    cn_profiles_long <- rbindlist(cn_profiles_long_list, use.names = FALSE, fill = FALSE, idcol = NULL)
    class(cn_profiles_long) <- "data.frame"
    simulation$neutral_variations$sample$cn_profiles_long <- cn_profiles_long
    return(simulation)
}
