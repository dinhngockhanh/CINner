#' @export
p3_internal_node_cn_profiles_long <- function(simulation) {
    cn_profiles_long <- simulation$sample$cn_profiles_long

    sample_genotype_internal_nodes <- simulation$sample_phylogeny$package_cell_phylogeny_hclust_extra$hclust_internal_genotypes
    sample_internal_node_ID <- simulation$sample$sample_internal_node_ID

    sample_genotype_internal_nodes_unique <- simulation$sample$sample_genotype_internal_nodes_unique
    sample_genotype_internal_nodes_unique_profile <- simulation$sample$sample_genotype_internal_nodes_unique_profile
    #-------------------Find the CN profiles for each node in the sample
    cn_profiles_long_list <- vector("list", length = (length(sample_genotype_internal_nodes) + 1))
    cn_profiles_long_list[[1]] <- cn_profiles_long
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
        cell_genotype_profile$cell_id <- sample_internal_node_ID[node]
        names(cell_genotype_profile) <- c("chr", "start", "end", "copy", "state", "Min", "Maj", "cell_id")
        #   Add record for this node to list
        cn_profiles_long_list[[node + 1]] <- cell_genotype_profile
    }
    #   Bind all cells' CN profiles together
    cn_profiles_long <- rbindlist(cn_profiles_long_list, use.names = FALSE, fill = FALSE, idcol = NULL)
    class(cn_profiles_long) <- "data.frame"
    simulation$sample$cn_profiles_long <- cn_profiles_long
    return(simulation)
}
