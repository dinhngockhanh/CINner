#==================================PLOT HEATMAP OF CN PROFILES IN SAMPLE
PLOT_copy_number_heatmap <- function(package_simulation,phylo,plotcol,reorderclusters,plottree){
#------------------------Find phylogeny and clustering for sampled cells
    if (phylo=='true'){
        package_sample                <- package_simulation[[2]]
        sample_genotype_profiles      <- package_sample[[1]]
        package_sample_phylogeny      <- package_simulation[[3]]
        phylogeny_clustering_truth    <- package_sample_phylogeny[[1]]
    }
#--------------------------------------------Plot heatmap of CN profiles

sample_genotype_profiles
phylogeny_clustering_truth

    plotHeatmap(sample_genotype_profiles,plotcol,clusters=phylogeny_clustering_truth$clustering,tree=phylogeny_clustering_truth$tree,reorderclusters,plottree)
}
