#==================================PLOT HEATMAP OF CN PROFILES IN SAMPLE
PLOT_copy_number_heatmap <- function(package_simulation,phylo,plotcol,reorderclusters,plottree){
#------------------------Find phylogeny and clustering for sampled cells
    if (phylo=='true'){
        package_sample                <- package_simulation[[2]]
        package_sample_phylogeny      <- package_simulation[[3]]

        sample_genotype_profiles      <- package_sample[[1]]
        phylogeny_clustering_truth    <- package_sample_phylogeny[[1]]
    }
#--------------------------------------------Plot heatmap of CN profiles
    if (plotcol=='total_copy'){
        # plotHeatmap(sample_genotype_profiles,plotcol="state",tree=phylogeny_clustering_truth$tree,reorderclusters=TRUE,plottree=TRUE)

        plotHeatmap(sample_genotype_profiles,plotcol="state",clusters=phylogeny_clustering_truth$clustering,tree=phylogeny_clustering_truth$tree,reorderclusters=TRUE,plottree=TRUE)
    }else{if (plotcol=='min_copy'){
        plotHeatmap(sample_genotype_profiles,plotcol="Min",clusters=phylogeny_clustering_truth$clustering,tree=phylogeny_clustering_truth$tree,reorderclusters=TRUE,plottree=TRUE)
    }}
}
