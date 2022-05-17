# ===================================PLOT CELL EVOLUTION AS PHYLOGENY TREE
plot_cell_phylo <- function(model = "",
                            n_simulations = 0,
                            width = 1000,
                            height = 500) {
    for (i in 1:n_simulations) {
        #------------------------------------------Input simulation file
        filename <- paste(model, "_simulation_", i, ".rda", sep = "")
        load(filename)
        #--------------------------------------Input cell phylogeny tree
        cell_phylogeny_phylo <- simulation$sample_phylogeny$phylogeny_clustering_truth$tree
        cell_phylogeny_deathtime <- simulation$sample_phylogeny$package_cell_phylogeny$phylogeny_deathtime

        # print("===============================================================")
        # print(cell_phylogeny_birthtime)
        # print("===============================================================")
        # print(cell_phylogeny_deathtime)
        # print("===============================================================")











        #-------------------------------------Find length of clonal edge
        l_clonal_edge <- cell_phylogeny_deathtime[1]

        # cell_phylogeny_labels <- simulation$sample_phylogeny$
        #---------------------------------------Plot cell phylogeny tree
        jpeg(paste(model, "_sim", i, "_cell_phylo", ".jpeg", sep = ""), width = width, height = height)

        # ape::plot.phylo(cell_phylogeny_phylo, show.tip.label = FALSE, root.edge = TRUE, use.edge.length = TRUE, direction = "downward")
        # # ape::plot.phylo(cell_phylogeny_phylo, direction = "downward", show.tip.label = FALSE, root.edge = TRUE)
        #
        # ape::axisPhylo()
        #
        # # ape::tiplabels()




        p <- ggtree(cell_phylogeny_phylo) + geom_rootedge(rootedge = l_clonal_edge)

        # cell_phylogeny_hclust <- simulation$sample_phylogeny$cell_phylogeny_hclust
        # print("+++++++++++++++++++++++++++++++++++++++++++++++++++++++")
        # print(cell_phylogeny_hclust$merge)
        # print("+++++++++++++++++++++++++++++++++++++++++++++++++++++++")
        # print(cell_phylogeny_hclust$height)
        # print("+++++++++++++++++++++++++++++++++++++++++++++++++++++++")
        # p <- ggtree(cell_phylogeny_hclust)



        print(p)

        dev.off()
    }
}
