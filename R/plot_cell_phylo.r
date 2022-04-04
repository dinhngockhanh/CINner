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

        cell_phylogeny_birthtime <- simulation$sample_phylogeny$package_cell_phylogeny$phylogeny_birthtime

        cell_phylogeny_deathtime <- simulation$sample_phylogeny$package_cell_phylogeny$phylogeny_deathtime

        # print(cell_phylogeny_birthtime)
        # print(cell_phylogeny_deathtime)

        # cell_phylogeny_labels <- simulation$sample_phylogeny$
        #---------------------------------------Plot cell phylogeny tree
        jpeg(paste(model, "_cell_phylo_", i, ".jpeg", sep = ""), width = width, height = height)





        # ape::plot.phylo(cell_phylogeny_phylo, show.tip.label = FALSE, root.edge = TRUE, use.edge.length = TRUE, direction = "downward")
        # # ape::plot.phylo(cell_phylogeny_phylo, direction = "downward", show.tip.label = FALSE, root.edge = TRUE)
        #
        # ape::axisPhylo()
        #
        # # ape::tiplabels()

        p <- ggtree(cell_phylogeny_phylo) + theme_tree2()
        # + geom_rootedge(rootedge = 3000)

        print(p)



        dev.off()
    }
}
