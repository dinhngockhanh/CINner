# ================================PLOT CLONAL EVOLUTION AS PHYLOGENY TREE
plot_clonal_phylo <- function(model = "",
                              n_simulations = 0,
                              width = 1000,
                              height = 500) {
    for (i in 1:n_simulations) {
        #------------------------------------------Input simulation file
        filename <- paste(model, "_simulation_", i, ".rda", sep = "")
        load(filename)
        #-------------------------------------Input clone phylogeny tree
        clone_phylogeny_phylo <- simulation$sample_phylogeny$clone_phylogeny_phylo


        print(clone_phylogeny_phylo)
        # print(clone_phylogeny_phylo$height)
        print(ape::node.height(clone_phylogeny_phylo))
        print(max(ape::node.height(clone_phylogeny_phylo)))
        print(80 - max(ape::node.height(clone_phylogeny_phylo)))

        clone_phylogeny_labels <- simulation$sample_phylogeny$package_clone_phylogeny$clone_phylogeny_labels
        #--------------------------------------Plot clone phylogeny tree
        if (length(clone_phylogeny_labels) > 1) {
            jpeg(paste(model, "_clonal_phylo_", i, ".jpeg", sep = ""), width = width, height = height)

            p <- ggtree(clone_phylogeny_phylo) + geom_tiplab() + theme_tree2() + geom_rootedge(rootedge = 80 - max(ape::node.height(clone_phylogeny_phylo)))
            # + layout_dendrogram()

            print(p)




            # # i <- ape::which.edge(clone_phylogeny_phylo, c("A"))
            # ape::plot.phylo(clone_phylogeny_phylo, direction = "downward", font = 3)
            # # ape::edgelabels(text='A',edge=i)
            # # ape::tiplabels()

            dev.off()
        } else {
            print("CANNOT PLOT CLONAL EVOLUTION FOR ONLY ONE CLONE")
        }
    }
}
