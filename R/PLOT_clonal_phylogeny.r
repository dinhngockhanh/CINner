#================================PLOT CLONAL EVOLUTION AS PHYLOGENY TREE
PLOT_clonal_phylogeny <- function(package_simulation){
    package_sample_phylogeny        <- package_simulation[[3]]
    package_clone_phylogeny         <- package_sample_phylogeny[[4]]

    clone_phylogeny_phylo           <- package_sample_phylogeny[[2]]
    clone_phylogeny_labels          <- package_clone_phylogeny[[1]]

    if (length(clone_phylogeny_labels)>1){
        i       <- ape::which.edge(clone_phylogeny_phylo, c("A"))

        ape::plot.phylo(clone_phylogeny_phylo,direction="downward")

        # ape::edgelabels(text='A',edge=i)

        # ape::tiplabels()
    }else{
        print('CANNOT PLOT CLONAL EVOLUTION FOR ONLY ONE CLONE')
    }

}
