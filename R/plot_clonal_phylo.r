# ================================PLOT CLONAL EVOLUTION AS PHYLOGENY TREE
#' @export
plot_clonal_phylo <- function(model = "",
                              n_simulations = 0,
                              folder_workplace = "",
                              width = 1000,
                              height = 500,
                              compute_parallel = TRUE,
                              n_cores = NULL) {
    library(ggtree)
    library(adephylo)
    library(ggplot2)
    if (compute_parallel == FALSE) {
        #-----------------------Plot clonal evolution in sequential mode
        for (iteration in 1:n_simulations) {
            plot_clonal_phylo_one_simulation(
                model,
                iteration,
                folder_workplace,
                width,
                height
            )
        }
    } else {
        #-------------------------Plot clonal evolution in parallel mode
        library(pbapply)
        #   Start parallel cluster
        if (is.null(n_cores)) {
            numCores <- detectCores()
        } else {
            numCores <- n_cores
        }
        cl <- makePSOCKcluster(numCores - 1)
        #   Prepare input parameters for plotting
        model <<- model
        folder_workplace <<- folder_workplace
        width <<- width
        height <<- height
        plot_clonal_phylo_one_simulation <<- plot_clonal_phylo_one_simulation
        clusterExport(cl, varlist = c(
            "plot_clonal_phylo_one_simulation",
            "model",
            "folder_workplace",
            "width",
            "height"
        ))
        clusterEvalQ(cl = cl, require(ggtree))
        clusterEvalQ(cl = cl, require(adephylo))
        clusterEvalQ(cl = cl, require(ggplot2))
        #   Plot in parallel
        pblapply(cl = cl, X = 1:n_simulations, FUN = function(iteration) {
            plot_clonal_phylo_one_simulation(
                model,
                iteration,
                folder_workplace,
                width,
                height
            )
        })
        #   Stop parallel cluster
        stopCluster(cl)
    }
}

plot_clonal_phylo_one_simulation <- function(model,
                                             iteration,
                                             folder_workplace,
                                             width,
                                             height) {
    #----------------------------------------------Input simulation file
    filename <- paste(folder_workplace, model, "_simulation_", iteration, ".rda", sep = "")
    load(filename)
    #----------------------------------------------Input clone evolution
    evolution_origin <- simulation$clonal_evolution$evolution_origin
    evolution_genotype_changes <- simulation$clonal_evolution$evolution_genotype_changes
    #-----------------------------------------Input clone phylogeny tree
    clone_phylogeny_phylo <- simulation$sample_phylogeny$clone_phylogeny_phylo

    clone_phylogeny_genotype <- simulation$sample_phylogeny$package_clone_phylogeny$clone_phylogeny_genotype
    clone_phylogeny_origin <- simulation$sample_phylogeny$package_clone_phylogeny$clone_phylogeny_origin
    #---------------------------------------Find labels of clonal leaves
    clone_phylogeny_labels <- simulation$sample_phylogeny$package_clone_phylogeny$clone_phylogeny_labels
    N_clones <- length(clone_phylogeny_labels)
    if (N_clones <= 1) {
        print("CURRENTLY CANNOT PLOT CLONAL EVOLUTION FOR ONLY ONE CLONE")
        return()
    }
    #-----------------------------------Choose color for each event type
    cols <- c(
        "Drivers" = "blue",
        "WGD" = "red",
        "Missegregation" = "green",
        "Arm_Missegregation" = "darkgoldenrod1",
        "Amplification" = "cadetblue1",
        "Deletion" = "burlywood4",
        "Interstitial_CNLOH" = "blueviolet",
        "Terminal_CNLOH" = "darkgray"
    )
    #--------------------------------------------Identify nodes in phylo
    #---Find list of leaves for each internal node in phylo
    list_leaves_phylo <- listTips(clone_phylogeny_phylo)
    for (node in 1:length(list_leaves_phylo)) {
        list_leaves_phylo[[node]] <- sort(list_leaves_phylo[[node]])
    }
    #---Find list of leaves for each internal node in our record
    list_leaves_record <- vector("list", length = (N_clones - 1))
    for (node_leaf in N_clones:(2 * N_clones - 1)) {
        # id_leaf <- clone_phylogeny_labels[node_leaf - N_clones + 1]
        node_index <- node_leaf - N_clones + 1
        node <- node_leaf
        while (node > 0) {
            if (node <= (N_clones - 1)) {
                list_leaves_record[[node]] <- c(list_leaves_record[[node]], node_index)
            }
            node <- clone_phylogeny_origin[node]
        }
    }
    for (node in 1:length(list_leaves_record)) {
        list_leaves_record[[node]] <- sort(list_leaves_record[[node]])
    }
    #---Find index in our record for every node in phylo
    vec_node_id <- rep(0, length = (2 * N_clones - 1))
    #   Assign indices for leaf nodes
    vec_node_id[1:N_clones] <- N_clones:(2 * N_clones - 1)
    #   Assign indices for internal nodes
    for (node_phylo in (N_clones + 1):(2 * N_clones - 1)) {
        node_internal <- node_phylo - N_clones
        leaves <- list_leaves_phylo[[node_internal]]
        node_record <- 0
        for (node in 1:length(list_leaves_record)) {
            if (length(list_leaves_record[[node]]) != length(leaves)) {
                next
            }
            if (all(leaves == list_leaves_record[[node]])) {
                node_record <- node
                break
            }
        }
        vec_node_id[node_phylo] <- node_record
    }
    #-------------------------------------Find length for tree root edge
    for (node_plot in 1:(2 * N_clones - 1)) {
        #   Genotype at the end of edge
        node_phylo <- vec_node_id[node_plot]
        node_genotype <- clone_phylogeny_genotype[node_phylo]
        #   Genotype at the beginning of edge
        node_mother_phylo <- clone_phylogeny_origin[node_phylo]
        if (node_mother_phylo == 0) {
            node_root <- node_phylo
            node_mother_genotype <- 0
            #   Count events happening on tree root edge
            L_root <- 0
            genotype_current <- node_genotype
            while ((genotype_current != node_mother_genotype) & (genotype_current > 0)) {
                L_root <- L_root + 0.5 * length(evolution_genotype_changes[[genotype_current]])
                genotype_current <- evolution_origin[genotype_current]
            }
        }
    }
    #------------------------------------------Plot clone phylogeny tree
    #   Initiate destination file for plot
    jpeg(paste(model, "_sim", iteration, "_clonal_phylo", ".jpeg", sep = ""), width = width, height = height)
    #   Plot clone phylogeny tree
    p <- ggtree(clone_phylogeny_phylo, branch.length = "none")
    #   Plot clone ID
    p <- p + geom_tiplab(as_ylab = TRUE, size = 100)
    #   Plot phylogeny tree root edge
    p <- p + geom_rootedge(rootedge = L_root)
    #-------------------------------------Find node length in dendrogram
    #   Find node depth in dendrogram
    clone_phylogeny_dend_depth <- rep(0, length(clone_phylogeny_origin))
    for (node_inner in (N_clones - 1):1) {
        nodes_daughter <- which(clone_phylogeny_origin == node_inner)
        if (length(nodes_daughter) == 0) {
            next
        }
        depth <- max(clone_phylogeny_dend_depth[nodes_daughter]) + 1
        clone_phylogeny_dend_depth[c(node_inner, nodes_daughter)] <- depth
    }
    #   Find node length in dendrogram
    clone_phylogeny_dend_length <- rep(0, length(clone_phylogeny_origin))
    clone_phylogeny_dend_length[N_clones:(2 * N_clones - 1)] <- clone_phylogeny_dend_depth[N_clones:(2 * N_clones - 1)]
    for (node in 1:(2 * N_clones - 1)) {
        nodes_daughter <- which(clone_phylogeny_origin == node)
        if (length(nodes_daughter) == 0) {
            next
        }
        clone_phylogeny_dend_length[node] <- clone_phylogeny_dend_depth[node] - clone_phylogeny_dend_depth[nodes_daughter][1]
    }
    clone_phylogeny_dend_length[node_root] <- L_root
    #-----------------------------------------Plot events in each branch
    #   Plot events in each branch
    for (node_plot in 1:(2 * N_clones - 1)) {
        #   Genotype at the end of edge
        node_phylo <- vec_node_id[node_plot]
        node_genotype <- clone_phylogeny_genotype[node_phylo]
        #   Genotype at the beginning of edge
        node_mother_phylo <- clone_phylogeny_origin[node_phylo]
        if (node_mother_phylo == 0) {
            node_mother_genotype <- 0
        } else {
            node_mother_genotype <- clone_phylogeny_genotype[node_mother_phylo]
        }
        #   List of events to plot
        all_event_types <- c()
        all_event_text <- c()
        genotype_current <- node_genotype
        while ((genotype_current != node_mother_genotype) & (genotype_current > 0)) {
            n_events <- length(evolution_genotype_changes[[genotype_current]])
            if (n_events == 0) {
                genotype_current <- evolution_origin[genotype_current]
                next
            }
            for (event in 1:n_events) {
                all_event_types <- c(all_event_types, evolution_genotype_changes[[genotype_current]][[event]][1])

                # print(evolution_genotype_changes[[genotype_current]][[event]])

                # all_event_text <- ??????????????????????????????
            }
            genotype_current <- evolution_origin[genotype_current]
        }
        #   Decide distances between events
        n_all_events <- length(all_event_types)
        dend_length <- clone_phylogeny_dend_length[node_phylo]
        if (node_mother_phylo == 0) {
            hjust_start <- L_root - 0.5
        } else {
            hjust_start <- dend_length / 2 - dend_length / (n_all_events + 1)
        }
        hjust_unit <- dend_length / (n_all_events + 1)
        #   Move on if no event happened within edge
        if (n_all_events == 0) {
            next
        }
        #   Plot each event
        for (event in 1:n_all_events) {
            #   Find which event type
            event_type <- all_event_types[event]
            Drivers <- 0
            WGD <- 0
            Missegregation <- 0
            Arm_Missegregation <- 0
            Amplification <- 0
            Deletion <- 0
            Interstitial_CNLOH <- 0
            Terminal_CNLOH <- 0
            if (event_type == "new-driver") {
                Drivers <- 1
            }
            if (event_type == "whole-genome-duplication") {
                WGD <- 1
            }
            if (event_type == "missegregation") {
                Missegregation <- 1
            }
            if (event_type == "chromosome-arm-missegregation") {
                Arm_Missegregation <- 1
            }
            if (event_type == "focal-amplification") {
                Amplification <- 1
            }
            if (event_type == "focal-deletion") {
                Deletion <- 1
            }
            if (event_type == "cnloh-interstitial") {
                Interstitial_CNLOH <- 1
            }
            if (event_type == "cnloh-terminal") {
                Terminal_CNLOH <- 1
            }
            #   Create pie chart for event accordingly
            dat <- data.frame(
                Drivers = Drivers,
                WGD = WGD,
                Missegregation = Missegregation,
                Arm_Missegregation = Arm_Missegregation,
                Amplification = Amplification,
                Deletion = Deletion,
                Interstitial_CNLOH = Interstitial_CNLOH,
                Terminal_CNLOH = Terminal_CNLOH
            )
            dat$node <- node_plot
            event_pie <- nodepie(dat, cols = 1:8)
            event_pie <- lapply(
                event_pie,
                function(g) g + scale_fill_manual(values = cols)
            )
            #   Attach pie chart onto phylogeny plot
            size_pie <- max(0.01, min(0.05, 0.05 * 10 / (N_clones)))
            p <- p + geom_inset(
                event_pie,
                width = size_pie, height = size_pie, hjust = (hjust_start - (event - 1) * hjust_unit), x = "branch"
            )
        }
    }
    #--------------------------------------------------Add in the legend
    x.range <- ggplot_build(p)$layout$panel_params[[1]]$x.range
    y.range <- ggplot_build(p)$layout$panel_params[[1]]$y.range
    for (row in 1:length(cols)) {
        p <- p + annotate("point", x.range[1], (y.range[2] - row * (y.range[2] - y.range[1]) / 20), size = 16, color = cols[row])
        p <- p + annotate("text", x.range[1] + (x.range[2] - x.range[1]) / 20, (y.range[2] - row * (y.range[2] - y.range[1]) / 20), size = 16, hjust = 0, label = names(cols)[row])
    }
    print(p)
    dev.off()
}
