# ================================PLOT CLONAL EVOLUTION AS PHYLOGENY TREE
plot_clonal_phylo <- function(model = "",
                              n_simulations = 0,
                              width = 1000,
                              height = 500) {
    for (i in 1:n_simulations) {
        #------------------------------------------Input simulation file
        filename <- paste(model, "_simulation_", i, ".rda", sep = "")
        load(filename)
        #------------------------------------------Input clone evolution
        evolution_origin <- simulation$clonal_evolution$evolution_origin
        evolution_genotype_changes <- simulation$clonal_evolution$evolution_genotype_changes
        #-------------------------------------Input clone phylogeny tree
        clone_phylogeny_phylo <- simulation$sample_phylogeny$clone_phylogeny_phylo

        clone_phylogeny_genotype <- simulation$sample_phylogeny$package_clone_phylogeny$clone_phylogeny_genotype
        clone_phylogeny_origin <- simulation$sample_phylogeny$package_clone_phylogeny$clone_phylogeny_origin
        clone_phylogeny_deathtime <- simulation$sample_phylogeny$package_clone_phylogeny$clone_phylogeny_deathtime
        #-----------------------------------Find labels of clonal leaves
        # l_clonal_edge <- clone_phylogeny_deathtime[1]
        clone_phylogeny_labels <- simulation$sample_phylogeny$package_clone_phylogeny$clone_phylogeny_labels
        N_clones <- length(clone_phylogeny_labels)






        #--------------------------------------Plot clone phylogeny tree
        if (N_clones > 1) {
            #   Initiate destination file for plot
            jpeg(paste(model, "_clonal_phylo_", i, ".jpeg", sep = ""), width = width, height = height)
            #   Plot clone phylogeny tree
            p <- ggtree(clone_phylogeny_phylo, branch.length = "none")
            #   Plot clone ID
            size_text <- floor(height / N_clones)
            p <- p + geom_tiplab(as_ylab = TRUE, size = size_text)
            #   Plot phylogeny tree root edge
            p <- p + geom_rootedge(rootedge = 1)
            # p <- p + geom_rootedge(rootedge = l_clonal_edge)
            #-----------------------------Find node length in dendrogram
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
            print("~~~~~~~~~~~~~~~~~~~~~")
            print(clone_phylogeny_dend_depth)
            print("~~~~~~~~~~~~~~~~~~~~~")
            print(clone_phylogeny_dend_depth[N_clones:(2 * N_clones - 1)])







            #---------------------------------Plot events in each branch
            #   Choose color for each event type
            cols <- setNames(palette()[1:8], c(
                "Drivers", "WGD", "Missegregation", "Arm_Missegregation",
                "Amplification", "Deletion", "Interstitial_CNLOH", "Terminal_CNLOH"
            ))

            p <- p + scale_color_manual(values = cols) + theme(legend.position = "right")





            for (node_plot in 1:(2 * N_clones - 1)) {
                #   Genotype at the end of edge
                if (node_plot <= N_clones) {
                    node_phylo <- N_clones - 1 + node_plot
                } else {
                    node_phylo <- node_plot - N_clones
                }


                dat <- data.frame(
                    Drivers = 1,
                    WGD = 0,
                    Missegregation = 0,
                    Arm_Missegregation = 0,
                    Amplification = 0,
                    Deletion = 0,
                    Interstitial_CNLOH = 0,
                    Terminal_CNLOH = 0
                )
                dat$node <- node_plot
                event_pie <- nodepie(dat, cols = 1:8)
                event_pie <- lapply(
                    event_pie,
                    function(g) g + scale_fill_manual(values = cols)
                )



                p <- p + geom_inset(
                    event_pie,
                    width = 0.08, height = 0.08, hjust = 0, x = "branch"
                )


                # p <- p + geom_text(
                #
                # )

                p <- p + geom_label(
                    mapping = NULL,
                    data = NULL
                )
            }
            # # ===========================================================
            # # ===========================================================
            # # ===========================================================
            #
            # for (node_plot in 1:(2 * N_clones - 1)) {
            #     #   Genotype at the end of edge
            #     if (node_plot <= N_clones) {
            #         node_phylo <- N_clones - 1 + node_plot
            #     } else {
            #         node_phylo <- node_plot - N_clones
            #     }
            #     node_genotype <- clone_phylogeny_genotype[node_phylo]
            #     #   Genotype at the beginning of edge
            #     node_mother_phylo <- clone_phylogeny_origin[node_phylo]
            #     if (node_mother_phylo == 0) {
            #         node_mother_genotype <- 0
            #     } else {
            #         node_mother_genotype <- clone_phylogeny_genotype[node_mother_phylo]
            #     }
            #     #   Move on if no event happened within edge
            #     if (node_genotype == node_mother_genotype) {
            #         next
            #     }
            #     #   List of events to plot
            #     n_all_events <- 0
            #     list_all_events <- c()
            #     list_all_notations <- c()
            #     genotype_current <- node_genotype
            #     while (genotype_current != node_mother_genotype) {
            #         n_events <- length(evolution_genotype_changes[[genotype_current]])
            #         if (n_events == 0) {
            #             genotype_current <- evolution_origin[genotype_current]
            #             next
            #         }
            #         for (event in 1:n_events) {
            #             n_all_events <- n_all_events + 1
            #             list_all_events[n_all_events] <- evolution_genotype_changes[[genotype_current]][[event]][1]
            #             if (list_all_events[n_all_events] != "whole-genome-duplication") {
            #                 list_all_notations[n_all_events] <- evolution_genotype_changes[[genotype_current]][[event]][2]
            #             } else {
            #                 list_all_notations[n_all_events] <- ""
            #             }
            #         }
            #         genotype_current <- evolution_origin[genotype_current]
            #     }
            #     print("===============================================")
            #     print(n_all_events)
            #     print("-----------------------------------------------")
            #     print(list_all_events)
            #     print("-----------------------------------------------")
            #     print(list_all_notations)
            #     #   Decide length between two events
            # }

            # ===========================================================
            # ===========================================================
            # ===========================================================
            # #   Plot events for each edge
            # for (node_plot in 1:(2 * N_clones - 1)) {
            #     #   Genotype at the end of edge
            #     if (node_plot <= N_clones) {
            #         node_phylo <- N_clones - 1 + node_plot
            #     } else {
            #         node_phylo <- node_plot - N_clones
            #     }
            #     node_genotype <- clone_phylogeny_genotype[node_phylo]
            #     #   Genotype at the beginning of edge
            #     node_mother_phylo <- clone_phylogeny_origin[node_phylo]
            #     if (node_mother_phylo == 0) {
            #         node_mother_genotype <- 0
            #     } else {
            #         node_mother_genotype <- clone_phylogeny_genotype[node_mother_phylo]
            #     }
            #     #   Move on if no event happened within edge
            #     if (node_genotype == node_mother_genotype) {
            #         next
            #     }
            #     #   Count number of each event types within edge
            #     n_drivers <- 0
            #     n_WGD <- 0
            #     n_missegregation <- 0
            #     n_chrom_arm_missegregation <- 0
            #     n_amplification <- 0
            #     n_deletion <- 0
            #     n_cnloh_interstitial <- 0
            #     n_cnloh_terminal <- 0
            #     genotype_current <- node_genotype
            #     while (genotype_current != node_mother_genotype) {
            #         if (length(evolution_genotype_changes[[genotype_current]]) == 0) {
            #             genotype_current <- evolution_origin[genotype_current]
            #             next
            #         }
            #         for (event in 1:length(evolution_genotype_changes[[genotype_current]])) {
            #             event_type <- evolution_genotype_changes[[genotype_current]][[event]][1]
            #             if (event_type == "new-driver") {
            #                 n_drivers <- n_drivers + 1
            #             } else {
            #                 if (event_type == "whole-genome-duplication") {
            #                     n_WGD <- n_WGD + 1
            #                 } else {
            #                     if (event_type == "missegregation") {
            #                         n_missegregation <- n_missegregation + 1
            #                     } else {
            #                         if (event_type == "chromosome-arm-missegregation") {
            #                             n_chrom_arm_missegregation <- n_chrom_arm_missegregation + 1
            #                         } else {
            #                             if (event_type == "focal-amplification") {
            #                                 n_amplification <- n_amplification + 1
            #                             } else {
            #                                 if (event_type == "focal-deletion") {
            #                                     n_deletion <- n_deletion + 1
            #                                 } else {
            #                                     if (event_type == "cnloh-interstitial") {
            #                                         n_cnloh_interstitial <- n_cnloh_interstitial + 1
            #                                     } else {
            #                                         if (event_type == "cnloh-terminal") {
            #                                             n_cnloh_terminal <- n_cnloh_terminal + 1
            #                                         }
            #                                     }
            #                                 }
            #                             }
            #                         }
            #                     }
            #                 }
            #             }
            #         }
            #         genotype_current <- evolution_origin[genotype_current]
            #     }
            #     #   Find total count of events
            #     n_total_events <- n_drivers + n_WGD + n_missegregation + n_chrom_arm_missegregation
            #     +n_amplification + n_deletion
            #     +n_cnloh_interstitial + n_cnloh_terminal
            #     #   Skip this edge if there is no event
            #     if (n_total_events == 0) {
            #         next
            #     }
            #     #   Normalize number of each event types
            #     n_drivers <- n_drivers / n_total_events
            #     n_WGD <- n_WGD / n_total_events
            #     n_missegregation <- n_missegregation / n_total_events
            #     n_chrom_arm_missegregation <- n_chrom_arm_missegregation / n_total_events
            #     n_amplification <- n_amplification / n_total_events
            #     n_deletion <- n_deletion / n_total_events
            #     n_cnloh_interstitial <- n_cnloh_interstitial / n_total_events
            #     n_cnloh_terminal <- n_cnloh_terminal / n_total_events
            #     #   Add plot of event counts for the edge
            #     dat <- data.frame(
            #         Drivers = n_drivers,
            #         WGD = n_WGD,
            #         Missegregation = n_missegregation,
            #         Arm_Missegregation = n_chrom_arm_missegregation,
            #         Amplification = n_amplification,
            #         Deletion = n_deletion,
            #         Interstitial_CNLOH = n_cnloh_interstitial,
            #         Terminal_CNLOH = n_cnloh_terminal
            #     )
            #     dat$node <- node_plot
            #     event_pie <- nodepie(dat, cols = 1:8)
            #     event_pie <- lapply(
            #         event_pie,
            #         function(g) g + scale_fill_manual(values = cols)
            #     )
            #
            #
            #     p <- p + geom_inset(
            #         event_pie,
            #         width = 0.1, height = 0.1, hjust = 0, x = "branch"
            #     )
            # }
            # ===========================================================
            # ===========================================================
            # ===========================================================

            # p <- p + scale_color_manual(values = cols)





            # p <- p + scale_color_manual(values = cols) + theme(legend.position = "right")
            # print(cols)

            # p <- p + theme_bw()



            # # + layout_dendrogram()


            print(p)

            dev.off()
        } else {
            print("CANNOT PLOT CLONAL EVOLUTION FOR ONLY ONE CLONE")
        }
    }
}
