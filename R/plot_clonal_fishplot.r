# =====================================PLOT CLONAL EVOLUTION AS FISH PLOT
plot_clonal_fishplot <- function(model = "",
                                 n_simulations = 0,
                                 vec_time = c(0),
                                 unit_time = "year",
                                 width = 1000,
                                 height = 500) {
    for (i in 1:n_simulations) {
        #------------------------------------------Input simulation file
        filename <- paste(model, "_simulation_", i, ".rda", sep = "")
        load(filename)
        #-----------------------------Transform time points if necessary
        if (unit_time == "day") {
            vec_time_simulation <- vec_time
        } else {
            if (unit_time == "week") {
                vec_time_simulation <- 7 * vec_time
            } else {
                if (unit_time == "month") {
                    vec_time_simulation <- 30 * vec_time
                } else {
                    if (unit_time == "year") {
                        vec_time_simulation <- 365 * vec_time
                    }
                }
            }
        }
        #---------------------------------------Input the sampling table
        Table_sampling <- simulation$sample$Table_sampling
        #-------------------------------------Input the clonal evolution
        evolution_origin <- simulation$clonal_evolution$evolution_origin
        evolution_traj_time <- simulation$clonal_evolution$evolution_traj_time
        evolution_traj_clonal_ID <- simulation$clonal_evolution$evolution_traj_clonal_ID
        evolution_traj_population <- simulation$clonal_evolution$evolution_traj_population
        #-------------------------------------Input the clonal phylogeny
        clone_phylogeny_labels <- simulation$sample_phylogeny$package_clone_phylogeny$clone_phylogeny_labels
        clone_phylogeny_ID <- simulation$sample_phylogeny$package_clone_phylogeny$clone_phylogeny_ID
        clone_phylogeny_origin <- simulation$sample_phylogeny$package_clone_phylogeny$clone_phylogeny_origin
        clone_phylogeny_genotype <- simulation$sample_phylogeny$package_clone_phylogeny$clone_phylogeny_genotype
        clone_hclust_nodes <- simulation$sample_phylogeny$package_clone_phylogeny$clone_hclust_nodes
        clone_hclust_merge <- simulation$sample_phylogeny$package_clone_phylogeny$clone_hclust_merge
        N_clones <- length(clone_phylogeny_labels)

        if (length(clone_phylogeny_origin) != (2 * N_clones - 1)) {
            print("LENGTH OF CLONAL RECORD IS NOT 2*N_CLONES-1: REVISIT LATER")
            next
        }

        #---------------------------Initialize the list of clonal labels
        vec_clonal_labels <- rep("", length = length(clone_phylogeny_genotype))
        vec_clonal_labels[(length(vec_clonal_labels) - N_clones + 1):length(vec_clonal_labels)] <- clone_phylogeny_labels
        #--------------------Build the genotype list for each clone node
        clone_phylogeny_all_genotypes <- vector("list", length = length(clone_phylogeny_genotype))
        #   Initialize genotype lists for clone leaves
        for (clone in length(clone_phylogeny_genotype):(length(clone_phylogeny_genotype) - N_clones + 1)) {
            all_genotypes <- clone_phylogeny_genotype[clone]
            while (all_genotypes[1] != 0) {
                ancestor_genotype <- evolution_origin[all_genotypes[1]]
                all_genotypes <- c(ancestor_genotype, all_genotypes)
            }
            clone_phylogeny_all_genotypes[[clone]] <- all_genotypes
        }
        #   Update genotype lists with information from clone hclust
        for (clone_hclust_mother_node in 1:nrow(clone_hclust_merge)) {
            #   Find mother clone node's index in phylogeny our style
            clone_phylogeny_mother_node <- which(is.element(clone_hclust_nodes, clone_hclust_mother_node))
            #   Find daughter clone nodes' indices in phylogeny our style
            vec_clone_hclust_daughter_nodes <- clone_hclust_merge[clone_hclust_mother_node, ]
            vec_clone_phylogeny_daughter_nodes <- which(is.element(clone_hclust_nodes, vec_clone_hclust_daughter_nodes))
            clone_phylogeny_daughter_node_1 <- vec_clone_phylogeny_daughter_nodes[[1]]
            clone_phylogeny_daughter_node_2 <- vec_clone_phylogeny_daughter_nodes[[2]]
            #   Update genotype lists for the 3 clone nodes
            mother_genotypes <- intersect(clone_phylogeny_all_genotypes[[clone_phylogeny_daughter_node_1]], clone_phylogeny_all_genotypes[[clone_phylogeny_daughter_node_2]])
            clone_phylogeny_all_genotypes[[clone_phylogeny_mother_node]] <- mother_genotypes
            clone_phylogeny_all_genotypes[[clone_phylogeny_daughter_node_1]] <- setdiff(clone_phylogeny_all_genotypes[[clone_phylogeny_daughter_node_1]], mother_genotypes)
            clone_phylogeny_all_genotypes[[clone_phylogeny_daughter_node_2]] <- setdiff(clone_phylogeny_all_genotypes[[clone_phylogeny_daughter_node_2]], mother_genotypes)
            #   Update clonal labels if necessary
            if (length(clone_phylogeny_all_genotypes[[clone_phylogeny_daughter_node_1]]) == 0) {
                vec_clonal_labels[clone_phylogeny_mother_node] <- vec_clonal_labels[clone_phylogeny_daughter_node_1]
                vec_clonal_labels[clone_phylogeny_daughter_node_1] <- ""
            }
            if (length(clone_phylogeny_all_genotypes[[clone_phylogeny_daughter_node_2]]) == 0) {
                vec_clonal_labels[clone_phylogeny_mother_node] <- vec_clonal_labels[clone_phylogeny_daughter_node_2]
                vec_clonal_labels[clone_phylogeny_daughter_node_2] <- ""
            }
        }
        #-------------------------Find clonal populations as time series
        table_clonal_populations <- matrix(0, nrow = length(clone_phylogeny_all_genotypes), ncol = length(vec_time_simulation))
        vec_total_populations <- rep(0, length = length(vec_time_simulation))
        for (col in 1:length(vec_time_simulation)) {
            time <- vec_time_simulation[col]
            loc <- which.min(abs(evolution_traj_time - time))
            vec_clonal_ID <- evolution_traj_clonal_ID[[loc]]
            vec_clonal_population <- evolution_traj_population[[loc]]
            total_clonal_population <- sum(vec_clonal_population)
            for (row in 1:length(clone_phylogeny_all_genotypes)) {
                vec_loc <- which(is.element(vec_clonal_ID, clone_phylogeny_all_genotypes[[row]]))
                if (length(vec_loc) == 0) {
                    next
                }
                table_clonal_populations[row, col] <- sum(vec_clonal_population[vec_loc])
            }
            vec_total_populations[col] <- total_clonal_population
        }
        #------------------------------------------Find clonal parentage
        vec_clonal_parentage <- clone_phylogeny_origin
        #--------------------------------Add a clone for other genotypes
        table_clonal_populations <- rbind(rep(0, length = length(vec_time_simulation)), table_clonal_populations)
        for (col in 1:length(vec_time_simulation)) {
            table_clonal_populations[1, col] <- vec_total_populations[col] - sum(table_clonal_populations[, col])
        }
        vec_clonal_parentage <- c(0, (vec_clonal_parentage + 1))
        vec_clonal_labels <- c("Others", vec_clonal_labels)
        #--------------------------------------Remove unnecessary clones
        #---------------------------------i.e. clones that are always 0%
        #   Find unnecessary clones
        vec_unnecessary_clones <- c()
        for (clone in 1:nrow(table_clonal_populations)) {
            if (all(table_clonal_populations[clone, ] == 0)) {
                vec_unnecessary_clones <- c(vec_unnecessary_clones, clone)
            }
        }
        #   Find remaining clones
        vec_remaining_clones <- setdiff(1:length(vec_clonal_parentage), vec_unnecessary_clones)
        #   Rewire clones to new mother clones
        for (clone_daughter in 1:length(vec_clonal_parentage)) {
            clone_mother <- vec_clonal_parentage[clone_daughter]
            while ((clone_mother != 0) & (is.element(clone_mother, vec_unnecessary_clones))) {
                clone_mother <- vec_clonal_parentage[clone_mother]
            }
            vec_clonal_parentage[clone_daughter] <- clone_mother
        }
        #   Remove unnecessary clones from existence
        table_clonal_populations <- table_clonal_populations[-vec_unnecessary_clones, ]
        vec_clonal_parentage <- vec_clonal_parentage[-vec_unnecessary_clones]
        vec_clonal_labels <- vec_clonal_labels[-vec_unnecessary_clones]
        #   Correct clone indices for remaining clones
        for (node in 1:length(vec_clonal_parentage)) {
            mother_old <- vec_clonal_parentage[node]
            if (mother_old == 0) {
                next
            }
            mother_new <- which(vec_remaining_clones == mother_old)
            vec_clonal_parentage[node] <- mother_new
        }
        #---------Scale the clonal populations to match total population
        max_total_population <- max(vec_total_populations)
        for (col in 1:length(vec_time_simulation)) {
            #   Total population size can only go up to 90% to have leeway with
            #   numerical errors
            table_clonal_populations[, col] <- 49 * table_clonal_populations[, col] / max_total_population
        }
        #-------Conform clonal populations to nested format of fish plot
        table_clonal_populations_tmp <- table_clonal_populations
        for (clone_daughter in 1:length(vec_clonal_parentage)) {
            clone_mother <- vec_clonal_parentage[clone_daughter]
            while (clone_mother > 0) {
                #   Add in 0.001 to offset potential numerical errors to make
                #   sure mother clone always has more than sum of daughter clones
                table_clonal_populations[clone_mother, ] <- table_clonal_populations[clone_mother, ] + table_clonal_populations_tmp[clone_daughter, ] + 0.001
                clone_mother <- vec_clonal_parentage[clone_mother]
            }
        }
        #--Final check to make sure total population doesn't exceed 100%
        for (row in 1:nrow(table_clonal_populations)) {
            vec_row <- table_clonal_populations[row, ]
            vec_fix <- which(vec_row > 100)
            vec_row[vec_fix] <- 100
            table_clonal_populations[row, ] <- vec_row
        }
        #--------------------------------------Plot the clonal evolution
        #   Create fish object
        # fish <- createFishObject(table_clonal_populations, vec_clonal_parentage, timepoints = vec_time)
        fish <- createFishObject(table_clonal_populations, vec_clonal_parentage, timepoints = vec_time, clone.labels = vec_clonal_labels)
        #   Set clonal colors
        qual_col_pals <- brewer.pal.info[brewer.pal.info$category == "qual", ]
        col_vector <- unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
        cols <- sample(col_vector, nrow(table_clonal_populations))
        fish <- setCol(fish, cols)


        #
        vlines_pos <- Table_sampling$Age_sample
        vlines_sample_id <- Table_sampling$Sample_ID
        vlines_tit <- paste(vlines_sample_id, " (T=", as.character(vlines_pos), ")", sep = "")
        vlines_pos <- c(0, vlines_pos)
        vlines_tit <- c("T=0", vlines_tit)

        #   Create fish plot for clonal evolution
        fish <- layoutClones(fish)
        filename <- paste(model, "_sim", i, "_clonal_fishplot", ".jpeg", sep = "")
        jpeg(file = filename, width = width, height = height)
        p <- fishPlot(fish, shape = "polygon", vlines = vlines_pos, vlab = vlines_tit, cex.vlab = 3, bg.type = "solid")

        drawLegend(fish, xpos = -20, ypos = 25, cex = 4)

        print(p)
        dev.off()
    }
}
