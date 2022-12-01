# =====================================PLOT CLONAL EVOLUTION AS FISH PLOT
#' @export
plot_clonal_fishplot <- function(model = "",
                                 n_simulations = 0,
                                 vec_time = NULL,
                                 unit_time = "year",
                                 width = 1000,
                                 height = 500,
                                 compute_parallel = TRUE,
                                 n_cores = NULL) {
    library(RColorBrewer)
    library(fishplot)
    library(wesanderson)
    if (compute_parallel == FALSE) {
        #-----------------------Plot clonal evolution in sequential mode
        for (iteration in 1:n_simulations) {
            plot_clonal_fishplot_one_simulation(
                model,
                iteration,
                vec_time,
                unit_time,
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
        width <<- width
        height <<- height
        vec_time <<- vec_time
        unit_time <<- unit_time
        plot_clonal_fishplot_one_simulation <<- plot_clonal_fishplot_one_simulation
        clusterExport(cl, varlist = c(
            "plot_clonal_fishplot_one_simulation",
            "model",
            "vec_time",
            "unit_time",
            "width",
            "height"
        ))
        clusterEvalQ(cl = cl, require(RColorBrewer))
        clusterEvalQ(cl = cl, require(fishplot))
        clusterEvalQ(cl = cl, require(wesanderson))
        #   Plot in parallel
        pblapply(cl = cl, X = 1:n_simulations, FUN = function(iteration) {
            plot_clonal_fishplot_one_simulation(
                model,
                iteration,
                vec_time,
                unit_time,
                width,
                height
            )
        })
        #   Stop parallel cluster
        stopCluster(cl)
    }
}

plot_clonal_fishplot_one_simulation <- function(model,
                                                iteration,
                                                vec_time,
                                                unit_time,
                                                width,
                                                height) {
    #----------------------------------------------Input simulation file
    filename <- paste(model, "_simulation_", iteration, ".rda", sep = "")
    load(filename)
    #---------------------------------Transform time points if necessary
    if (is.null(vec_time)) {
        evolution_traj_time <- simulation$clonal_evolution$evolution_traj_time
        T_0 <- evolution_traj_time[1]
        T_end <- evolution_traj_time[length(evolution_traj_time)]
        vec_time_simulation <- seq(T_0, T_end, by = ((T_end - T_0) / 1000))
        if (unit_time == "day") {
            vec_time <- vec_time_simulation
        } else if (unit_time == "week") {
            vec_time <- vec_time_simulation / 7
        } else if (unit_time == "month") {
            vec_time <- vec_time_simulation / 30
        } else if (unit_time == "year") {
            vec_time <- vec_time_simulation / 365
        }
    } else {
        if (unit_time == "day") {
            vec_time_simulation <- vec_time
        } else if (unit_time == "week") {
            vec_time_simulation <- 7 * vec_time
        } else if (unit_time == "month") {
            vec_time_simulation <- 30 * vec_time
        } else if (unit_time == "year") {
            vec_time_simulation <- 365 * vec_time
        }
    }
    #-------------------------------------------Input the sampling table
    Table_sampling <- simulation$sample$Table_sampling
    #-----------------------------------------Input the clonal evolution
    evolution_origin <- simulation$clonal_evolution$evolution_origin
    evolution_traj_time <- simulation$clonal_evolution$evolution_traj_time
    evolution_traj_clonal_ID <- simulation$clonal_evolution$evolution_traj_clonal_ID
    evolution_traj_population <- simulation$clonal_evolution$evolution_traj_population
    #-----------------------------------------Input the clonal phylogeny
    clone_phylogeny_labels <- simulation$sample_phylogeny$package_clone_phylogeny$clone_phylogeny_labels
    clone_phylogeny_origin <- simulation$sample_phylogeny$package_clone_phylogeny$clone_phylogeny_origin
    clone_phylogeny_genotype <- simulation$sample_phylogeny$package_clone_phylogeny$clone_phylogeny_genotype
    clone_hclust_nodes <- simulation$sample_phylogeny$package_clone_phylogeny$clone_hclust_nodes
    clone_hclust_merge <- simulation$sample_phylogeny$package_clone_phylogeny$clone_hclust_merge
    N_clones <- length(clone_phylogeny_labels)
    if (N_clones <= 1) {
        print("CURRENTLY CANNOT PLOT CLONAL EVOLUTION FOR ONLY ONE CLONE")
        return()
    }
    #-------------------------------Initialize the list of clonal labels
    vec_clonal_labels <- rep("", length = length(clone_phylogeny_genotype))
    vec_clonal_labels[(length(vec_clonal_labels) - N_clones + 1):length(vec_clonal_labels)] <- clone_phylogeny_labels
    #------------------------Build the genotype list for each clone node
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
    #----------------------------------------------Find clonal parentage
    vec_clonal_parentage <- clone_phylogeny_origin
    #----------------------------------Add a clone for MRCA if necessary
    if (length(clone_phylogeny_all_genotypes[[1]]) > 2) {
        tmp <- clone_phylogeny_all_genotypes[[1]]
        clone_phylogeny_all_genotypes[[1]] <- tmp[3:length(tmp)]
        clone_phylogeny_all_genotypes <- c(list(tmp[1:2]), clone_phylogeny_all_genotypes)
        vec_clonal_parentage <- c(0, (vec_clonal_parentage + 1))
        vec_clonal_labels <- c("MRCA", vec_clonal_labels)
    }
    #-----------------------------Find clonal populations as time series
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
    #------------------------------------Add a clone for other genotypes
    table_clonal_populations <- rbind(rep(0, length = length(vec_time_simulation)), table_clonal_populations)
    for (col in 1:length(vec_time_simulation)) {
        table_clonal_populations[1, col] <- vec_total_populations[col] - sum(table_clonal_populations[, col])
    }
    vec_clonal_parentage <- c(0, (vec_clonal_parentage + 1))
    vec_clonal_labels <- c("Others", vec_clonal_labels)
    #------------------------------------------Remove unnecessary clones
    #-------------------------------------i.e. clones that are always 0%
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
    if (length(vec_unnecessary_clones) > 0) {
        table_clonal_populations <- table_clonal_populations[-vec_unnecessary_clones, ]
        vec_clonal_parentage <- vec_clonal_parentage[-vec_unnecessary_clones]
        vec_clonal_labels <- vec_clonal_labels[-vec_unnecessary_clones]
    }
    #   Correct clone indices for remaining clones
    for (node in 1:length(vec_clonal_parentage)) {
        mother_old <- vec_clonal_parentage[node]
        if (mother_old == 0) {
            next
        }
        mother_new <- which(vec_remaining_clones == mother_old)
        vec_clonal_parentage[node] <- mother_new
    }
    #-------------Scale the clonal populations to match total population
    max_total_population <- max(vec_total_populations)
    for (col in 1:length(vec_time_simulation)) {
        #   Total population size can only go up to 90% to have leeway with
        #   numerical errors
        table_clonal_populations[, col] <- 49 * table_clonal_populations[, col] / max_total_population
    }
    #-----------Conform clonal populations to nested format of fish plot
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
    #------Final check to make sure total population doesn't exceed 100%
    for (row in 1:nrow(table_clonal_populations)) {
        vec_row <- table_clonal_populations[row, ]
        vec_fix <- which(vec_row > 100)
        vec_row[vec_fix] <- 100
        table_clonal_populations[row, ] <- vec_row
    }
    #------------------------------------------Name the ancestral clones
    vec_ancestral_clones <- which(vec_clonal_labels == "")
    for (j in 1:length(vec_ancestral_clones)) {
        ancestral_clone <- vec_ancestral_clones[j]
        vec_children_ind <- which(vec_clonal_parentage == ancestral_clone)
        vec_children <- vec_clonal_labels[vec_children_ind]
        while ("" %in% vec_children) {
            vec_loc <- which(vec_children == "")
            vec_ind <- vec_children_ind[vec_loc]
            vec_children_ind <- vec_children_ind[-vec_loc]
            for (k in 1:length(vec_ind)) {
                vec_children_ind <- c(vec_children_ind, which(vec_clonal_parentage == vec_ind[k]))
            }
            vec_children <- vec_clonal_labels[vec_children_ind]
        }
        label <- paste("Ancestor(", vec_children[1], sep = "")
        if (length(vec_children) >= 2) {
            for (k in 2:length(vec_children)) {
                label <- paste(label, ", ", vec_children[k], sep = "")
            }
        }
        label <- paste(label, ")", sep = "")
        vec_clonal_labels[ancestral_clone] <- label
    }
    #-------------------------------------------------Reorder the clones
    #   Reorder labels alphabetically
    vec_clonal_labels_new <- c(clone_phylogeny_labels, sort(setdiff(vec_clonal_labels, clone_phylogeny_labels)))
    #   Reorder table of clonal percentages
    table_clonal_populations_new <- matrix(0, nrow = nrow(table_clonal_populations), ncol = ncol(table_clonal_populations))
    for (j in 1:length(vec_clonal_labels_new)) {
        label <- vec_clonal_labels_new[j]
        loc <- which(vec_clonal_labels == label)
        table_clonal_populations_new[j, ] <- table_clonal_populations[loc, ]
    }
    #   Reorder vector of clonal parentage
    vec_clonal_parentage_new <- rep(0, length = length(vec_clonal_parentage))
    for (j in 1:length(vec_clonal_labels_new)) {
        label <- vec_clonal_labels_new[j]
        loc <- which(vec_clonal_labels == label)
        parent_old <- vec_clonal_parentage[loc]
        if (parent_old == 0) {
            parent_new <- 0
        } else {
            parent_new <- which(vec_clonal_labels_new == vec_clonal_labels[parent_old])
        }
        vec_clonal_parentage_new[j] <- parent_new
    }
    #   Reorder everything
    vec_clonal_labels <- vec_clonal_labels_new
    table_clonal_populations <- table_clonal_populations_new
    vec_clonal_parentage <- vec_clonal_parentage_new
    #------------------------------------------Plot the clonal evolution
    filename <- paste(model, "_sim", iteration, "_clonal_fishplot", ".jpeg", sep = "")
    jpeg(file = filename, width = width, height = height)
    #   Find clonal order of birth times
    vec_clonal_birth_order <- rep(0, length(vec_clonal_labels))
    for (clone in 1:nrow(table_clonal_populations)) {
        vec_clonal_birth_order[clone] <- which(table_clonal_populations[clone, ] > 0)[1]
    }
    vec_clonal_order <- sort(vec_clonal_birth_order, index.return = TRUE)$ix
    #   Define clonal colors by order of birth times
    vec_color_order <- wes_palette("Zissou1", length(vec_clonal_order), type = "continuous")
    vec_cols <- rep(0, length(vec_color_order))
    for (i in 1:length(vec_clonal_order)) {
        vec_cols[vec_clonal_order[i]] <- vec_color_order[i]
    }

    vec_cols[which(vec_clonal_labels == "Others")] <- "gray"
    #---Create fish object
    fish <- createFishObject(table_clonal_populations, vec_clonal_parentage, timepoints = vec_time, col = vec_cols)
    #---Set vertical time lines
    vlines_pos <- Table_sampling$Age_sample
    vlines_tit <- paste(Table_sampling$Sample_ID, " (T=", as.character(vlines_pos), ")", sep = "")
    vlines_pos <- c(0, vlines_pos)
    vlines_tit <- c("T=0", vlines_tit)
    #---Create fish plot for clonal evolution
    fish <- layoutClones(fish)
    p <- fishPlot(fish,
        shape = "polygon",
        pad.left = 0,
        border = 0.1,
        col.border = "black",
        vlines = vlines_pos,
        vlab = vlines_tit,
        col.vline = "black",
        cex.vlab = 3,
        bg.type = "solid",
        bg.col = "white"
    )
    #---Draw legend
    clone_phylogeny_labels <- clone_phylogeny_labels[vec_clonal_order]

    ind_legend_leaf <- which(vec_clonal_labels %in% clone_phylogeny_labels)
    ind_legend_others <- setdiff(1:length(vec_clonal_labels), ind_legend_leaf)

    vec_legend_leaf_unsorted <- vec_clonal_labels[ind_legend_leaf]
    vec_legend_leaf <- clone_phylogeny_labels[which(clone_phylogeny_labels %in% vec_legend_leaf_unsorted)]

    vec_legend_others <- vec_clonal_labels[ind_legend_others]
    #   Draw legend for leaf clones
    n_clones_per_row <- 20
    x_start <- 0
    x_space <- ((vec_time[length(vec_time)] - x_start) / n_clones_per_row)
    x_pad <- 0.25 * x_space
    y_start <- 20
    y_space <- 5
    for (ind in 1:length(vec_legend_leaf)) {
        lab <- vec_legend_leaf[ind]
        col <- vec_cols[which(vec_clonal_labels == lab)]
        column <- ind %% n_clones_per_row
        if (column == 0) {
            column <- n_clones_per_row
            row <- ind %/% n_clones_per_row
        } else {
            row <- ind %/% n_clones_per_row + 1
        }
        x <- x_start + (column - 1) * x_space
        y <- y_start - (row - 1) * y_space
        p <- p + points(x = x, y = y, pch = 21, col = "black", bg = col, cex = 5)
        p <- p + text(x = x + x_pad, y = y, labels = lab, cex = 2, adj = c(0, NA))
    }
    #   Draw legend for other clones
    n_row_max <- 3 + y_start / y_space - row

    n_clones_per_row <- 5

    x_start <- 0
    x_space <- ((vec_time[length(vec_time)] - x_start) / n_clones_per_row)
    y_start <- y - y_space
    y_space <- 5
    for (ind in 1:length(vec_legend_others)) {
        lab <- vec_legend_others[ind]
        col <- vec_cols[which(vec_clonal_labels == lab)]
        row <- ind %% n_row_max
        if (row == 0) {
            row <- n_row_max
            column <- ind %/% n_row_max
        } else {
            column <- ind %/% n_row_max + 1
        }
        x <- x_start + (column - 1) * x_space
        y <- y_start - (row - 1) * y_space
        p <- p + points(x = x, y = y, pch = 21, col = "black", bg = col, cex = 5)
        p <- p + text(x = x + x_pad, y = y, labels = lab, cex = 2, adj = c(0, NA))
    }
    #---Print complete plot
    tmp <- print(p)
    dev.off()
}
