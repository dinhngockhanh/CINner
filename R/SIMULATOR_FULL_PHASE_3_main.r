# =============================PHASE 3: COPY-NUMBER PROFILES OF A SAMPLE
#' @export
SIMULATOR_FULL_PHASE_3_main <- function(package_clonal_evolution, package_sample, report_progress) {
    #-----------------------------------------Input the clonal evolution
    T_final <- package_clonal_evolution$T_current
    genotype_list_ploidy_chrom <- package_clonal_evolution$genotype_list_ploidy_chrom
    genotype_list_ploidy_block <- package_clonal_evolution$genotype_list_ploidy_block
    genotype_list_ploidy_allele <- package_clonal_evolution$genotype_list_ploidy_allele
    evolution_traj_time <- package_clonal_evolution$evolution_traj_time
    evolution_traj_divisions <- package_clonal_evolution$evolution_traj_divisions
    evolution_traj_clonal_ID <- package_clonal_evolution$evolution_traj_clonal_ID
    evolution_traj_population <- package_clonal_evolution$evolution_traj_population
    evolution_origin <- package_clonal_evolution$evolution_origin
    #---------------------------------------------------Input the sample
    sample_cell_ID <- package_sample$sample_cell_ID
    sample_clone_ID <- package_sample$sample_clone_ID
    sample_clone_ID_letters <- package_sample$sample_clone_ID_letters
    table_clone_ID_vs_letters <- package_sample$table_clone_ID_vs_letters
    sample_time <- package_sample$sample_time

    Table_sampling <- package_sample$Table_sampling
    all_sample_ID <- package_sample$all_sample_ID

    N_sample <- length(sample_cell_ID)

    if (standard_time_unit == "day") {
        Table_sampling$T_sample_phylo <- Table_sampling$T_sample_real
    } else if (standard_time_unit == "week") {
        T_final <- T_final / 7
        evolution_traj_time <- evolution_traj_time / 7
        Table_sampling$T_sample_phylo <- Table_sampling$T_sample_real / 7
    } else if (standard_time_unit == "month") {
        T_final <- T_final / 30
        evolution_traj_time <- evolution_traj_time / 30
        Table_sampling$T_sample_phylo <- Table_sampling$T_sample_real / 30
    } else if (standard_time_unit == "year") {
        T_final <- T_final / 365
        evolution_traj_time <- evolution_traj_time / 365
        Table_sampling$T_sample_phylo <- Table_sampling$T_sample_real / 365
    }
    #-------------------------------Initialize phylogeny in hclust style
    #   Initialize information to build phylogeny in hclust style
    hclust_row <- 0
    hclust_nodes <- rep(0, 1, 2 * N_sample - 1)
    hclust_nodes[N_sample:(2 * N_sample - 1)] <- (-1:-N_sample)
    hclust_labels <- sample_cell_ID
    #   Initialize extra information - genotypes of internal nodes
    hclust_internal_genotypes <- rep(-1, 1, N_sample - 1)
    #   Initialize extra information - table of CN events within nodes
    hclust_CN_events <- data.frame(matrix(ncol = 4, nrow = 0))
    colnames(hclust_CN_events) <- c("Mother_hclust", "Daughter_hclust", "Mother_genotype", "Daughter_genotype")
    #   Initialize actual phylogeny in hclust style
    hclust_merge <- matrix(0, nrow = N_sample - 1, ncol = 2)
    hclust_height <- rep(0, 1, N_sample - 1)
    #----------------------------------Initialize phylogeny in our style
    phylogeny_origin <- rep(0, length = 2 * N_sample - 1)
    phylogeny_elapsed_gens <- rep(0, length = 2 * N_sample - 1)
    phylogeny_elapsed_genotypes <- vector("list", length = 2 * N_sample - 1)
    phylogeny_genotype <- rep(0, length = 2 * N_sample - 1)
    phylogeny_birthtime <- rep(0, length = 2 * N_sample - 1)
    phylogeny_deathtime <- rep(0, length = 2 * N_sample - 1)
    #   Initialize data for leaves of sample phylogeny
    phylogeny_elapsed_gens[N_sample:(2 * N_sample - 1)] <- 1
    for (node in N_sample:(2 * N_sample - 1)) {
        phylogeny_elapsed_genotypes[[node]] <- sample_clone_ID[node - N_sample + 1]
        sample_ID <- all_sample_ID[node - N_sample + 1]
        phylogeny_deathtime[node] <- Table_sampling$T_sample_phylo[which(Table_sampling$Sample_ID == sample_ID)]
    }
    phylogeny_genotype[N_sample:(2 * N_sample - 1)] <- sample_clone_ID
    #-------------------------------Build the sample cell phylogeny tree
    node_mother_next <- N_sample - 1
    if (report_progress == TRUE) {
        pb <- txtProgressBar(
            min = 1, max = length(evolution_traj_divisions),
            style = 3, width = 50, char = "="
        )
    }
    for (i in seq(length(evolution_traj_divisions), 1, -1)) {
        if (report_progress == TRUE) {
            setTxtProgressBar(pb, length(evolution_traj_divisions) - i + 1)
        }
        #   Get time point
        time <- evolution_traj_time[i]
        #   Get list of eligible nodes
        current_node_list <- which((phylogeny_deathtime > time) & (phylogeny_birthtime <= time))
        current_node_genotype <- rep(0, length = length(current_node_list))
        for (j in 1:length(current_node_list)) {
            node <- current_node_list[j]
            current_node_genotype[j] <- phylogeny_elapsed_genotypes[[node]][1]
        }
        #   Get current clonal populations in total population
        eligible_clonal_ID <- evolution_traj_clonal_ID[[i + 1]]
        eligible_clonal_total_population <- evolution_traj_population[[i + 1]]
        #   Get current clonal populations in sample
        eligible_clonal_sample_population <- rep(0, length(eligible_clonal_ID))
        for (clone in 1:length(eligible_clonal_ID)) {
            clone_ID <- eligible_clonal_ID[clone]
            eligible_clonal_sample_population[clone] <- length(which(current_node_genotype == clone_ID))
        }
        #   Translate next clonal populations in total population as
        #   thresholds that clonal populations in sample cannot exceed
        if (i == 1) {
            limit_clonal_total_population <- rep(Inf, length(eligible_clonal_ID))
        } else {
            limit_clonal_total_population <- rep(0, length(eligible_clonal_ID))
            eligible_clonal_ID_tmp <- evolution_traj_clonal_ID[[i]]
            eligible_clonal_total_population_tmp <- evolution_traj_population[[i]]
            for (clone in 1:length(eligible_clonal_ID)) {
                clone_ID <- eligible_clonal_ID[clone]
                loc_tmp <- which(eligible_clonal_ID_tmp == clone_ID)
                if (length(loc_tmp) != 0) {
                    limit_clonal_total_population[clone] <- eligible_clonal_total_population_tmp[loc_tmp]
                }
            }
        }
        # =======Sanity tests
        if (sum(eligible_clonal_sample_population) != length(current_node_genotype)) {
            stop("Clonal populations in sample do not add up", call. = FALSE)
        }
        if (any(eligible_clonal_sample_population > eligible_clonal_total_population)) {
            stop("Clonal populations in sample are larger than in total cell population", call. = FALSE)
        }
        # ==Get list of divisions occurring in total population
        #   Column 1:       number of divisions
        #   Column 2:       genotype mother
        #   Column 3:       genotype daughter 1
        #   Column 4:       genotype daughter 2
        mat_division_total_population <- evolution_traj_divisions[[i]]
        if (length(mat_division_total_population) == 0) {
            next
        }
        # =======One huge loop to make sure clonal populations in sample are correct
        logic_correct <- 0
        while (logic_correct == 0) {
            #---Simulate identities of all divisions occurring in sample
            #   Row:            corresponding to mat_division_total_population
            #   Column 1:       node indices undergoing division as daughter 1
            #   Column 2:       node indices undergoing division as daughter 2
            #   Column 3:       division indices for corresponding nodes on column 1
            #   Column 4:       division indices for corresponding nodes on column 2
            #---Translation for R: linearized into vector:
            #   Entries (1:4)   corresponds to row 1
            #   Entries (5:8)   corresponds to row 2
            #   ...
            mat_division_sample <- vector("list", length = (4 * nrow(mat_division_total_population)))
            for (clone in 1:length(eligible_clonal_ID)) {
                #   For every clone found in the total population...
                #   Find its clonal ID
                clonal_ID <- eligible_clonal_ID[clone]
                #   Find its population in total population
                clonal_total_population <- eligible_clonal_total_population[clone]
                #   Find its population in sample's eligible nodes
                clonal_sample_population <- eligible_clonal_sample_population[clone]
                if (clonal_sample_population <= 0) {
                    next
                }
                #   Find all division roles that this clone plays in total population
                #   Row 1:      division index (= row in mat_division_total_population)
                #   Row 2:      daughter position (= 1 or 2)
                #   Row 3:      Cell count for this division/position in total population
                #   Row 4:      Node count for this division/position in sample ------> to be done in next sections
                mat_division_sample_clone <- c()
                for (daughter in 1:2) {
                    vec_division_genotype_daughter <- mat_division_total_population[, daughter + 2]
                    vec_division <- which(vec_division_genotype_daughter == clonal_ID)
                    mat_division_sample_clone_new <- rbind(
                        matrix(vec_division, nrow = 1),
                        matrix(rep(daughter, length(vec_division)), nrow = 1),
                        matrix(mat_division_total_population[vec_division, 1], nrow = 1)
                    )
                    mat_division_sample_clone <- cbind(mat_division_sample_clone, mat_division_sample_clone_new)
                }
                if (length(mat_division_sample_clone) == 0) {
                    next
                }
                #   Find total node count of this clone to undergo divisions
                #   of each type, i.e. row 4 in mat_division_sample_clone
                division_index_all <- mat_division_sample_clone[1, ]
                count_nodes_each_max <- mat_division_total_population[division_index_all, 1]
                freq <- sum(mat_division_sample_clone[3, ]) / clonal_total_population
                ########################################################
                ########################################################
                ########################################################
                #   Find total count of nodes to undergo divisions of all types
                while (1) {
                    count_nodes_all <- rbinom(n = 1, size = sum(mat_division_sample_clone[3, ]), clonal_sample_population / clonal_total_population)
                    if (count_nodes_all <= clonal_sample_population) break
                }
                #   Divide total count of nodes among different division types
                count_nodes_each <- rep(-1, ncol(mat_division_sample_clone))
                while (min(count_nodes_each) < 0) {
                    count_nodes_each <- mat_division_sample_clone[3, ] - rmultinom(n = 1, size = (sum(mat_division_sample_clone[3, ]) - count_nodes_all), mat_division_sample_clone[3, ] / sum(mat_division_sample_clone[3, ]))
                }
                count_nodes_each <- matrix(count_nodes_each, nrow = 1)
                mat_division_sample_clone <- rbind(mat_division_sample_clone, count_nodes_each)
                ########################################################
                ########################################################
                ########################################################
                # #   Find total count of nodes to undergo divisions of all types
                # count_nodes_all <- rbinom(n = 1, size = clonal_sample_population, prob = freq)
                # #   Divide total count of nodes among different division types
                # count_nodes_each <- rmultinom(n = 1, size = count_nodes_all, mat_division_sample_clone[3, ] / sum(mat_division_sample_clone[3, ]))
                # count_nodes_each <- matrix(count_nodes_each, nrow = 1)
                # mat_division_sample_clone <- rbind(mat_division_sample_clone, count_nodes_each)
                ########################################################
                ########################################################
                ########################################################
                #   Check that node count in each position doesn't exceed limit in total population
                if (any(count_nodes_each > count_nodes_each_max)) {
                    logic_correct <- -1
                    break
                }
                #   Jump to next clone if there is no division to perform
                if (max(mat_division_sample_clone[4, ]) == 0) {
                    next
                }
                #   Simulate which nodes undergo each division type
                #   i.e. columns 1 & 2 in mat_division_sample
                eligible_nodes <- current_node_list[which(current_node_genotype == clonal_ID)]
                if (length(eligible_nodes) == 1) {
                    node_indices_all <- eligible_nodes
                } else {
                    node_indices_all <- sample(x = eligible_nodes, size = count_nodes_all, replace = FALSE)
                }
                for (division_type in 1:ncol(mat_division_sample_clone)) {
                    row <- mat_division_sample_clone[1, division_type]
                    col <- mat_division_sample_clone[2, division_type]
                    count <- mat_division_sample_clone[4, division_type]
                    if (count > 0) {
                        ind <- (row - 1) * 4 + col
                        mat_division_sample[[ind]] <- node_indices_all[1:count]
                        node_indices_all <- node_indices_all[-(1:count)]
                    }
                }
                #   Simulate the division indices for each division type
                #   i.e. columns 3 & 4 in mat_division_sample
                for (division_type in 1:ncol(mat_division_sample_clone)) {
                    row <- mat_division_sample_clone[1, division_type]
                    col <- mat_division_sample_clone[2, division_type]
                    count <- mat_division_sample_clone[4, division_type]
                    if (count > 0) {
                        col <- col + 2
                        ind <- (row - 1) * 4 + col
                        count_divisions_total <- mat_division_total_population[row, 1]
                        mat_division_sample[[ind]] <- sample(x = count_divisions_total, size = count, replace = FALSE)
                    }
                }
            }
            #   Redo the whole process if some node count in some position exceeded limit in total population
            if (logic_correct == -1) {
                logic_correct <- 0
                next
            }
            #---Update phylogeny tree according to the division identities
            #   Save the current phylogeny in case new changes are wrong

            node_mother_next_tmp <- node_mother_next

            hclust_row_tmp <- hclust_row
            hclust_nodes_tmp <- hclust_nodes
            hclust_merge_tmp <- hclust_merge
            hclust_height_tmp <- hclust_height
            hclust_internal_genotypes_tmp <- hclust_internal_genotypes
            hclust_CN_events_tmp <- hclust_CN_events

            phylogeny_origin_tmp <- phylogeny_origin
            phylogeny_elapsed_gens_tmp <- phylogeny_elapsed_gens
            phylogeny_elapsed_genotypes_tmp <- phylogeny_elapsed_genotypes
            phylogeny_genotype_tmp <- phylogeny_genotype
            phylogeny_birthtime_tmp <- phylogeny_birthtime
            phylogeny_deathtime_tmp <- phylogeny_deathtime

            current_node_list_tmp <- current_node_list
            current_node_genotype_tmp <- current_node_genotype
            #   Update phylogeny according to every division
            for (division_type in 1:nrow(mat_division_total_population)) {
                genotype_mother <- mat_division_total_population[division_type, 2]
                #   Get list of nodes in positions of daughter 1 and daughter 2
                ind_1 <- (division_type - 1) * 4 + 1
                vec_nodes_daughter_1 <- mat_division_sample[[ind_1]]
                ind_2 <- (division_type - 1) * 4 + 2
                vec_nodes_daughter_2 <- mat_division_sample[[ind_2]]
                if ((length(vec_nodes_daughter_1) == 0) && (length(vec_nodes_daughter_2) == 0)) {
                    next
                }
                #   Get list of division indices of daughter 1 and daughter 2
                ind_1 <- (division_type - 1) * 4 + 3
                vec_div_indices_1 <- mat_division_sample[[ind_1]]
                ind_2 <- (division_type - 1) * 4 + 4
                vec_div_indices_2 <- mat_division_sample[[ind_2]]
                vec_div_indices_all <- unique(c(vec_div_indices_1, vec_div_indices_2))
                #   Perform each division
                for (division in 1:length(vec_div_indices_all)) {
                    div_index <- vec_div_indices_all[division]
                    loc_1 <- which(vec_div_indices_1 == div_index)
                    loc_2 <- which(vec_div_indices_2 == div_index)
                    if ((length(loc_1) != 0) && (length(loc_2) != 0)) {
                        #   Nodes 1 and 2 are mergning...
                        node_1 <- vec_nodes_daughter_1[loc_1]
                        node_2 <- vec_nodes_daughter_2[loc_2]

                        node_mother <- node_mother_next
                        node_mother_next <- node_mother_next - 1
                        #   Update phylogeny in hclust style
                        hclust_row <- hclust_row + 1
                        hclust_nodes[node_mother] <- hclust_row
                        hclust_merge[hclust_row, ] <- c(hclust_nodes[node_1], hclust_nodes[node_2])
                        hclust_height[hclust_row] <- T_final - time
                        #   Update internal nodes' genotypes in hclust style
                        hclust_internal_genotypes[hclust_row] <- genotype_mother
                        #   Update table of CN events for internal nodes (if necessary)
                        genotype_daughter_1 <- phylogeny_genotype[node_1]
                        genotype_daughter_2 <- phylogeny_genotype[node_2]
                        if (genotype_daughter_1 != genotype_mother) {
                            hclust_CN_events[nrow(hclust_CN_events) + 1, ] <- c(hclust_row, hclust_nodes[node_1], genotype_mother, genotype_daughter_1)
                        }
                        if (genotype_daughter_2 != genotype_mother) {
                            hclust_CN_events[nrow(hclust_CN_events) + 1, ] <- c(hclust_row, hclust_nodes[node_2], genotype_mother, genotype_daughter_2)
                        }
                        #   Update phylogeny in our style
                        phylogeny_origin[node_1] <- node_mother
                        phylogeny_origin[node_2] <- node_mother
                        phylogeny_elapsed_gens[node_mother] <- 1
                        phylogeny_elapsed_genotypes[[node_mother]] <- c(genotype_mother)
                        phylogeny_genotype[node_mother] <- genotype_mother
                        phylogeny_birthtime[node_1] <- time
                        phylogeny_birthtime[node_2] <- time
                        phylogeny_deathtime[node_mother] <- time
                        #   Update phylogeny records in our style
                        pos_delete <- c(which(current_node_list == node_1), which(current_node_list == node_2))
                        current_node_list <- current_node_list[-pos_delete]
                        current_node_list <- c(node_mother, current_node_list)
                        current_node_genotype <- current_node_genotype[-pos_delete]
                        current_node_genotype <- c(genotype_mother, current_node_genotype)
                    } else {
                        #   Either node 1 or node 2 has one more division...
                        if (length(loc_1) != 0) {
                            node_daughter <- vec_nodes_daughter_1[loc_1]
                        } else {
                            node_daughter <- vec_nodes_daughter_2[loc_2]
                        }
                        #   Update phylogeny in our style
                        phylogeny_elapsed_gens[node_daughter] <- phylogeny_elapsed_gens[node_daughter] + 1
                        phylogeny_elapsed_genotypes[[node_daughter]] <- c(genotype_mother, phylogeny_elapsed_genotypes[[node_daughter]])
                        #   Update phylogeny records in our style
                        loc_daughter <- which(current_node_list == node_daughter)
                        current_node_genotype[loc_daughter] <- genotype_mother
                    }
                }
            }
            #---Check if the clonal populations in sample satisfy conditions
            #   Find clonal populations in sample after new divisions
            tmp_clonal_sample_population <- rep(0, length(eligible_clonal_ID))
            for (clone in 1:length(eligible_clonal_ID)) {
                clone_ID <- eligible_clonal_ID[clone]
                tmp_clonal_sample_population[clone] <- length(which(current_node_genotype == clone_ID))
            }
            #   Redo this whole step if clonal populations violate thresholds
            if (any(tmp_clonal_sample_population > limit_clonal_total_population)) {
                node_mother_next <- node_mother_next_tmp

                hclust_row <- hclust_row_tmp
                hclust_nodes <- hclust_nodes_tmp
                hclust_merge <- hclust_merge_tmp
                hclust_height <- hclust_height_tmp
                hclust_internal_genotypes <- hclust_internal_genotypes_tmp
                hclust_CN_events <- hclust_CN_events_tmp

                phylogeny_origin <- phylogeny_origin_tmp
                phylogeny_elapsed_gens <- phylogeny_elapsed_gens_tmp
                phylogeny_elapsed_genotypes <- phylogeny_elapsed_genotypes_tmp
                phylogeny_genotype <- phylogeny_genotype_tmp
                phylogeny_birthtime <- phylogeny_birthtime_tmp
                phylogeny_deathtime <- phylogeny_deathtime_tmp

                current_node_list <- current_node_list_tmp
                current_node_genotype <- current_node_genotype_tmp
            } else {
                logic_correct <- 1
            }
        }
    }
    if (report_progress == TRUE) {
        cat("\n")
    }
    #----------------------------------------Complete the unmerged nodes
    #-------------------------i.e. there is more than one ancestral cell
    #   Find all unmerged nodes
    list_unmerged_nodes <- which(phylogeny_origin == 0 & hclust_nodes != 0)
    list_unnecessary_nodes <- which(phylogeny_origin == 0 & hclust_nodes == 0)
    N_unnecessary_nodes <- length(list_unnecessary_nodes)
    #---Complete the phylogeny in hclust style
    node_anchor <- list_unmerged_nodes[1]
    hclust_node_anchor <- hclust_nodes[node_anchor]
    #   Merge all unmerged nodes together at first time point
    if (length(list_unmerged_nodes) >= 2) {
        for (i in 2:length(list_unmerged_nodes)) {
            node <- list_unmerged_nodes[i]
            hclust_row <- hclust_row + 1
            hclust_merge[hclust_row, ] <- c(hclust_node_anchor, hclust_nodes[node])
            hclust_node_anchor <- hclust_node_anchor + 1
            hclust_height[hclust_row] <- T_final
            hclust_internal_genotypes[hclust_row] <- 0
        }
    }
    #---Complete the phylogeny in our style
    #   Merge all unmerged nodes together at first time point
    phylogeny_birthtime[list_unmerged_nodes] <- evolution_traj_time[1]
    #   Delete unnecessary nodes
    phylogeny_origin <- phylogeny_origin - N_unnecessary_nodes
    phylogeny_origin[list_unmerged_nodes] <- 0
    if (length(list_unnecessary_nodes) > 0) {
        hclust_nodes <- hclust_nodes[-list_unnecessary_nodes]

        phylogeny_origin <- phylogeny_origin[-list_unnecessary_nodes]
        phylogeny_elapsed_gens <- phylogeny_elapsed_gens[-list_unnecessary_nodes]
        phylogeny_elapsed_genotypes <- phylogeny_elapsed_genotypes[-list_unnecessary_nodes]
        phylogeny_genotype <- phylogeny_genotype[-list_unnecessary_nodes]
        phylogeny_birthtime <- phylogeny_birthtime[-list_unnecessary_nodes]
        phylogeny_deathtime <- phylogeny_deathtime[-list_unnecessary_nodes]
    }
    #----------------------Add CN events for root node(s) (if necessary)
    list_hclust_roots_index <- hclust_nodes[list_unmerged_nodes - N_unnecessary_nodes]
    list_hclust_roots_genotype <- phylogeny_genotype[list_unmerged_nodes - N_unnecessary_nodes]
    for (root in 1:length(list_hclust_roots_index)) {
        root_hclust <- list_hclust_roots_index[root]
        genotype <- list_hclust_roots_genotype[root]
        mother_genotype <- evolution_origin[genotype]
        if (mother_genotype != 0) {
            grandmother_genotype <- evolution_origin[mother_genotype]
            while (grandmother_genotype != 0) {
                mother_genotype <- grandmother_genotype
                grandmother_genotype <- evolution_origin[grandmother_genotype]
            }
            hclust_CN_events[nrow(hclust_CN_events) + 1, ] <- c(0, root_hclust, mother_genotype, genotype)
        }
    }
    #-------------------------------------Reorder the nodes for plotting
    list_roots <- list_unmerged_nodes - N_unnecessary_nodes
    #---Find an order on all nodes of the phylogeny in our style
    #   Find number of progeny of each node
    progeny_count <- rep(0, length(phylogeny_origin))
    end <- length(progeny_count)
    progeny_count[(end - N_sample + 1):end] <- 1
    for (node in length(progeny_count):1) {
        mother_node <- phylogeny_origin[node]
        if (mother_node > 0) {
            progeny_count[mother_node] <- progeny_count[mother_node] + progeny_count[node]
        }
    }
    #   Reorder the sample phylogeny tree based on progeny counts
    phylogeny_order <- rep(0, length(phylogeny_origin))
    phylogeny_order[list_roots] <- 1
    for (node in 0:length(progeny_count)) {
        vec_daughter_nodes <- which(phylogeny_origin == node)
        if (length(vec_daughter_nodes) == 0) {
            next
        }
        vec_progeny_counts <- progeny_count[vec_daughter_nodes]
        tmp <- sort(vec_progeny_counts, index.return = TRUE)
        vec_progeny_counts <- tmp$x
        vec_order <- tmp$ix
        vec_daughter_nodes <- vec_daughter_nodes[vec_order]
        for (i in 1:length(vec_daughter_nodes)) {
            daughter_node <- vec_daughter_nodes[i]
            if (i > 1) {
                progeny_count_extra <- sum(vec_progeny_counts[1:i - 1])
            } else {
                progeny_count_extra <- 0
            }
            if (node == 0) {
                phylogeny_order[daughter_node] <- phylogeny_order[daughter_node] + progeny_count_extra
            } else {
                phylogeny_order[daughter_node] <- phylogeny_order[node] + progeny_count_extra
            }
        }
    }
    #---Extract the order for phylogeny in hclust style
    hclust_order_inverse <- phylogeny_order[(length(phylogeny_order) - N_sample + 1):length(phylogeny_order)]
    hclust_order <- rep(0, N_sample)
    for (i_cell in 1:N_sample) {
        loc <- hclust_order_inverse[i_cell]
        hclust_order[loc] <- i_cell
    }
    #--------------------------------------------Create clustering table
    hclust_clustering <- data.frame(sample_cell_ID, sample_clone_ID_letters)
    names(hclust_clustering) <- c("cell_id", "clone_id")
    #----------------------------Create phylogeny object in hclust style
    hclust_height <- 2 * hclust_height
    #   Create phylogeny object in hclust style
    phylogeny_hclust <- list()
    phylogeny_hclust$merge <- hclust_merge
    phylogeny_hclust$height <- hclust_height
    phylogeny_hclust$order <- hclust_order
    phylogeny_hclust$labels <- sample_cell_ID
    class(phylogeny_hclust) <- "hclust"
    #-----------------------------Create phylogeny object in phylo style
    #   Create phylogeny object in phylo style
    phylogeny_phylo <- ape::as.phylo(phylogeny_hclust, use.labels = TRUE)
    hclust_height <- hclust_height / 2
    #---------Adjust the leaf lengths in phylo according to sample times
    #   Find the sample time for all leaves
    vec_sample_id <- sub("-.*", "", phylogeny_phylo$tip.label)
    vec_T_sample <- Table_sampling$Age_sample[match(vec_sample_id, Table_sampling$Sample_ID)]
    #   Find the location of leaf mergings in phylo
    vec_leaf_row <- match(1:length(phylogeny_phylo$tip.label), phylogeny_phylo$edge[, 1])
    vec_leaf_row_2 <- match(1:length(phylogeny_phylo$tip.label), phylogeny_phylo$edge[, 2])
    vec_leaf_row[which(is.na(vec_leaf_row))] <- vec_leaf_row_2[which(is.na(vec_leaf_row))]
    #   Substract the duration from length of leaf mergings in phylo
    phylogeny_phylo$edge.length[vec_leaf_row] <- phylogeny_phylo$edge.length[vec_leaf_row] - (T_final - vec_T_sample)
    # for (leaf in 1:length(phylogeny_phylo$tip.label)) {
    #     #   Find the sample time for this leaf
    #     sample_id <- strsplit(phylogeny_phylo$tip.label[leaf], "-")[[1]][1]
    #     T_sample <- Table_sampling$Age_sample[which(Table_sampling$Sample_ID == sample_id)]
    #     #   Find the location of leaf merging in phylo
    #     leaf_row <- 0
    #     for (row in 1:nrow(phylogeny_phylo$edge)) {
    #         if (phylogeny_phylo$edge[row, 1] == leaf | phylogeny_phylo$edge[row, 2] == leaf) {
    #             leaf_row <- row
    #         }
    #     }
    #     #   Substract the duration from length of leaf merging in phylo
    #     phylogeny_phylo$edge.length[leaf_row] <- phylogeny_phylo$edge.length[leaf_row] - (T_final - T_sample)
    # }
    #   Create object containing both phylo-style tree and clustering
    phylogeny_clustering_truth <- list()
    phylogeny_clustering_truth$tree <- phylogeny_phylo
    phylogeny_clustering_truth$clustering <- hclust_clustering
    #-------------------------------------Build the clone phylogeny tree
    clone_phylogeny_labels <- table_clone_ID_vs_letters$Clone_ID_letter
    clone_phylogeny_ID <- table_clone_ID_vs_letters$Clone_ID_number
    N_clones <- length(clone_phylogeny_labels)
    #---Initialize clone phylogeny in hclust style
    #   Initialize information to build clone phylogeny in hclust style
    clone_hclust_row <- 0
    clone_hclust_nodes <- rep(0, 1, 2 * N_clones - 1)
    clone_hclust_nodes[N_clones:(2 * N_clones - 1)] <- (-1:-N_clones)
    clone_hclust_labels <- clone_phylogeny_labels
    #   Initialize actual clone phylogeny in hclust style
    clone_hclust_merge <- matrix(0, nrow = N_clones - 1, ncol = 2)
    clone_hclust_height <- rep(0, 1, N_clones - 1)
    #---Initialize clone phylogeny in our style
    clone_phylogeny_origin <- rep(0, length = 2 * N_clones - 1)
    clone_phylogeny_genotype <- rep(0, length = 2 * N_clones - 1)
    clone_phylogeny_birthtime <- rep(0, length = 2 * N_clones - 1)
    clone_phylogeny_deathtime <- rep(0, length = 2 * N_clones - 1)
    #   Initialize the current list of node genotypes
    clone_current_node_genotype <- clone_phylogeny_ID
    #   Initialize the current list of nodes in the clone phylogeny
    clone_current_node_list <- N_clones:(2 * N_clones - 1)
    #   Initialize data for leaves of clone phylogeny
    for (node in N_clones:(2 * N_clones - 1)) {
        clone_phylogeny_genotype[node] <- clone_phylogeny_ID[node - N_clones + 1]
        clone_phylogeny_deathtime[node] <- T_final
    }
    #   Build the clone phylogeny tree
    for (hclust_mother_cell_node in 1:nrow(hclust_merge)) {
        #   Get daughter cells' indices in hclust style
        hclust_daughter_cell_nodes <- hclust_merge[hclust_mother_cell_node, ]
        #   Translate into daughter cells' indices in our style
        phylogeny_daughter_cell_nodes <- which(is.element(hclust_nodes, hclust_daughter_cell_nodes))

        if (length(phylogeny_daughter_cell_nodes) == 1) {
            next
        }

        #   Find daughter cells' genotypes
        genotype_daughter_cell_nodes <- phylogeny_genotype[phylogeny_daughter_cell_nodes]
        #   Find daughter cells' indices in clone hclust
        clone_phylogeny_daughter_nodes <- rep(0, length(genotype_daughter_cell_nodes))
        loc_1 <- which(clone_current_node_genotype == genotype_daughter_cell_nodes[1])
        clone_phylogeny_daughter_nodes[1] <- clone_current_node_list[loc_1]
        loc_2 <- which(clone_current_node_genotype == genotype_daughter_cell_nodes[2])
        clone_phylogeny_daughter_nodes[2] <- clone_current_node_list[loc_2]
        #   Update clone phylogeny...
        cell_node_1 <- phylogeny_daughter_cell_nodes[1]
        cell_node_2 <- phylogeny_daughter_cell_nodes[2]
        cell_node_mother <- phylogeny_origin[cell_node_1]
        clone_node_1 <- clone_phylogeny_daughter_nodes[1]
        clone_node_2 <- clone_phylogeny_daughter_nodes[2]
        if (clone_node_1 == clone_node_2) {
            #   If the cell merging happens within the same clone...
            if (cell_node_mother > 0) {
                genotype_mother <- phylogeny_genotype[cell_node_mother]
            } else {
                genotype_mother <- phylogeny_genotype[cell_node_1]
            }
            clone_node_mother <- clone_node_1
            #   Update collection of genotypes for this clone
            daughter_elapsed_genotypes <- unique(c(phylogeny_elapsed_genotypes[[cell_node_1]], phylogeny_elapsed_genotypes[[cell_node_2]]))

            clone_current_node_genotype[loc_1] <- genotype_mother
        } else {
            #   If the cell merging happens between different clones...
            if (cell_node_mother > 0) {
                genotype_mother <- phylogeny_genotype[cell_node_mother]
            } else {
                genotype_mother <- 0
            }
            # genotype_mother <- phylogeny_genotype[cell_node_mother]
            clone_node_mother <- min(clone_current_node_list) - 1
            #   Update clone phylogeny in hclust style
            clone_hclust_row <- clone_hclust_row + 1
            clone_hclust_nodes[clone_node_mother] <- clone_hclust_row
            clone_hclust_merge[clone_hclust_row, ] <- c(clone_hclust_nodes[clone_node_1], clone_hclust_nodes[clone_node_2])
            clone_hclust_height[clone_hclust_row] <- hclust_height[hclust_mother_cell_node]
            #   Update phylogeny in our style
            clone_phylogeny_origin[clone_node_1] <- clone_node_mother
            clone_phylogeny_origin[clone_node_2] <- clone_node_mother

            clone_phylogeny_genotype[clone_node_mother] <- genotype_mother

            clone_phylogeny_birthtime[clone_node_1] <- T_final - hclust_height[hclust_mother_cell_node]
            clone_phylogeny_birthtime[clone_node_2] <- T_final - hclust_height[hclust_mother_cell_node]
            clone_phylogeny_deathtime[clone_node_mother] <- T_final - hclust_height[hclust_mother_cell_node]
            #   Update clone phylogeny records in our style

            pos_delete <- c(which(clone_current_node_list == clone_node_1), which(clone_current_node_list == clone_node_2))

            clone_current_node_list <- clone_current_node_list[-pos_delete]
            clone_current_node_list <- c(clone_node_mother, clone_current_node_list)

            clone_current_node_genotype <- clone_current_node_genotype[-pos_delete]
            clone_current_node_genotype <- c(genotype_mother, clone_current_node_genotype)
        }
        #   Create another clone merging if there are repeated genotypes in the current records
        if (length(unique(clone_current_node_genotype)) < length(clone_current_node_genotype)) {
            #   Find the genotype that has to be resolved
            unique_genotypes_current <- unique(clone_current_node_genotype)
            for (i in 1:length(unique_genotypes_current)) {
                if (length(which(clone_current_node_genotype == unique_genotypes_current[i])) > 1) {
                    genotype_resolve <- unique_genotypes_current[i]
                    vec_clone_nodes <- which(clone_current_node_genotype == genotype_resolve)
                    clone_node_1 <- clone_current_node_list[vec_clone_nodes[1]]
                    clone_node_2 <- clone_current_node_list[vec_clone_nodes[2]]

                    clone_node_mother <- min(clone_current_node_list) - 1

                    clone_hclust_row <- clone_hclust_row + 1
                    clone_hclust_nodes[clone_node_mother] <- clone_hclust_row
                    clone_hclust_merge[clone_hclust_row, ] <- c(clone_hclust_nodes[clone_node_1], clone_hclust_nodes[clone_node_2])
                    clone_hclust_height[clone_hclust_row] <- hclust_height[hclust_mother_cell_node]

                    clone_phylogeny_origin[clone_node_1] <- clone_node_mother
                    clone_phylogeny_origin[clone_node_2] <- clone_node_mother

                    clone_phylogeny_genotype[clone_node_mother] <- genotype_resolve

                    clone_phylogeny_birthtime[clone_node_1] <- T_final - hclust_height[hclust_mother_cell_node]
                    clone_phylogeny_birthtime[clone_node_2] <- T_final - hclust_height[hclust_mother_cell_node]
                    clone_phylogeny_deathtime[clone_node_mother] <- T_final - hclust_height[hclust_mother_cell_node]

                    pos_delete <- c(which(clone_current_node_list == clone_node_1), which(clone_current_node_list == clone_node_2))

                    clone_current_node_list <- clone_current_node_list[-pos_delete]
                    clone_current_node_list <- c(clone_node_mother, clone_current_node_list)

                    clone_current_node_genotype <- clone_current_node_genotype[-pos_delete]
                    clone_current_node_genotype <- c(genotype_mother, clone_current_node_genotype)

                    break
                }
            }
        }
    }
    #----------------------------------------Complete the unmerged nodes
    #------------------------i.e. there is more than one ancestral clone
    #   Merge all unmerged nodes together at first time point
    list_unmerged_nodes <- which(clone_phylogeny_origin == 0 & clone_hclust_nodes != 0)
    clone_node_1 <- list_unmerged_nodes[1]
    # clone_hclust_node_anchor <- clone_hclust_nodes[clone_node_1]
    if (length(list_unmerged_nodes) >= 2) {
        for (i in 2:length(list_unmerged_nodes)) {
            clone_node_2 <- list_unmerged_nodes[i]

            clone_node_mother <- min(clone_current_node_list) - 1

            clone_hclust_row <- clone_hclust_row + 1
            clone_hclust_nodes[clone_node_mother] <- clone_hclust_row
            clone_hclust_merge[clone_hclust_row, ] <- c(clone_hclust_nodes[clone_node_1], clone_hclust_nodes[clone_node_2])
            clone_hclust_height[clone_hclust_row] <- T_final

            clone_phylogeny_origin[clone_node_1] <- clone_node_mother
            clone_phylogeny_origin[clone_node_2] <- clone_node_mother

            clone_phylogeny_genotype[clone_node_mother] <- min(clone_phylogeny_origin[clone_node_1], clone_phylogeny_genotype[clone_node_2])

            clone_phylogeny_birthtime[clone_node_1] <- evolution_traj_time[1]
            clone_phylogeny_birthtime[clone_node_2] <- evolution_traj_time[1]
            clone_phylogeny_deathtime[clone_node_mother] <- evolution_traj_time[1]

            pos_delete <- c(which(clone_current_node_list == clone_node_2))

            clone_current_node_list <- clone_current_node_list[-pos_delete]
            clone_current_node_list <- c(clone_node_mother, clone_current_node_list)

            clone_node_1 <- clone_node_mother
        }
        # for (i in 2:length(list_unmerged_nodes)) {
        #     clone_node <- list_unmerged_nodes[i]
        #     clone_hclust_row <- clone_hclust_row + 1
        #     clone_hclust_merge[clone_hclust_row, ] <- c(clone_hclust_node_anchor, clone_hclust_nodes[clone_node])
        #     clone_hclust_node_anchor <- clone_hclust_row
        #     clone_hclust_height[clone_hclust_row] <- T_final
        # }
    }
    # #   Delete unnecessary nodes
    # list_unnecessary_nodes <- which(clone_phylogeny_origin == 0 & clone_hclust_nodes == 0)
    # N_unnecessary_nodes <- length(list_unnecessary_nodes)
    #
    # clone_phylogeny_origin <- clone_phylogeny_origin - N_unnecessary_nodes
    # clone_phylogeny_origin[list_unmerged_nodes] <- 0
    # if (length(list_unnecessary_nodes) > 0) {
    #     clone_hclust_nodes <- clone_hclust_nodes[-list_unnecessary_nodes]
    #
    #     clone_phylogeny_origin <- clone_phylogeny_origin[-list_unnecessary_nodes]
    #     clone_phylogeny_genotype <- clone_phylogeny_genotype[-list_unnecessary_nodes]
    #     clone_phylogeny_birthtime <- clone_phylogeny_birthtime[-list_unnecessary_nodes]
    #     clone_phylogeny_deathtime <- clone_phylogeny_deathtime[-list_unnecessary_nodes]
    # }
    #-----------------------Create clone phylogeny object in phylo style
    if (N_clones > 1) {
        clone_hclust_height <- 2 * clone_hclust_height
        #       Create clone phylogeny object in hclust style
        clone_phylogeny_hclust <- list()
        clone_phylogeny_hclust$merge <- clone_hclust_merge
        clone_phylogeny_hclust$height <- clone_hclust_height
        clone_phylogeny_hclust$order <- 1:N_clones
        # clone_phylogeny_hclust$order                        <- clone_hclust_order
        clone_phylogeny_hclust$labels <- clone_hclust_labels
        class(clone_phylogeny_hclust) <- "hclust"
        #       Create clone phylogeny object in phylo style
        clone_phylogeny_phylo <- ape::as.phylo(clone_phylogeny_hclust, use.labels = TRUE)
        clone_hclust_height <- clone_hclust_height / 2
    } else {
        clone_phylogeny_phylo <- list()
        clone_phylogeny_hclust <- list()
    }
    #-----------------------------Output package of data from simulation
    package_cell_phylogeny_hclust_extra <- list()
    package_cell_phylogeny_hclust_extra$hclust_internal_genotypes <- hclust_internal_genotypes
    package_cell_phylogeny_hclust_extra$hclust_CN_events <- hclust_CN_events

    package_cell_phylogeny <- list()
    package_cell_phylogeny$phylogeny_origin <- phylogeny_origin
    package_cell_phylogeny$phylogeny_elapsed_gens <- phylogeny_elapsed_gens
    package_cell_phylogeny$phylogeny_elapsed_genotypes <- phylogeny_elapsed_genotypes
    package_cell_phylogeny$phylogeny_genotype <- phylogeny_genotype
    package_cell_phylogeny$phylogeny_birthtime <- phylogeny_birthtime
    package_cell_phylogeny$phylogeny_deathtime <- phylogeny_deathtime
    package_cell_phylogeny$phylogeny_order <- phylogeny_order

    package_clone_phylogeny <- list()
    package_clone_phylogeny$clone_phylogeny_labels <- clone_phylogeny_labels
    package_clone_phylogeny$clone_phylogeny_ID <- clone_phylogeny_ID
    package_clone_phylogeny$clone_phylogeny_origin <- clone_phylogeny_origin
    package_clone_phylogeny$clone_phylogeny_genotype <- clone_phylogeny_genotype
    package_clone_phylogeny$clone_phylogeny_birthtime <- clone_phylogeny_birthtime
    package_clone_phylogeny$clone_phylogeny_deathtime <- clone_phylogeny_deathtime
    package_clone_phylogeny$clone_hclust_nodes <- clone_hclust_nodes
    package_clone_phylogeny$clone_hclust_merge <- clone_hclust_merge
    package_clone_phylogeny$clone_hclust_height <- clone_hclust_height

    output <- list()
    output$phylogeny_clustering_truth <- phylogeny_clustering_truth
    output$cell_phylogeny_hclust <- phylogeny_hclust
    output$clone_phylogeny_hclust <- clone_phylogeny_hclust
    output$package_cell_phylogeny_hclust_extra <- package_cell_phylogeny_hclust_extra
    output$clone_phylogeny_phylo <- clone_phylogeny_phylo
    output$package_cell_phylogeny <- package_cell_phylogeny
    output$package_clone_phylogeny <- package_clone_phylogeny

    return(output)
}
