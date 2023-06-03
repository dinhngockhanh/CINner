# =========================================PLOT AVERAGE PLOIDY OVER TIME
#' @export
plot_average_ploidy <- function(model = "",
                                n_simulations = 0,
                                folder_workplace = "",
                                folder_plots = "",
                                vec_time = NULL,
                                unit_time = "year",
                                width = 1000,
                                height = 500,
                                compute_parallel = TRUE,
                                n_cores = NULL,
                                R_libPaths = NULL) {
    library(ggplot2)
    if (compute_parallel == FALSE) {
        #-------------------------Plot average ploidy in sequential mode
        for (iteration in 1:n_simulations) {
            plot_average_ploidy_one_simulation(
                model,
                iteration,
                folder_workplace,
                folder_plots,
                vec_time,
                unit_time,
                width,
                height
            )
        }
    } else {
        #---------------------------Plot average ploidy in parallel mode
        library(parallel)
        library(pbapply)
        #   Start parallel cluster
        if (is.null(n_cores)) {
            numCores <- detectCores()
        } else {
            numCores <- n_cores
        }
        cl <- makePSOCKcluster(numCores - 1)
        if (is.null(R_libPaths) == FALSE) {
            R_libPaths <<- R_libPaths
            clusterExport(cl, varlist = c(
                "R_libPaths"
            ))
            clusterEvalQ(cl = cl, .libPaths(R_libPaths))
        }
        clusterEvalQ(cl = cl, library(RColorBrewer))
        clusterEvalQ(cl = cl, library(fishplot))
        clusterEvalQ(cl = cl, library(wesanderson))
        clusterEvalQ(cl = cl, library(ggplot2))
        #   Prepare input parameters for plotting
        model <<- model
        folder_workplace <<- folder_workplace
        folder_plots <<- folder_plots
        width <<- width
        height <<- height
        vec_time <<- vec_time
        unit_time <<- unit_time
        find_ploidy <<- find_ploidy
        clusterExport(cl, varlist = c(
            "plot_average_ploidy_one_simulation",
            "find_ploidy",
            "model",
            "folder_workplace",
            "folder_plots",
            "vec_time",
            "unit_time",
            "width",
            "height"
        ))
        #   Plot in parallel
        pblapply(cl = cl, X = 1:n_simulations, FUN = function(iteration) {
            plot_average_ploidy_one_simulation(
                model,
                iteration,
                folder_workplace,
                folder_plots,
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

#' @export
plot_average_ploidy_one_simulation <- function(model,
                                               iteration,
                                               folder_workplace,
                                               folder_plots,
                                               vec_time,
                                               unit_time,
                                               width,
                                               height) {
    #----------------------------------------------Input simulation file
    filename <- paste(folder_workplace, model, "_simulation_", iteration, ".rda", sep = "")
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
    #-----------------------------------------Input the clonal evolution
    evolution_traj_time <- simulation$clonal_evolution$evolution_traj_time
    evolution_traj_clonal_ID <- simulation$clonal_evolution$evolution_traj_clonal_ID
    evolution_traj_population <- simulation$clonal_evolution$evolution_traj_population
    genotype_list_ploidy_chrom <- simulation$clonal_evolution$genotype_list_ploidy_chrom
    genotype_list_ploidy_block <- simulation$clonal_evolution$genotype_list_ploidy_block
    #-----------------------Find average ploidy dynamics from simulation
    N_chromosomes <- length(genotype_list_ploidy_chrom[[1]])
    vec_CN_block_no <- rep(0, N_chromosomes)
    for (chrom in 1:N_chromosomes) {
        vec_CN_block_no[chrom] <- length(genotype_list_ploidy_block[[1]][[chrom]][[1]])
    }
    #   Find average ploidy at desired time points
    library_clone <- c()
    library_ploidy <- c()
    vec_average_ploidy <- rep(0, length(vec_time))
    df_distribution_ploidy <- data.frame(matrix(0, nrow = 0, ncol = 3))
    colnames(df_distribution_ploidy) <- c("time", "ploidy", "percentage")
    for (j in 1:length(vec_time)) {
        time_plot <- vec_time_simulation[j]
        loc <- which.min(abs(evolution_traj_time - time_plot))
        clonal_ID <- evolution_traj_clonal_ID[[loc]]
        clonal_population <- evolution_traj_population[[loc]]
        clonal_ploidy <- rep(0, length(clonal_ID))
        for (i_clone in 1:length(clonal_ID)) {
            clone <- clonal_ID[i_clone]
            if (clone %in% library_clone) {
                ploidy <- library_ploidy[which(library_clone == clone)]
            } else {
                ploidy <- find_ploidy(
                    genotype_list_ploidy_chrom[[clone]],
                    genotype_list_ploidy_block[[clone]],
                    N_chromosomes, vec_CN_block_no
                )
                library_clone <- c(library_clone, clone)
                library_ploidy <- c(library_ploidy, ploidy)
            }
            clonal_ploidy[i_clone] <- find_ploidy(
                genotype_list_ploidy_chrom[[clone]],
                genotype_list_ploidy_block[[clone]],
                N_chromosomes, vec_CN_block_no
            )
        }
        vec_average_ploidy[j] <- sum(clonal_ploidy * clonal_population) / sum(clonal_population)
        clonal_ploidy_rounded <- round(clonal_ploidy)
        for (ploidy in 1:4) {
            ploidy_perc <- 100 * sum(clonal_population[which(clonal_ploidy_rounded == ploidy)]) / sum(clonal_population)
            df_distribution_ploidy[nrow(df_distribution_ploidy) + 1, ] <- c(time_plot, ploidy, ploidy_perc)
        }
    }
    #---------------------------------------Plot average ploidy dynamics
    df <- data.frame(vec_time, vec_average_ploidy)
    filename <- paste(folder_plots, model, "_sim", iteration, "_ploidy_average", ".jpeg", sep = "")
    jpeg(file = filename, width = width, height = height)
    p <- ggplot(df, aes(vec_time)) +
        annotate("rect", fill = "blue3", alpha = 0.5, xmin = vec_time[1], xmax = vec_time[length(vec_time)], ymin = 0.5, ymax = 1.5) +
        annotate("rect", fill = "azure3", alpha = 0.5, xmin = vec_time[1], xmax = vec_time[length(vec_time)], ymin = 1.5, ymax = 2.5) +
        annotate("rect", fill = "forestgreen", alpha = 0.5, xmin = vec_time[1], xmax = vec_time[length(vec_time)], ymin = 2.5, ymax = 3.5) +
        annotate("rect", fill = "firebrick2", alpha = 0.5, xmin = vec_time[1], xmax = vec_time[length(vec_time)], ymin = 3.5, ymax = 4.5) +
        geom_line(aes(y = vec_average_ploidy), colour = "black", size = 2) +
        xlim(vec_time[1], vec_time[length(vec_time)]) +
        scale_x_continuous(limits = c(vec_time[1], vec_time[length(vec_time)]), expand = c(0, 0)) +
        scale_y_continuous(limits = c(0, 5), expand = c(0, 0)) +
        xlab("Age") +
        ylab("Average ploidy") +
        theme(panel.background = element_rect(fill = "white", colour = "grey50")) +
        theme(text = element_text(size = 20))
    print(p)
    dev.off()
    #----------------------------------Plot ploidy distribution dynamics
    df_distribution_ploidy$ploidy <- paste0(df_distribution_ploidy$ploidy)
    filename <- paste(folder_plots, model, "_sim", iteration, "_ploidy_distribution", ".jpeg", sep = "")
    jpeg(file = filename, width = width, height = height)
    p <- ggplot(df_distribution_ploidy, aes(x = time, y = percentage, fill = ploidy)) +
        geom_area() +
        scale_x_continuous(expand = c(0, 0)) +
        scale_y_continuous(expand = c(0, 0)) +
        xlab("Age") +
        ylab("Percentage") +
        scale_fill_manual(values = c("blue3", "azure3", "forestgreen", "firebrick2")) +
        theme(panel.background = element_rect(fill = "white", colour = "grey50"), text = element_text(size = 20), legend.position = "top", legend.justification = "left", legend.direction = "horizontal", legend.key.width = unit(2.5, "cm"))
    print(p)
    dev.off()



    print(df_distribution_ploidy)
}

find_ploidy <- function(ploidy_chrom, ploidy_block,
                        N_chromosomes, vec_CN_block_no) {
    vec_CN_all <- c()
    for (chrom in 1:N_chromosomes) {
        vec_CN <- rep(0, vec_CN_block_no[chrom])
        no_strands <- ploidy_chrom[chrom]
        if (no_strands > 0) {
            for (strand in 1:no_strands) {
                vec_CN <- vec_CN + ploidy_block[[chrom]][[strand]]
            }
        }
        vec_CN_all <- c(vec_CN_all, vec_CN)
    }
    mean_ploidy <- mean(vec_CN_all)
    return(mean_ploidy)
}
