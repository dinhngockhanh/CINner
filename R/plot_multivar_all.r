plot_multivar_all <- function(model_prefix,
                              var1_name, var1_vals,
                              var2_name, var2_vals,
                              n_simulations_per_batch, stage_final) {
    #-----------------------------------------Input matrix of statistics
    load(file = paste(model_prefix, "_statistics.rda", sep = ""))
    #----------------------------Create data frame for plotting heatmaps
    if (stage_final >= 1) {
        #   Get statistics from each simulation
        mat_true_n_clones_whole_pop <- array(0, c(length(var1_vals), length(var2_vals), n_simulations_per_batch))
        for (row in 1:length(var1_vals)) {
            for (col in 1:length(var2_vals)) {
                batch_statistics <- mat_simulation_statistics[[row]][[col]]
                for (sim in 1:n_simulations_per_batch) {
                    true_n_clones_whole_pop_values <- batch_statistics[[sim]]$true_n_clones_whole_pop_values
                    mat_true_n_clones_whole_pop[row, col, sim] <- true_n_clones_whole_pop_values[length(true_n_clones_whole_pop_values)]
                }
            }
        }
        #   Create data frame for all statistics
        dataframe_plot <- data.frame(matrix(nrow = 0, ncol = 3))
        for (row in 1:length(var1_vals)) {
            for (col in 1:length(var2_vals)) {
                dataframe_plot[nrow(dataframe_plot) + 1, ] <-
                    c(
                        var1_vals[row], var2_vals[col],
                        mean(mat_true_n_clones_whole_pop[row, col, ])
                    )
            }
        }
        colnames(dataframe_plot) <- c(
            "x", "y",
            "true_n_clones_whole_pop_mean"
        )
        #   Order on x- and y- axes
        dataframe_plot$x <- as.factor(dataframe_plot$x)
        dataframe_plot$x <- factor(dataframe_plot$x, levels = var1_vals)
        dataframe_plot$y <- as.factor(dataframe_plot$y)
        dataframe_plot$y <- factor(dataframe_plot$y, levels = var2_vals)
    }
    if (stage_final >= 2) {
        #   Get statistics from each simulation
        mat_true_n_clones_sample <- array(0, c(length(var1_vals), length(var2_vals), n_simulations_per_batch))
        mat_cn_based_n_clones_sample <- array(0, c(length(var1_vals), length(var2_vals), n_simulations_per_batch))
        mat_ratio_n_clones_sample <- array(0, c(length(var1_vals), length(var2_vals), n_simulations_per_batch))
        mat_diploid_perc_sample <- array(0, c(length(var1_vals), length(var2_vals), n_simulations_per_batch))
        mat_driver_perc_sample <- array(0, c(length(var1_vals), length(var2_vals), n_simulations_per_batch))
        mat_WGD_perc_sample <- array(0, c(length(var1_vals), length(var2_vals), n_simulations_per_batch))
        mat_missegregation_perc_sample <- array(0, c(length(var1_vals), length(var2_vals), n_simulations_per_batch))

        mat_miss_before_WGD <- array(0, c(length(var1_vals), length(var2_vals), n_simulations_per_batch))
        mat_miss_after_WGD <- array(0, c(length(var1_vals), length(var2_vals), n_simulations_per_batch))

        for (row in 1:length(var1_vals)) {
            for (col in 1:length(var2_vals)) {
                batch_statistics <- mat_simulation_statistics[[row]][[col]]
                for (sim in 1:n_simulations_per_batch) {
                    true_n_clones_whole_pop_values <- batch_statistics[[sim]]$true_n_clones_whole_pop_values
                    mat_true_n_clones_sample[row, col, sim] <- batch_statistics[[sim]]$stats_sample$true_n_clones[1]
                    mat_cn_based_n_clones_sample[row, col, sim] <- batch_statistics[[sim]]$stats_sample$total_cn_based_n_clones[1]
                    mat_ratio_n_clones_sample[row, col, sim] <- 100 * batch_statistics[[sim]]$stats_sample$total_cn_based_n_clones[1] / batch_statistics[[sim]]$stats_sample$true_n_clones[1]
                    mat_diploid_perc_sample[row, col, sim] <- batch_statistics[[sim]]$stats_sample$perc_initial
                    mat_driver_perc_sample[row, col, sim] <- batch_statistics[[sim]]$stats_sample$perc_driver_mut
                    mat_WGD_perc_sample[row, col, sim] <- batch_statistics[[sim]]$stats_sample$perc_WGD
                    mat_missegregation_perc_sample[row, col, sim] <- batch_statistics[[sim]]$stats_sample$perc_missegregation

                    mat_miss_before_WGD[row, col, sim] <- batch_statistics[[sim]]$stats_sample$mean_miss_before_WGD
                    mat_miss_after_WGD[row, col, sim] <- batch_statistics[[sim]]$stats_sample$mean_miss_after_WGD
                }
            }
        }
        #   Update data frame for all statistics
        ind <- 0
        dataframe_plot$true_n_clones_sample_mean <- 0
        dataframe_plot$cn_based_n_clones_sample_mean <- 0
        dataframe_plot$ratio_n_clones_sample_mean <- 0
        dataframe_plot$diploid_perc_sample_mean <- 0
        dataframe_plot$driver_perc_sample_mean <- 0
        dataframe_plot$WGD_perc_sample_mean <- 0
        dataframe_plot$missegregation_perc_sample_mean <- 0

        dataframe_plot$count_missegregations_before_WGD_mean <- 0
        dataframe_plot$count_missegregations_after_WGD_mean <- 0

        for (row in 1:length(var1_vals)) {
            for (col in 1:length(var2_vals)) {
                ind <- ind + 1
                dataframe_plot$true_n_clones_sample_mean[ind] <- mean(mat_true_n_clones_sample[row, col, ])
                dataframe_plot$cn_based_n_clones_sample_mean[ind] <- mean(mat_cn_based_n_clones_sample[row, col, ])
                dataframe_plot$ratio_n_clones_sample_mean[ind] <- 100 * mean(mat_cn_based_n_clones_sample[row, col, ]) / mean(mat_true_n_clones_sample[row, col, ])
                dataframe_plot$diploid_perc_sample_mean[ind] <- mean(mat_diploid_perc_sample[row, col, ])
                dataframe_plot$driver_perc_sample_mean[ind] <- mean(mat_driver_perc_sample[row, col, ])
                dataframe_plot$WGD_perc_sample_mean[ind] <- mean(mat_WGD_perc_sample[row, col, ])
                dataframe_plot$missegregation_perc_sample_mean[ind] <- mean(mat_missegregation_perc_sample[row, col, ])


                mat_miss_before_WGD_ind <- mat_miss_before_WGD[row, col, ]
                mat_miss_before_WGD_ind <- mat_miss_before_WGD_ind[is.na(mat_miss_before_WGD_ind) == FALSE]
                dataframe_plot$count_missegregations_before_WGD_mean[ind] <- mean(mat_miss_before_WGD_ind)
                mat_miss_after_WGD_ind <- mat_miss_after_WGD[row, col, ]
                mat_miss_after_WGD_ind <- mat_miss_after_WGD_ind[is.na(mat_miss_after_WGD_ind) == FALSE]
                dataframe_plot$count_missegregations_after_WGD_mean[ind] <- mean(mat_miss_after_WGD_ind)
            }
        }
        #   Order on x- and y- axes
        dataframe_plot$x <- as.factor(dataframe_plot$x)
        dataframe_plot$x <- factor(dataframe_plot$x, levels = var1_vals)
        dataframe_plot$y <- as.factor(dataframe_plot$y)
        dataframe_plot$y <- factor(dataframe_plot$y, levels = var2_vals)
    }
    #-------------------------Plot true clone counts in whole population
    if (stage_final >= 1) {
        filename <- paste(model_prefix, "_pop_Nclones_true.jpeg", sep = "")
        jpeg(file = filename, width = 1200, height = 1000)
        #   Plot heatmap
        p <- ggplot(dataframe_plot, aes(x, y, fill = true_n_clones_whole_pop_mean)) +
            geom_tile() +
            geom_text(aes(label = round(true_n_clones_whole_pop_mean, 1))) +
            scale_fill_distiller(palette = "RdPu") +
            theme_ipsum() +
            labs(title = "True clone count in whole population", x = var1_name, y = var2_name, fill = "Clone count") +
            theme_classic(
                base_size = 20
            )
        #   Plot averages
        h_average <- dataframe_plot %>%
            group_by(y) %>%
            summarise(true_n_clones_whole_pop_mean = sum(true_n_clones_whole_pop_mean) / length(var2_vals)) %>%
            mutate(x = "Average")
        v_average <- dataframe_plot %>%
            group_by(x) %>%
            summarise(true_n_clones_whole_pop_mean = sum(true_n_clones_whole_pop_mean) / length(var1_vals)) %>%
            mutate(y = "Average")
        p <- p +
            # geom_point(
            #     data = h_average,
            #     aes(color = true_n_clones_whole_pop_mean),
            #     size = 10,
            #     shape = 19
            # ) +
            # geom_point(
            #     data = v_average,
            #     aes(color = true_n_clones_whole_pop_mean),
            #     size = 10,
            #     shape = 19
            # ) +
            geom_text(data = h_average, size = 10, aes(label = round(true_n_clones_whole_pop_mean, 1))) +
            geom_text(data = v_average, size = 10, aes(label = round(true_n_clones_whole_pop_mean, 1))) +
            theme(
                panel.grid.major.x = element_blank(),
                panel.grid.minor.x = element_blank(),
                panel.grid.major.y = element_blank(),
                panel.grid.minor.y = element_blank(),
                panel.background = element_rect(fill = "white"),
                axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 0.5, size = 20, face = NULL),
                axis.text.y = element_text(size = 20, face = NULL),
                plot.title = element_text(size = 20, face = "bold"),
                legend.title = element_text(face = "bold", size = 20)
            )
        print(p)
        dev.off()
    }
    #-----------------------------------Plot true clone counts in sample
    if (stage_final >= 2) {
        filename <- paste(model_prefix, "_sample_Nclones_true.jpeg", sep = "")
        jpeg(file = filename, width = 1200, height = 1000)
        #   Plot heatmap
        p <- ggplot(dataframe_plot, aes(x, y, fill = true_n_clones_sample_mean)) +
            geom_tile() +
            geom_text(aes(label = round(true_n_clones_sample_mean, 1))) +
            scale_fill_distiller(palette = "RdPu") +
            theme_ipsum() +
            labs(title = "True clone count in sampled cells", x = var1_name, y = var2_name, fill = "Clone count") +
            theme_classic(
                base_size = 20
            )
        #   Plot averages
        h_average <- dataframe_plot %>%
            group_by(y) %>%
            summarise(true_n_clones_sample_mean = sum(true_n_clones_sample_mean) / length(var2_vals)) %>%
            mutate(x = "Average")
        v_average <- dataframe_plot %>%
            group_by(x) %>%
            summarise(true_n_clones_sample_mean = sum(true_n_clones_sample_mean) / length(var1_vals)) %>%
            mutate(y = "Average")
        p <- p +
            # geom_point(
            #     data = h_average,
            #     aes(color = true_n_clones_sample_mean),
            #     size = 10,
            #     shape = 19
            # ) +
            # geom_point(
            #     data = v_average,
            #     aes(color = true_n_clones_sample_mean),
            #     size = 10,
            #     shape = 19
            # ) +
            geom_text(data = h_average, size = 10, aes(label = round(true_n_clones_sample_mean, 1))) +
            geom_text(data = v_average, size = 10, aes(label = round(true_n_clones_sample_mean, 1))) +
            theme(
                panel.grid.major.x = element_blank(),
                panel.grid.minor.x = element_blank(),
                panel.grid.major.y = element_blank(),
                panel.grid.minor.y = element_blank(),
                panel.background = element_rect(fill = "white"),
                axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 0.5, size = 20, face = NULL),
                axis.text.y = element_text(size = 20, face = NULL),
                plot.title = element_text(size = 20, face = "bold"),
                legend.title = element_text(face = "bold", size = 20)
            )
        print(p)
        dev.off()
    }
    #-------------------------Plot total-CN-based clone counts in sample
    if (stage_final >= 2) {
        filename <- paste(model_prefix, "_sample_Nclones_cn_based.jpeg", sep = "")
        jpeg(file = filename, width = 1200, height = 1000)
        #   Plot heatmap
        p <- ggplot(dataframe_plot, aes(x, y, fill = cn_based_n_clones_sample_mean)) +
            geom_tile() +
            geom_text(aes(label = round(cn_based_n_clones_sample_mean, 1))) +
            # scale_fill_viridis(discrete = FALSE) +
            scale_fill_distiller(palette = "RdPu") +
            theme_ipsum() +
            labs(title = "Total-CN-based clone count in sampled cells", x = var1_name, y = var2_name, fill = "Clone count") +
            theme_classic(
                base_size = 20
            )
        #   Plot averages
        h_average <- dataframe_plot %>%
            group_by(y) %>%
            summarise(cn_based_n_clones_sample_mean = sum(cn_based_n_clones_sample_mean) / length(var2_vals)) %>%
            mutate(x = "Average")
        v_average <- dataframe_plot %>%
            group_by(x) %>%
            summarise(cn_based_n_clones_sample_mean = sum(cn_based_n_clones_sample_mean) / length(var1_vals)) %>%
            mutate(y = "Average")
        p <- p +
            # geom_point(
            #     data = h_average,
            #     aes(color = cn_based_n_clones_sample_mean),
            #     size = 10,
            #     shape = 19
            # ) +
            # geom_point(
            #     data = v_average,
            #     aes(color = cn_based_n_clones_sample_mean),
            #     size = 10,
            #     shape = 19
            # ) +
            geom_text(data = h_average, size = 10, aes(label = round(cn_based_n_clones_sample_mean, 1))) +
            geom_text(data = v_average, size = 10, aes(label = round(cn_based_n_clones_sample_mean, 1))) +
            theme(
                panel.grid.major.x = element_blank(),
                panel.grid.minor.x = element_blank(),
                panel.grid.major.y = element_blank(),
                panel.grid.minor.y = element_blank(),
                panel.background = element_rect(fill = "white"),
                axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 0.5, size = 20, face = NULL),
                axis.text.y = element_text(size = 20, face = NULL),
                plot.title = element_text(size = 20, face = "bold"),
                legend.title = element_text(face = "bold", size = 20)
            )
        print(p)
        dev.off()
    }
    #--------Plot ratio of total-CN-based to true clone counts in sample
    if (stage_final >= 2) {
        filename <- paste(model_prefix, "_sample_Nclones_ratio_cn_based_over_true.jpeg", sep = "")
        jpeg(file = filename, width = 1200, height = 1000)
        #   Plot heatmap
        p <- ggplot(dataframe_plot, aes(x, y, fill = ratio_n_clones_sample_mean)) +
            geom_tile() +
            geom_text(aes(label = round(ratio_n_clones_sample_mean, 1))) +
            # scale_fill_viridis(discrete = FALSE) +
            scale_fill_distiller(palette = "RdPu") +
            theme_ipsum() +
            labs(title = "Ratio of total-CN-based clone count over true clone count in sampled cells", x = var1_name, y = var2_name, fill = "%") +
            theme_classic(
                base_size = 20
            )
        #   Plot averages
        h_average <- dataframe_plot %>%
            group_by(y) %>%
            summarise(ratio_n_clones_sample_mean = sum(ratio_n_clones_sample_mean) / length(var2_vals)) %>%
            mutate(x = "Average")
        v_average <- dataframe_plot %>%
            group_by(x) %>%
            summarise(ratio_n_clones_sample_mean = sum(ratio_n_clones_sample_mean) / length(var1_vals)) %>%
            mutate(y = "Average")
        p <- p +
            # geom_point(
            #     data = h_average,
            #     aes(color = ratio_n_clones_sample_mean),
            #     size = 10,
            #     shape = 19
            # ) +
            # geom_point(
            #     data = v_average,
            #     aes(color = ratio_n_clones_sample_mean),
            #     size = 10,
            #     shape = 19
            # ) +
            geom_text(data = h_average, size = 10, aes(label = round(ratio_n_clones_sample_mean, 1))) +
            geom_text(data = v_average, size = 10, aes(label = round(ratio_n_clones_sample_mean, 1))) +
            theme(
                panel.grid.major.x = element_blank(),
                panel.grid.minor.x = element_blank(),
                panel.grid.major.y = element_blank(),
                panel.grid.minor.y = element_blank(),
                panel.background = element_rect(fill = "white"),
                axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 0.5, size = 20, face = NULL),
                axis.text.y = element_text(size = 20, face = NULL),
                plot.title = element_text(size = 20, face = "bold"),
                legend.title = element_text(face = "bold", size = 20)
            )
        print(p)
        dev.off()
    }
    #----------------------------------Plot diploid percentage in sample
    if (stage_final >= 2) {
        filename <- paste(model_prefix, "_sample_perc_diploid.jpeg", sep = "")
        jpeg(file = filename, width = 1200, height = 1000)
        #   Plot heatmap
        p <- ggplot(dataframe_plot, aes(x, y, fill = diploid_perc_sample_mean)) +
            geom_tile() +
            geom_text(aes(label = round(diploid_perc_sample_mean, 1))) +
            # scale_fill_viridis(discrete = FALSE) +
            scale_fill_distiller(palette = "RdPu") +
            theme_ipsum() +
            labs(title = "Percentage of sampled cells that are diploid", x = var1_name, y = var2_name, fill = "%") +
            theme_classic(
                base_size = 20
            )
        #   Plot averages
        h_average <- dataframe_plot %>%
            group_by(y) %>%
            summarise(diploid_perc_sample_mean = sum(diploid_perc_sample_mean) / length(var2_vals)) %>%
            mutate(x = "Average")
        v_average <- dataframe_plot %>%
            group_by(x) %>%
            summarise(diploid_perc_sample_mean = sum(diploid_perc_sample_mean) / length(var1_vals)) %>%
            mutate(y = "Average")
        p <- p +
            # geom_point(
            #     data = h_average,
            #     aes(color = diploid_perc_sample_mean),
            #     size = 10,
            #     shape = 19
            # ) +
            # geom_point(
            #     data = v_average,
            #     aes(color = diploid_perc_sample_mean),
            #     size = 10,
            #     shape = 19
            # ) +
            geom_text(data = h_average, size = 10, aes(label = round(diploid_perc_sample_mean, 1))) +
            geom_text(data = v_average, size = 10, aes(label = round(diploid_perc_sample_mean, 1))) +
            theme(
                panel.grid.major.x = element_blank(),
                panel.grid.minor.x = element_blank(),
                panel.grid.major.y = element_blank(),
                panel.grid.minor.y = element_blank(),
                panel.background = element_rect(fill = "white"),
                axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 0.5, size = 20, face = NULL),
                axis.text.y = element_text(size = 20, face = NULL),
                plot.title = element_text(size = 20, face = "bold"),
                legend.title = element_text(face = "bold", size = 20)
            )
        print(p)
        dev.off()
    }
    #-----------------------------------Plot driver percentage in sample
    if (stage_final >= 2) {
        filename <- paste(model_prefix, "_sample_perc_driver.jpeg", sep = "")
        jpeg(file = filename, width = 1200, height = 1000)
        #   Plot heatmap
        p <- ggplot(dataframe_plot, aes(x, y, fill = driver_perc_sample_mean)) +
            geom_tile() +
            geom_text(aes(label = round(driver_perc_sample_mean, 1))) +
            # scale_fill_viridis(discrete = FALSE) +
            scale_fill_distiller(palette = "RdPu") +
            theme_ipsum() +
            labs(title = "Percentage of sampled cells that contain driver mutations", x = var1_name, y = var2_name, fill = "%") +
            theme_classic(
                base_size = 20
            )
        #   Plot averages
        h_average <- dataframe_plot %>%
            group_by(y) %>%
            summarise(driver_perc_sample_mean = sum(driver_perc_sample_mean) / length(var2_vals)) %>%
            mutate(x = "Average")
        v_average <- dataframe_plot %>%
            group_by(x) %>%
            summarise(driver_perc_sample_mean = sum(driver_perc_sample_mean) / length(var1_vals)) %>%
            mutate(y = "Average")
        p <- p +
            # geom_point(
            #     data = h_average,
            #     aes(color = driver_perc_sample_mean),
            #     size = 10,
            #     shape = 19
            # ) +
            # geom_point(
            #     data = v_average,
            #     aes(color = driver_perc_sample_mean),
            #     size = 10,
            #     shape = 19
            # ) +
            geom_text(data = h_average, size = 10, aes(label = round(driver_perc_sample_mean, 1))) +
            geom_text(data = v_average, size = 10, aes(label = round(driver_perc_sample_mean, 1))) +
            theme(
                panel.grid.major.x = element_blank(),
                panel.grid.minor.x = element_blank(),
                panel.grid.major.y = element_blank(),
                panel.grid.minor.y = element_blank(),
                panel.background = element_rect(fill = "white"),
                axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 0.5, size = 20, face = NULL),
                axis.text.y = element_text(size = 20, face = NULL),
                plot.title = element_text(size = 20, face = "bold"),
                legend.title = element_text(face = "bold", size = 20)
            )
        print(p)
        dev.off()
    }
    #--------------------------------------Plot WGD percentage in sample
    if (stage_final >= 2) {
        filename <- paste(model_prefix, "_sample_perc_WGD.jpeg", sep = "")
        jpeg(file = filename, width = 1200, height = 1000)
        #   Plot heatmap
        p <- ggplot(dataframe_plot, aes(x, y, fill = WGD_perc_sample_mean)) +
            geom_tile() +
            geom_text(aes(label = round(WGD_perc_sample_mean, 1))) +
            # scale_fill_viridis(discrete = FALSE) +
            scale_fill_distiller(palette = "RdPu") +
            theme_ipsum() +
            labs(title = "Percentage of sampled cells that underwent WGD", x = var1_name, y = var2_name, fill = "%") +
            theme_classic(
                base_size = 20
            )
        #   Plot averages
        h_average <- dataframe_plot %>%
            group_by(y) %>%
            summarise(WGD_perc_sample_mean = sum(WGD_perc_sample_mean) / length(var2_vals)) %>%
            mutate(x = "Average")
        v_average <- dataframe_plot %>%
            group_by(x) %>%
            summarise(WGD_perc_sample_mean = sum(WGD_perc_sample_mean) / length(var1_vals)) %>%
            mutate(y = "Average")
        p <- p +
            # geom_point(
            #     data = h_average,
            #     aes(color = WGD_perc_sample_mean),
            #     size = 10,
            #     shape = 19
            # ) +
            # geom_point(
            #     data = v_average,
            #     aes(color = WGD_perc_sample_mean),
            #     size = 10,
            #     shape = 19
            # ) +
            geom_text(data = h_average, size = 10, aes(label = round(WGD_perc_sample_mean, 1))) +
            geom_text(data = v_average, size = 10, aes(label = round(WGD_perc_sample_mean, 1))) +
            theme(
                panel.grid.major.x = element_blank(),
                panel.grid.minor.x = element_blank(),
                panel.grid.major.y = element_blank(),
                panel.grid.minor.y = element_blank(),
                panel.background = element_rect(fill = "white"),
                axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 0.5, size = 20, face = NULL),
                axis.text.y = element_text(size = 20, face = NULL),
                plot.title = element_text(size = 20, face = "bold"),
                legend.title = element_text(face = "bold", size = 20)
            )
        print(p)
        dev.off()
    }
    #---------------------------Plot missegregation percentage in sample
    if (stage_final >= 2) {
        filename <- paste(model_prefix, "_sample_perc_missegregation.jpeg", sep = "")
        jpeg(file = filename, width = 1200, height = 1000)
        #   Plot heatmap
        p <- ggplot(dataframe_plot, aes(x, y, fill = missegregation_perc_sample_mean)) +
            geom_tile() +
            geom_text(aes(label = round(missegregation_perc_sample_mean, 1))) +
            # scale_fill_viridis(discrete = FALSE) +
            scale_fill_distiller(palette = "RdPu") +
            theme_ipsum() +
            labs(title = "Percentage of sampled cells that underwent missegregations", x = var1_name, y = var2_name, fill = "%") +
            theme_classic(
                base_size = 20
            )
        #   Plot averages
        h_average <- dataframe_plot %>%
            group_by(y) %>%
            summarise(missegregation_perc_sample_mean = sum(missegregation_perc_sample_mean) / length(var2_vals)) %>%
            mutate(x = "Average")
        v_average <- dataframe_plot %>%
            group_by(x) %>%
            summarise(missegregation_perc_sample_mean = sum(missegregation_perc_sample_mean) / length(var1_vals)) %>%
            mutate(y = "Average")
        p <- p +
            # geom_point(
            #     data = h_average,
            #     aes(color = missegregation_perc_sample_mean),
            #     size = 10,
            #     shape = 19
            # ) +
            # geom_point(
            #     data = v_average,
            #     aes(color = missegregation_perc_sample_mean),
            #     size = 10,
            #     shape = 19
            # ) +
            geom_text(data = h_average, size = 10, aes(label = round(missegregation_perc_sample_mean, 1))) +
            geom_text(data = v_average, size = 10, aes(label = round(missegregation_perc_sample_mean, 1))) +
            theme(
                panel.grid.major.x = element_blank(),
                panel.grid.minor.x = element_blank(),
                panel.grid.major.y = element_blank(),
                panel.grid.minor.y = element_blank(),
                panel.background = element_rect(fill = "white"),
                axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 0.5, size = 20, face = NULL),
                axis.text.y = element_text(size = 20, face = NULL),
                plot.title = element_text(size = 20, face = "bold"),
                legend.title = element_text(face = "bold", size = 20)
            )
        print(p)
        dev.off()
    }
    #------------Plot mean count of missegregations before WGD in sample
    if (stage_final >= 2) {
        filename <- paste(model_prefix, "_count_missegregations_before_WGD.jpeg", sep = "")
        jpeg(file = filename, width = 1200, height = 1000)
        #   Plot heatmap
        p <- ggplot(dataframe_plot, aes(x, y, fill = count_missegregations_before_WGD_mean)) +
            geom_tile() +
            geom_text(aes(label = round(count_missegregations_before_WGD_mean, 1))) +
            # scale_fill_viridis(discrete = FALSE) +
            scale_fill_distiller(palette = "RdPu") +
            theme_ipsum() +
            labs(title = "Count of missegregations before WGD (average weighted by cell counts)", x = var1_name, y = var2_name, fill = "Count") +
            theme_classic(
                base_size = 20
            )
        #   Plot averages
        h_average <- dataframe_plot %>%
            group_by(y) %>%
            summarise(count_missegregations_before_WGD_mean = sum(count_missegregations_before_WGD_mean) / length(var2_vals)) %>%
            mutate(x = "Average")
        v_average <- dataframe_plot %>%
            group_by(x) %>%
            summarise(count_missegregations_before_WGD_mean = sum(count_missegregations_before_WGD_mean) / length(var1_vals)) %>%
            mutate(y = "Average")
        p <- p +
            # geom_point(
            #     data = h_average,
            #     aes(color = count_missegregations_before_WGD_mean),
            #     size = 10,
            #     shape = 19
            # ) +
            # geom_point(
            #     data = v_average,
            #     aes(color = count_missegregations_before_WGD_mean),
            #     size = 10,
            #     shape = 19
            # ) +
            geom_text(data = h_average, size = 10, aes(label = round(count_missegregations_before_WGD_mean, 2))) +
            geom_text(data = v_average, size = 10, aes(label = round(count_missegregations_before_WGD_mean, 2))) +
            theme(
                panel.grid.major.x = element_blank(),
                panel.grid.minor.x = element_blank(),
                panel.grid.major.y = element_blank(),
                panel.grid.minor.y = element_blank(),
                panel.background = element_rect(fill = "white"),
                axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 0.5, size = 20, face = NULL),
                axis.text.y = element_text(size = 20, face = NULL),
                plot.title = element_text(size = 20, face = "bold"),
                legend.title = element_text(face = "bold", size = 20)
            )
        print(p)
        dev.off()
    }



    #------------Plot mean count of missegregations after WGD in sample
    if (stage_final >= 2) {
        filename <- paste(model_prefix, "_count_missegregations_after_WGD.jpeg", sep = "")
        jpeg(file = filename, width = 1200, height = 1000)
        #   Plot heatmap
        p <- ggplot(dataframe_plot, aes(x, y, fill = count_missegregations_after_WGD_mean)) +
            geom_tile() +
            geom_text(aes(label = round(count_missegregations_after_WGD_mean, 1))) +
            # scale_fill_viridis(discrete = FALSE) +
            scale_fill_distiller(palette = "RdPu") +
            theme_ipsum() +
            labs(title = "Count of missegregations after WGD (average weighted by cell counts)", x = var1_name, y = var2_name, fill = "Count") +
            theme_classic(
                base_size = 20
            )
        #   Plot averages
        h_average <- dataframe_plot %>%
            group_by(y) %>%
            summarise(count_missegregations_after_WGD_mean = sum(count_missegregations_after_WGD_mean) / length(var2_vals)) %>%
            mutate(x = "Average")
        v_average <- dataframe_plot %>%
            group_by(x) %>%
            summarise(count_missegregations_after_WGD_mean = sum(count_missegregations_after_WGD_mean) / length(var1_vals)) %>%
            mutate(y = "Average")
        p <- p +
            # geom_point(
            #     data = h_average,
            #     aes(color = count_missegregations_after_WGD_mean),
            #     size = 10,
            #     shape = 19
            # ) +
            # geom_point(
            #     data = v_average,
            #     aes(color = count_missegregations_after_WGD_mean),
            #     size = 10,
            #     shape = 19
            # ) +
            geom_text(data = h_average, size = 10, aes(label = round(count_missegregations_after_WGD_mean, 2))) +
            geom_text(data = v_average, size = 10, aes(label = round(count_missegregations_after_WGD_mean, 2))) +
            theme(
                panel.grid.major.x = element_blank(),
                panel.grid.minor.x = element_blank(),
                panel.grid.major.y = element_blank(),
                panel.grid.minor.y = element_blank(),
                panel.background = element_rect(fill = "white"),
                axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 0.5, size = 20, face = NULL),
                axis.text.y = element_text(size = 20, face = NULL),
                plot.title = element_text(size = 20, face = "bold"),
                legend.title = element_text(face = "bold", size = 20)
            )
        print(p)
        dev.off()
    }
}
