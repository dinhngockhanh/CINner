analysis_dlp_cn <- function(CNbins, table_gc, filename_prefix) {
    CNbins_DLP <- CNbins
    #-----------------------------------Set parameters for fitting sigma
    #   Slope of GC-readcount linear relation
    gc_slope <- 1
    #   Intercept of GC-readcount linear relation
    gc_int <- 0
    #   Fitting is limited to bins of CN=fit_ploidy in cells of ploidy=fit_ploidy
    fit_ploidy <- 2
    #   Lower threshold of reads for bins to be included in fitting data
    min_reads <- 100
    #   Number of trials from prior distribution for ABC
    fit_ABC_count <- 100000
    #   Lower threshold for sigma in ABC prior distribution
    fit_ABC_sigma_min <- 0.001
    #   Upper threshold for sigma in ABC prior distribution
    fit_ABC_sigma_max <- 0.1
    #   Tolerance for accepting sigma in ABC-rejection
    fit_ABC_tol <- 0.01
    #   Number of bins to simulater per ABC trial
    fit_bin_count <- 1000
    #--------------------------Find distribution of total reads per cell
    cat("Find distribution of total reads in CNbins...\n")
    pb <- txtProgressBar(
        min = 0, max = length(list_cell_id),
        style = 3, width = 50, char = "="
    )
    list_cell_id <- unique(CNbins_DLP$cell_id)
    vec_total_reads <- rep(0, length = length(list_cell_id))
    for (i in 1:length(list_cell_id)) {
        setTxtProgressBar(pb, i)
        vec_total_reads[i] <- sum(CNbins_DLP$reads[CNbins_DLP$cell_id == list_cell_id[i]])
    }
    cat("\n")
    #-------------------------Supplement the CN data with GC/mappability
    cat("Add GC/mappability data to CNbins...\n")
    pb <- txtProgressBar(
        min = 0, max = nrow(table_gc),
        style = 3, width = 50, char = "="
    )
    CNbins_DLP$gc <- -1
    CNbins_DLP$map <- -1
    for (row in 1:nrow(table_gc)) {
        setTxtProgressBar(pb, row)
        vec_loc <- which(CNbins_DLP$chr == table_gc$chr[row] & CNbins_DLP$start_point == table_gc$start[row] & CNbins_DLP$end_point == table_gc$end[row])
        CNbins_DLP$gc[vec_loc] <- table_gc$gc[row]
        CNbins_DLP$map[vec_loc] <- table_gc$map[row]
    }
    vec_delete <- which(CNbins_DLP$gc <= 0 | CNbins_DLP$map <= 0 | CNbins_DLP$state <= 0)
    if (length(vec_delete) > 0) {
        CNbins_DLP <- CNbins_DLP[-vec_delete, ]
    }
    cat("\n")
    #-----------------------------------Compute cell ploidy = average CN
    list_cell_id <- unique(CNbins_DLP$cell_id)
    CNbins_DLP$ploidy <- 0
    cat("Compute cell ploidy...\n")
    pb <- txtProgressBar(
        min = 0, max = length(list_cell_id),
        style = 3, width = 50, char = "="
    )
    for (j in 1:length(list_cell_id)) {
        setTxtProgressBar(pb, j)
        vec_loc <- which(CNbins_DLP$cell_id == list_cell_id[j])
        CNbins_DLP$ploidy[vec_loc] <- round(mean(CNbins_DLP$state[vec_loc]))
    }
    cat("\n")
    CNbins_DLP$multiplier <- CNbins_DLP$state / CNbins_DLP$ploidy
    #--------------------------------------------Filter data for fitting
    #   Filters:    bins of CN=2 in cells of ploidy=2
    #               readcounts >= 100
    CNbins_fit <- CNbins_DLP[which(CNbins_DLP$state == fit_ploidy & CNbins_DLP$ploidy == fit_ploidy & CNbins_DLP$reads >= min_reads), ]
    CNbins_fit$reads_normalized <- 0
    #   Normalize reads for GC and total reads per cell
    list_cell_id <- unique(CNbins_fit$cell_id)
    for (i in 1:length(list_cell_id)) {
        vec_loc <- which(CNbins_fit$cell_id == list_cell_id[i])
        #   Normalize reads for GC
        CNbins_fit$reads_normalized[vec_loc] <- CNbins_fit$reads[vec_loc] / (fit_ploidy * (gc_slope * CNbins_fit$gc[vec_loc] + gc_int))
        #   Normalize reads for total reads per cell
        CNbins_fit$reads_normalized[vec_loc] <- CNbins_fit$reads_normalized[vec_loc] / mean(CNbins_fit$reads_normalized[vec_loc])
    }
    #---------------------------------------------------Fit sigma by ABC
    #   Find total reads for bins in simulations
    num_reads <- fit_bin_count * sum(CNbins_fit$reads) / nrow(CNbins_fit)
    #   Find target variance from normalized DLP data
    var_DLP <- var(CNbins_fit$reads_normalized)
    #   Define objective function for ABC fitting
    func_ABC <- function(fit_sigma) {
        #   Use Adam's noisy readcount model
        sims_gc <- 0.5
        sims_observed_CN <- fit_ploidy * (gc_slope * sims_gc + gc_int)
        sims_noisy_CN <- rgamma(n = fit_bin_count, shape = sims_observed_CN / fit_sigma, scale = fit_sigma)
        sims_noisy_CN_pval <- sims_noisy_CN / sum(sims_noisy_CN)
        sims_reads <- rmultinom(n = 1, size = num_reads, prob = sims_noisy_CN_pval)
        #   Normalize noisy readcount the same way that DLP data was
        sims_reads_normalized <- sims_reads / (fit_ploidy * (gc_slope * sims_gc + gc_int))
        sims_reads_normalized <- sims_reads_normalized / mean(sims_reads_normalized)
        #   Statistics = variance of normalized simulated readcounts
        stat <- var(sims_reads_normalized)
        return(stat)
    }
    #   Simulate sigma for ABC
    vec_sigma <- runif(fit_ABC_count, min = fit_ABC_sigma_min, max = fit_ABC_sigma_max)
    #   Make simulations for each sigma
    vec_stat <- sapply(vec_sigma, func_ABC)
    #   Perform ABC
    result <- abc(target = c(var_DLP), param = vec_sigma, sumstat = vec_stat, tol = fit_ABC_tol, method = "rejection")
    vec_sigma_posterior <- result$unadj.values
    mean_sigma <- mean(vec_sigma_posterior)
    #-----------------Create simulated readcounts for DLP for validation
    CNbins_DLP_validation <- CNbins_DLP[which(((CNbins_DLP$multiplier %% 0.5) == 0) & (CNbins_DLP$multiplier <= 1.5)), ]
    CNbins_DLP_validation$reads_simulated <- 0
    list_cell_id <- unique(CNbins_DLP_validation$cell_id)
    cat("Simulate readcounts for validation...\n")
    pb <- txtProgressBar(
        min = 0, max = length(list_cell_id),
        style = 3, width = 50, char = "="
    )
    for (i in 1:length(list_cell_id)) {
        setTxtProgressBar(pb, i)
        vec_loc <- which(CNbins_DLP_validation$cell_id == list_cell_id[i])
        num_reads <- sum(CNbins_DLP_validation$reads[vec_loc])
        sims_observed_CN <- CNbins_DLP_validation$state[vec_loc] * (gc_slope * CNbins_DLP_validation$gc[vec_loc] + gc_int)
        sims_noisy_CN <- rgamma(n = length(vec_loc), shape = sims_observed_CN / mean_sigma, scale = mean_sigma)
        sims_noisy_CN_pval <- sims_noisy_CN / sum(sims_noisy_CN)
        sims_reads <- rmultinom(n = 1, size = num_reads, prob = sims_noisy_CN_pval)
        CNbins_DLP_validation$reads_simulated[vec_loc] <- sims_reads
    }
    cat("\n")
    #--------------------------Plot distribution of total reads per cell
    filename <- paste(filename_prefix, "_total_reads.jpeg", sep = "")
    jpeg(file = filename, width = 1000, height = 1000)
    p <- ggplot(data.frame(total_reads = vec_total_reads), aes(x = total_reads)) +
        geom_histogram() +
        geom_histogram(color = "darkblue", fill = "lightblue") +
        geom_density(alpha = .2, fill = "#FF6666") +
        geom_vline(aes(xintercept = mean(total_reads)),
            color = "blue", linetype = "dashed", size = 1
        ) +
        xlab("Total reads per cell") +
        ylab("Count") +
        theme(panel.background = element_rect(fill = "white", colour = "grey50")) +
        theme(text = element_text(size = 20)) +
        scale_x_continuous(expand = c(0, 0)) +
        scale_y_continuous(expand = c(0, 0))
    print(p)
    dev.off()
    #-------------------------------Plot normalized DLP data for fitting
    filename <- paste(filename_prefix, "_normalized_reads_for_fitting.jpeg", sep = "")
    jpeg(file = filename, width = 1000, height = 1000)
    p <- ggplot(CNbins_fit, aes(x = gc, y = reads_normalized)) +
        geom_point(color = "cyan") +
        # geom_smooth(na.rm = TRUE, method = "lm", se = TRUE)+
        xlab("GC content") +
        ylab("Normalized readcounts") +
        theme(panel.background = element_rect(fill = "white", colour = "grey50")) +
        theme(text = element_text(size = 20)) +
        scale_x_continuous(expand = c(0, 0)) +
        scale_y_continuous(expand = c(0, 0))
    print(p)
    dev.off()
    #------------------------------Plot posterior distribution for sigma
    filename <- paste(filename_prefix, "_sigma_posterior.jpeg", sep = "")
    jpeg(file = filename, width = 1000, height = 1000)
    p <- ggplot(data.frame(vec_sigma_posterior = vec_sigma_posterior), aes(x = vec_sigma_posterior)) +
        geom_histogram() +
        geom_histogram(color = "coral", fill = "lightpink") +
        geom_vline(aes(xintercept = mean(vec_sigma_posterior)),
            color = "red", linetype = "dashed", size = 1
        ) +
        xlab("Sigma") +
        ylab("Posterior distribution") +
        theme(panel.background = element_rect(fill = "white", colour = "grey50")) +
        theme(text = element_text(size = 20)) +
        scale_x_continuous(expand = c(0, 0)) +
        scale_y_continuous(expand = c(0, 0))
    print(p)
    dev.off()
    #---------------------------Plot validation of readcount-GC per bins
    #   Find standardized residual errors for both DLP and simulated reads
    validation_groups <- sort(unique(CNbins_DLP_validation$multiplier), decreasing = TRUE)
    residual_DLP <- rep(0, length = length(validation_groups))
    residual_SIM <- rep(0, length = length(validation_groups))
    x_final <- 0.6
    y_final_DLP <- rep(0, length = length(validation_groups))
    y_final_SIM <- rep(0, length = length(validation_groups))
    for (i in 1:length(validation_groups)) {
        #   Perform linear model for DLP/simulated reads for the CN group
        vec_loc <- which(CNbins_DLP_validation$multiplier == validation_groups[i])
        x <- CNbins_DLP_validation$gc[vec_loc]
        y <- CNbins_DLP_validation$reads[vec_loc]
        lm.result.DLP <- lm(y ~ x)
        y <- CNbins_DLP_validation$reads_simulated[vec_loc]
        lm.result.SIM <- lm(y ~ x)
        #   Get the standardized residual errors for DLP/simulated reads
        residual_DLP[i] <- sigma(lm.result.DLP)
        residual_SIM[i] <- sigma(lm.result.SIM)
        #   Find the height of the final position
        y_final_DLP[i] <- (lm.result.DLP$coefficients[1] + x_final * lm.result.DLP$coefficients[2])
        y_final_SIM[i] <- (lm.result.SIM$coefficients[1] + x_final * lm.result.SIM$coefficients[2])
    }
    #   Plot validation of readcount-GC per bins
    CNbins_DLP_validation$CN <- paste(CNbins_DLP_validation$multiplier, "x ploidy", sep = "")
    CNbins_DLP_validation$CN <- factor(CNbins_DLP_validation$CN, levels = paste(sort(unique(CNbins_DLP_validation$multiplier), decreasing = TRUE), "x ploidy", sep = ""))
    filename <- paste(filename_prefix, "_readcount_model_validation.jpeg", sep = "")
    jpeg(file = filename, width = 2000, height = 1000)
    p1 <- ggplot(CNbins_DLP_validation) +
        geom_point(aes(x = gc, y = reads, col = CN), alpha = 0.01) +
        geom_smooth(aes(x = gc, y = reads, col = CN), na.rm = TRUE, method = "lm", se = TRUE) +
        geom_smooth(aes(x = gc, y = reads_simulated, col = CN), na.rm = TRUE, method = "lm", se = TRUE, linetype = "dashed") +
        labs(title = "DLP data", x = "GC content", y = "Readcounts") +
        theme(panel.background = element_rect(fill = "white", colour = "grey50")) +
        theme(text = element_text(size = 20)) +
        scale_x_continuous(expand = c(0, 0)) +
        scale_y_continuous(limits = c(0, 4000), expand = c(0, 0)) +
        theme(legend.position = "bottom")
    for (i in 1:length(validation_groups)) {
        y_final <- y_final_DLP[i]
        lab <- paste("Standardized residual error = ", round(residual_DLP[i]), sep = "")
        p1 <- p1 + annotate("text", x = x_final, y = y_final + 100, label = lab, size = 10, hjust = 1)
    }
    p2 <- ggplot(CNbins_DLP_validation) +
        geom_point(aes(x = gc, y = reads_simulated, col = CN), alpha = 0.01) +
        geom_smooth(aes(x = gc, y = reads, col = CN), na.rm = TRUE, method = "lm", se = TRUE) +
        geom_smooth(aes(x = gc, y = reads_simulated, col = CN), na.rm = TRUE, method = "lm", se = TRUE, linetype = "dashed") +
        labs(title = "Simulated data", x = "GC content", y = "Readcounts") +
        theme(panel.background = element_rect(fill = "white", colour = "grey50")) +
        theme(text = element_text(size = 20)) +
        scale_x_continuous(expand = c(0, 0)) +
        scale_y_continuous(limits = c(0, 4000), expand = c(0, 0))
    for (i in 1:length(validation_groups)) {
        y_final <- y_final_DLP[i]
        lab <- paste("Standardized residual error = ", round(residual_SIM[i]), sep = "")
        p2 <- p2 + annotate("text", x = x_final, y = y_final + 100, label = lab, size = 10, hjust = 1)
    }
    g_legend <- function(a.gplot) {
        tmp <- ggplot_gtable(ggplot_build(a.gplot))
        leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
        legend <- tmp$grobs[[leg]]
        return(legend)
    }
    common_legend <- g_legend(p1)
    p <- grid.arrange(arrangeGrob(p1 + theme(legend.position = "none"),
        p2 + theme(legend.position = "none"),
        nrow = 1
    ),
    common_legend,
    nrow = 2,
    heights = c(10, 1)
    )
    print(p)
    dev.off()
}
