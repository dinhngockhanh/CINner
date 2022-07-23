analysis_dlp_cn <- function(CNbins, vec_gc, filename_prefix) {
    CNbins_DLP <- CNbins
    #-----------------------------------Set parameters for fitting sigma
    gc_slope <- 1
    gc_int <- 0
    fit_ploidy <- 2
    min_reads <- 100

    fit_ABC_count <- 100000
    fit_ABC_sigma_min <- 0.001
    fit_ABC_sigma_max <- 0.1
    fit_ABC_tol <- 0.01
    fit_bin_count <- 1000
    #--------------------------Find distribution of total reads per cell
    list_cell_id <- unique(CNbins_DLP$cell_id)
    vec_total_reads <- rep(0, length = length(list_cell_id))
    for (i in 1:length(list_cell_id)) {
        vec_total_reads[i] <- sum(CNbins_DLP$reads[CNbins_DLP$cell_id == list_cell_id[i]])
    }
    #-------------------------Supplement the CN data with GC/mappability
    cat("Add vec_gc/mappability data to CNbins...\n")
    pb <- txtProgressBar(
        min = 0, max = nrow(vec_gc),
        style = 3, width = 50, char = "="
    )
    CNbins_DLP$gc <- -1
    CNbins_DLP$map <- -1
    for (row in 1:nrow(vec_gc)) {
        setTxtProgressBar(pb, row)
        vec_loc <- which(CNbins_DLP$chr == vec_gc$chr[row] & CNbins_DLP$start_point == vec_gc$start[row] & CNbins_DLP$end_point == vec_gc$end[row])
        CNbins_DLP$gc[vec_loc] <- vec_gc$gc[row]
        CNbins_DLP$map[vec_loc] <- vec_gc$map[row]
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
    #   Normalize reads for vec_gc and total reads per cell
    list_cell_id <- unique(CNbins_fit$cell_id)
    for (i in 1:length(list_cell_id)) {
        vec_loc <- which(CNbins_fit$cell_id == list_cell_id[i])
        #   Normalize reads for vec_gc
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
    #-----------------------Plot validation of readcount-vec_gc per bins
    CNbins_DLP_validation$CN <- paste(CNbins_DLP_validation$multiplier, "x ploidy", sep = "")
    CNbins_DLP_validation$CN <- factor(CNbins_DLP_validation$CN, levels = paste(sort(unique(CNbins_DLP_validation$multiplier), decreasing = TRUE), "x ploidy", sep = ""))
    filename <- paste(filename_prefix, "_readcount_model_validation.jpeg", sep = "")
    jpeg(file = filename, width = 1000, height = 1000)
    p <- ggplot(CNbins_DLP_validation) +
        geom_point(aes(x = gc, y = reads, col = CN), alpha = 0.01) +
        # geom_point(aes(x = gc, y = reads, col = CN), shape = ".", alpha = 0.01) +
        geom_smooth(aes(x = gc, y = reads, col = CN), na.rm = TRUE, method = "lm", se = TRUE) +
        geom_smooth(aes(x = gc, y = reads_simulated, col = CN), na.rm = TRUE, method = "lm", se = TRUE, linetype = "dashed") +
        xlab("GC content") +
        ylab("Readcounts") +
        theme(panel.background = element_rect(fill = "white", colour = "grey50")) +
        theme(text = element_text(size = 20)) +
        scale_x_continuous(expand = c(0, 0)) +
        scale_y_continuous(limits = c(0, 4000), expand = c(0, 0))
    print(p)
    dev.off()
}
