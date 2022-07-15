#' @export
plot_gainloss <- function(copynumber_sims,
                          copynumber_PCAWG) {
    plotcol <- "state"
    fillna <- TRUE
    cutoff <- 2
    #---------------------------Get gain/loss consensus from simulations
    CNbins_list_sims <- vector("list", length = length(copynumber_sims))
    for (iteration in 1:length(copynumber_sims)) {
        simulation <- copynumber_sims[[iteration]]
        all_sample_genotype <- simulation$sample$all_sample_genotype
        sample_cell_ID <- simulation$sample$sample_cell_ID
        cn_profiles_long <- simulation$sample$cn_profiles_long
        #   Get the genotype with the highest clonal percentage in sample
        tmp <- as.data.frame(table(all_sample_genotype))
        max_freq <- max(tmp$Freq)
        vec_loc <- which(tmp$Freq == max_freq)
        if (length(vec_loc) > 1) {
            loc <- sample(vec_loc, 1)
        } else {
            loc <- vec_loc
        }
        max_genotype <- tmp$all_sample_genotype[loc]
        #   Get a cell from this genotype
        max_cell <- which(all_sample_genotype == max_genotype)[1]
        max_cell_id <- sample_cell_ID[max_cell]
        #   Add record for this cell to list
        CNbins_iteration <- cn_profiles_long[which(cn_profiles_long$cell_id == max_cell_id), ]
        CNbins_iteration$cell_id <- paste("SIMULATION", iteration, "-Library-1-1", sep = "")
        CNbins_list_sims[[iteration]] <- CNbins_iteration
    }
    CNbins_sims <- rbindlist(CNbins_list_sims, use.names = FALSE, fill = FALSE, idcol = NULL)
    class(CNbins_sims) <- "data.frame"
    #---Transform CN table into copynumber_sims object
    copynumber_sims <- createCNmatrix(CNbins_sims,
        field = plotcol, wholegenome = TRUE,
        fillnaplot = fillna, centromere = FALSE
    )
    #---Normalize ploidy of each sample to 2
    copynumber_sims <- normalize_cell_ploidy(copynumber_sims)
    #---Get genome coordinates
    copynumber_coordinates <- copynumber_sims[, 1:4]
    #---Find gain_sims/loss_sims maps
    n_samples <- ncol(copynumber_sims) - 4
    f1_sims <- rowSums(copynumber_sims[, 5:ncol(copynumber_sims)] > cutoff) / n_samples
    f2_sims <- -rowSums(copynumber_sims[, 5:ncol(copynumber_sims)] < cutoff) / n_samples
    attr(f1_sims, "names") <- NULL
    attr(f2_sims, "names") <- NULL
    #---------------------------------Get gain/loss consensus from PCAWG
    samplelist_data <- unique(copynumber_PCAWG$donor_unique_id)
    CNbins_list_data <- vector("list", length = length(samplelist_data))
    NA_list_data <- vector("list", length = length(samplelist_data))
    CNbin_length <- copynumber_coordinates$end[1] - copynumber_coordinates$start[1] + 1
    for (iteration in 1:length(samplelist_data)) {
        sample_id <- samplelist_data[iteration]
        #   Get CN data for the sample
        cn_profiles_data <- copynumber_PCAWG[which(copynumber_PCAWG$donor_unique_id == sample_id), ]
        #   Get uniform-bin CN data from the data
        CNbins_iteration <- data.frame(chr = copynumber_coordinates$chr, start = copynumber_coordinates$start, end = copynumber_coordinates$end)
        CNbins_iteration$copy <- 0
        CNbins_iteration$state <- 0
        CNbins_iteration$Min <- 0
        CNbins_iteration$Maj <- 0
        CNbins_iteration$cell_id <- sample_id
        NA_iteration <- rep(1, length = nrow(CNbins_iteration))
        for (row in 1:nrow(cn_profiles_data)) {
            chr <- cn_profiles_data$chromosome[row]
            start <- cn_profiles_data$start[row]
            end <- cn_profiles_data$end[row]
            total_cn <- cn_profiles_data$total_cn[row]
            major_cn <- cn_profiles_data$major_cn[row]
            minor_cn <- cn_profiles_data$minor_cn[row]
            bin_start <- floor(start / CNbin_length) * CNbin_length + 1
            bin_end <- ceiling(start / CNbin_length) * CNbin_length
            vec_loc <- which(CNbins_iteration$chr == chr & CNbins_iteration$start >= bin_start & CNbins_iteration$end <= end)
            CNbins_iteration$copy[vec_loc] <- total_cn
            CNbins_iteration$state[vec_loc] <- total_cn
            CNbins_iteration$Min[vec_loc] <- minor_cn
            CNbins_iteration$Maj[vec_loc] <- major_cn
            NA_iteration[vec_loc] <- 0
        }
        #   Add record for this sample to list
        CNbins_list_data[[iteration]] <- CNbins_iteration
        NA_list_data[[iteration]] <- NA_iteration
    }
    CNbins_data <- rbindlist(CNbins_list_data, use.names = FALSE, fill = FALSE, idcol = NULL)
    class(CNbins_data) <- "data.frame"
    #---Transform CN table into copynumber_sims object
    copynumber_data <- createCNmatrix(CNbins_data,
        field = plotcol, wholegenome = TRUE,
        fillnaplot = fillna, centromere = FALSE
    )
    #---Normalize ploidy of each sample to 2
    copynumber_data <- normalize_cell_ploidy(copynumber_data)
    #---Find gain_sims/loss_sims maps
    n_samples <- ncol(copynumber_data) - 4
    f1_data <- rowSums(copynumber_data[, 5:ncol(copynumber_data)] > cutoff) / n_samples
    f2_data <- -rowSums(copynumber_data[, 5:ncol(copynumber_data)] < cutoff) / n_samples
    attr(f1_data, "names") <- NULL
    attr(f2_data, "names") <- NULL
    # #----------------------------------------------Rework the data frame
    # #   Conform data frame into format required by signals
    # CNbins_sims <- copynumber_sims_PCAWG[c("chromosome", "start", "end", "total_cn", "minor_cn", "major_cn", "donor_unique_id")]
    #
    # # CNbins_sims <- CNbins_sims[which(CNbins_sims$donor_unique_id == CNbins_sims$donor_unique_id[1]), ]
    #
    # CNbins_sims$state <- CNbins_sims$total_cn
    # CNbins_sims <- CNbins_sims[, c(1, 2, 3, 4, 8, 5, 6, 7)]
    # colnames(CNbins_sims) <- c("chr", "start", "end", "copy", "state", "Min", "Maj", "cell_id")
    # #
    #



    # tmp <- unique(CNbins_sims$cell_id)
    # for (i in 1:length(tmp)) {
    #     print(mean(CNbins_sims$copy[which(CNbins_sims$cell_id == tmp[i])]))
    #     print(var(CNbins_sims$copy[which(CNbins_sims$cell_id == tmp[i])]))
    # }


    #------------------------------------Make genome-wide gain_sims/loss_sims plot
    #-------------------------Make genome-wide gain_sims/loss_sims plot comparison
    df_plot <- data.frame(
        x = 1:length(f1_sims),
        gain_sims = f1_sims,
        loss_sims = f2_sims,
        gain_data = f1_data,
        loss_data = f2_data
    )
    #---Plot gain_sims/loss_sims map for simulations
    p_sims <- ggplot(data = df_plot) +
        geom_bar(aes(x = x, y = gain_sims), stat = "identity", colour = "#E34A33", fill = "#E34A33") +
        geom_bar(aes(x = x, y = loss_sims), stat = "identity", colour = "#3182BD", fill = "#3182BD") +
        geom_hline(yintercept = 0, color = "antiquewhite3") +
        theme(panel.background = element_rect(fill = "white", colour = "grey50")) +
        theme(text = element_text(size = 20)) +
        scale_x_continuous(limits = c(0, nrow(df_plot)), expand = c(0, 0)) +
        scale_y_continuous(limits = c(-1, 1), expand = c(0, 0)) +
        ylab("Simulations") +
        xlab("")
    #---Plot gain_sims/loss_sims map for data
    p_data <- ggplot(data = df_plot) +
        geom_bar(aes(x = x, y = gain_data), stat = "identity", colour = "#E34A33", fill = "#E34A33") +
        geom_bar(aes(x = x, y = loss_data), stat = "identity", colour = "#3182BD", fill = "#3182BD") +
        geom_hline(yintercept = 0, color = "antiquewhite3") +
        theme(panel.background = element_rect(fill = "white", colour = "grey50")) +
        theme(text = element_text(size = 20)) +
        scale_x_continuous(limits = c(0, nrow(df_plot)), expand = c(0, 0)) +
        scale_y_continuous(limits = c(-1, 1), expand = c(0, 0)) +
        ylab("PCAWG") +
        xlab("")
    #---Plot chromosomes
    #   Find chromosome information
    chr_id <- unique(copynumber_coordinates$chr)
    chr_bin_count <- rep(0, length = length(chr_id))
    chr_ticks <- rep(0, length = length(chr_id))
    chr_ends <- rep(0, length = length(chr_id))
    for (chr in 1:length(chr_id)) {
        chr_bin_count[chr] <- length(which(copynumber_coordinates$chr == chr_id[chr]))
        if (chr == 1) {
            chr_ticks[chr] <- chr_bin_count[chr] / 2
        } else {
            chr_ticks[chr] <- sum(chr_bin_count[1:(chr - 1)]) + chr_bin_count[chr] / 2
        }
        chr_ends[chr] <- sum(chr_bin_count[1:chr])
    }
    #   Plot chromosome separations
    for (chr in 1:length(chr_id)) {
        p_sims <- p_sims + geom_vline(xintercept = chr_ends[chr], color = "antiquewhite3")
        p_data <- p_data + geom_vline(xintercept = chr_ends[chr], color = "antiquewhite3")
    }
    #   Plot chromosome names
    p_sims <- p_sims + scale_x_continuous(
        breaks = chr_ticks,
        labels = chr_id,
        limits = c(0, nrow(df_plot)), expand = c(0, 0)
    ) +
        theme(axis.text.x = element_text(size = 10))
    p_data <- p_data + scale_x_continuous(
        breaks = chr_ticks,
        labels = chr_id,
        limits = c(0, nrow(df_plot)), expand = c(0, 0)
    ) +
        theme(axis.text.x = element_text(size = 10))

    grid.arrange(
        p_sims,
        p_data,
        ncol = 1
    )

    p <- arrangeGrob(p_sims, p_data, ncol = 1)

    print(p)
}

calc_state_mode <- function(states) {
    state_levels <- unique(states)
    state_mode <- state_levels[
        which.max(tabulate(match(states, state_levels)))
    ]
    if (!is.finite(state_mode)) {
        state_mode <- 2
    }
    return(state_mode)
}

normalize_cell_ploidy <- function(copynumber_sims) {
    cell_ids <- colnames(copynumber_sims)
    cell_ids <- cell_ids[!(cell_ids %in% c("chr", "start", "end", "width"))]

    for (cell_id in cell_ids) {
        state_mode <- calc_state_mode(copynumber_sims[[cell_id]])
        copynumber_sims[[cell_id]] <- as.integer(ceiling(
            copynumber_sims[[cell_id]] / (state_mode / 2)
        ))
        copynumber_sims[[cell_id]][copynumber_sims[[cell_id]] > 11] <- 11
    }
    return(copynumber_sims)
}
