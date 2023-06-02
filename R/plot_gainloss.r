#' @export
plot_gainloss <- function(copynumber_sims,
                          copynumber_DATA,
                          ploidy_normalization_sims = FALSE,
                          ploidy_normalization_DATA = TRUE,
                          type_sample_DATA = "individual",
                          title = NULL,
                          filename,
                          arm_level = FALSE,
                          pos_centromeres = c(),
                          height = 500,
                          width = 2000) {
    plotcol <- "state"
    fillna <- TRUE
    cutoff <- 2
    #---------------------------Get gain/loss consensus from simulations
    delta_sims <- gainloss_SIMS(
        copynumber_sims,
        ploidy_normalization = ploidy_normalization_sims,
        use_rbindlist = TRUE,
        get_coordinates = TRUE
    )
    f1_sims <- delta_sims$delta_gain
    f2_sims <- delta_sims$delta_loss
    copynumber_coordinates <- delta_sims$copynumber_coordinates
    #---------------------------------Get gain/loss consensus from PCAWG
    if (type_sample_DATA == "individual") {
        delta_PCAWG <- gainloss_DATA(
            copynumber_DATA,
            copynumber_coordinates,
            ploidy_normalization = ploidy_normalization_DATA,
            use_rbindlist = TRUE,
            arm_level = arm_level, pos_centromeres = pos_centromeres
        )
        f1_data <- delta_PCAWG$delta_gain
        f2_data <- delta_PCAWG$delta_loss
    } else if (type_sample_DATA == "average" & arm_level == TRUE) {
        list_chr <- unique(copynumber_coordinates$chr)
        tmp <- copynumber_coordinates
        tmp$delta_gain <- 0
        tmp$delta_loss <- 0
        CNbin_length <- tmp$end[1] - tmp$start[1] + 1

        for (i in 1:nrow(copynumber_DATA)) {
            ID <- copynumber_DATA$Arm[i]
            gain <- copynumber_DATA$Amp_freq_all[i]
            loss <- -copynumber_DATA$Del_freq_all[i]
            chr <- substr(ID, 1, nchar(ID) - 1)
            arm <- substr(ID, nchar(ID), nchar(ID))
            centromere <- pos_centromeres$Centromere_location[which(pos_centromeres$Chromosome == chr)]
            if (arm == "p") {
                vec_rows <- which((tmp$chr == chr) & (tmp$end <= centromere * CNbin_length))
            } else if (arm == "q") {
                vec_rows <- which((tmp$chr == chr) & (tmp$start > centromere * CNbin_length))
            }
            tmp$delta_gain[vec_rows] <- gain
            tmp$delta_loss[vec_rows] <- loss
        }
        f1_data <- tmp$delta_gain
        f2_data <- tmp$delta_loss
    }
    #-----------Find Spearman correlations between PCAWG and simulations
    stat_gain <- cor.test(f1_sims, f1_data, method = "spearman", exact = FALSE)
    rho_gain <- stat_gain$estimate[["rho"]]
    pval_gain <- stat_gain$p.value
    stat_loss <- cor.test(f2_sims, f2_data, method = "spearman", exact = FALSE)
    rho_loss <- stat_loss$estimate[["rho"]]
    pval_loss <- stat_loss$p.value
    #---------------Make genome-wide gain_sims/loss_sims plot comparison
    df_plot <- data.frame(
        x = 1:length(f1_sims),
        gain_sims = f1_sims,
        loss_sims = f2_sims,
        gain_data = f1_data,
        loss_data = f2_data
    )
    if (!is.null(title)) {
        p_title <- ggplot() +
            annotate("text", x = 0, y = 0, size = 12, colour = "black", label = title) +
            theme_void()
    }
    #---Plot gain_sims/loss_sims map for simulations
    p_sims <- ggplot(data = df_plot) +
        geom_bar(aes(x = x, y = gain_sims), stat = "identity", colour = "#E34A33", fill = "#E34A33") +
        geom_bar(aes(x = x, y = loss_sims), stat = "identity", colour = "#3182BD", fill = "#3182BD") +
        geom_hline(yintercept = 0, color = "antiquewhite3") +
        theme(panel.background = element_rect(fill = "white", colour = "grey50")) +
        theme(text = element_text(size = 30)) +
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
        theme(text = element_text(size = 30)) +
        scale_x_continuous(limits = c(0, nrow(df_plot)), expand = c(0, 0)) +
        scale_y_continuous(limits = c(-1, 1), expand = c(0, 0)) +
        ylab("Data") +
        xlab("")
    #---Plot Spearman correlation scores
    p_spearman_gain <- ggplot() +
        annotate("text", x = 0, y = 0, size = 12, colour = "#E34A33", label = paste("GAIN: rho = ", round(rho_gain, 2), sep = "")) +
        theme_void()
    p_spearman_loss <- ggplot() +
        annotate("text", x = 0, y = 0, size = 12, colour = "#3182BD", label = paste("LOSS: rho = ", round(rho_loss, 2), sep = "")) +
        theme_void()
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
        theme(axis.text.x = element_text(size = 20)) +
        theme(plot.margin = unit(c(1, 0.5, 0, 0.5), "cm"))
    p_data <- p_data + scale_x_continuous(
        breaks = chr_ticks,
        labels = chr_id,
        limits = c(0, nrow(df_plot)), expand = c(0, 0)
    ) +
        theme(axis.text.x = element_text(size = 20)) +
        theme(plot.margin = unit(c(1, 0.5, 0, 0.5), "cm"))
    #---Print the mini plots
    jpeg(filename, width = width, height = height)
    if (!is.null(title)) {
        p <- grid.arrange(p_title, p_sims, arrangeGrob(p_spearman_gain, p_spearman_loss, nrow = 1), p_data, heights = c(1, 5, 1, 5), ncol = 1)
    } else {
        p <- grid.arrange(p_sims, arrangeGrob(p_spearman_gain, p_spearman_loss, nrow = 1), p_data, heights = c(5, 1, 5), ncol = 1)
    }
    print(p)
    dev.off()
}

#' @export
gainloss_DATA <- function(copynumber_DATA,
                          copynumber_coordinates,
                          ploidy_normalization = FALSE,
                          use_rbindlist = FALSE,
                          arm_level = FALSE,
                          state_mode = NULL,
                          round = TRUE,
                          pos_centromeres) {
    library(signals)
    plotcol <- "state"
    fillna <- TRUE
    cutoff <- 2
    #---------------------------------Get gain/loss consensus from PCAWG
    samplelist_data <- unique(copynumber_DATA$donor_unique_id)
    CNbins_list_data <- vector("list", length = length(samplelist_data))
    NA_list_data <- vector("list", length = length(samplelist_data))
    CNbin_length <- copynumber_coordinates$end[1] - copynumber_coordinates$start[1] + 1
    for (iteration in 1:length(samplelist_data)) {
        sample_id <- samplelist_data[iteration]
        #   Get CN data for the sample
        cn_profiles_data <- copynumber_DATA[which(copynumber_DATA$donor_unique_id == sample_id), ]
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
    if (use_rbindlist == TRUE) {
        CNbins_data <- rbindlist(CNbins_list_data, use.names = FALSE, fill = FALSE, idcol = NULL)
        class(CNbins_data) <- "data.frame"
    }
    #---Transform CN table into copynumber_sims object
    if (use_rbindlist == TRUE) {
        copynumber_data <- createCNmatrix(CNbins_data,
            field = plotcol, wholegenome = FALSE,
            fillnaplot = fillna, centromere = FALSE
        )
    } else {
        copynumber_data <- CNbins_list_data[[1]]
        copynumber_data <- copynumber_data[, 1:3]
        copynumber_data$width <- copynumber_data$end - copynumber_data$start + 1
        for (iteration in 1:length(CNbins_list_data)) {
            CNbins_iteration <- CNbins_list_data[[iteration]]
            cell_id <- CNbins_iteration$cell_id[1]
            copynumber_data[[cell_id]] <- CNbins_iteration$copy
        }
    }
    #---Condense CN to arm level
    if (arm_level == TRUE) {
        list_chr <- unique(copynumber_data$chr)
        for (i in 1:length(list_chr)) {
            chr <- list_chr[i]
            centromere <- pos_centromeres$Centromere_location[which(pos_centromeres$Chromosome == chr)]
            length <- pos_centromeres$Bin_count[which(pos_centromeres$Chromosome == chr)]
            for (arm in 1:2) {
                if (arm == 1) {
                    vec_rows <- which((copynumber_data$chr == chr) & (copynumber_data$end <= centromere * CNbin_length))
                } else {
                    vec_rows <- which((copynumber_data$chr == chr) & (copynumber_data$start > centromere * CNbin_length))
                }
                for (col in 5:ncol(copynumber_data)) {
                    copynumber_data[vec_rows, col] <- calc_state_mode(copynumber_data[vec_rows, col])
                }
            }
        }
    }
    #---Normalize ploidy of each sample to 2
    if (ploidy_normalization == TRUE) {
        copynumber_data <- normalize_cell_ploidy(copynumber_data, state_mode, round)
    }
    #---Replace previously NA in each sample with ploidy 2
    if (arm_level != TRUE) {
        for (sample in 1:length(samplelist_data)) {
            sample_id <- samplelist_data[sample]
            NA_iteration <- NA_list_data[[sample]]
            vec_loc <- which(NA_iteration == 1)
            if (length(vec_loc) > 0) {
                tmp <- copynumber_data[sample_id]
                tmp[vec_loc, ] <- 2
                copynumber_data[sample_id] <- tmp
            }
        }
    }
    #---Find gain_sims/loss_sims maps
    n_samples <- ncol(copynumber_data) - 4
    f1_data <- rowSums(copynumber_data[, 5:ncol(copynumber_data)] > cutoff) / n_samples
    f2_data <- -rowSums(copynumber_data[, 5:ncol(copynumber_data)] < cutoff) / n_samples
    attr(f1_data, "names") <- NULL
    attr(f2_data, "names") <- NULL
    #------------------------------Output gain/loss consensus from PCAWG
    output <- list()
    output$delta_gain <- f1_data
    output$delta_loss <- f2_data
    output$copynumber_data <- copynumber_data
    return(output)
}

#' @export
gainloss_SIMS <- function(copynumber_sims,
                          ploidy_normalization = FALSE,
                          use_rbindlist = FALSE,
                          get_coordinates = FALSE,
                          state_mode = NULL,
                          get_CN = FALSE,
                          round = TRUE,
                          get_WGD_status = FALSE) {
    plotcol <- "state"
    fillna <- TRUE
    cutoff <- 2
    #---------------------------Get gain/loss consensus from simulations
    CNbins_list_sims <- vector("list", length = length(copynumber_sims))
    if (get_WGD_status) {
        wgd_status_sims <- rep("", length(copynumber_sims))
    }
    for (iteration in 1:length(copynumber_sims)) {
        simulation <- copynumber_sims[[iteration]]
        all_sample_genotype <- simulation$sample$all_sample_genotype
        sample_genotype_unique <- simulation$sample$sample_genotype_unique
        sample_genotype_unique_profile <- simulation$sample$sample_genotype_unique_profile
        #   Get the genotype with the highest clonal percentage in sample
        freqs <- rep(0, length(sample_genotype_unique))
        for (i in 1:length(sample_genotype_unique)) {
            freqs[i] <- length(which(all_sample_genotype == sample_genotype_unique[i]))
        }
        max_freq <- max(freqs)
        vec_loc <- which(freqs == max_freq)
        if (length(vec_loc) > 1) {
            max_loc <- sample(vec_loc, 1)
        } else {
            max_loc <- vec_loc
        }
        max_genotype <- sample_genotype_unique[max_loc]
        #   Add record for this cell to list
        CNbins_iteration <- sample_genotype_unique_profile[[max_loc]]
        CNbins_iteration$cell_id <- paste("SIMULATION", iteration, "-Library-1-1", sep = "")
        CNbins_list_sims[[iteration]] <- CNbins_iteration
        #   Find WGD status if inquired
        if (get_WGD_status) {
            evolution_genotype_changes <- simulation$clonal_evolution$evolution_genotype_changes
            evolution_origin <- simulation$clonal_evolution$evolution_origin
            WGD_status <- "no_wgd"
            genotype <- max_genotype
            while (genotype > 0) {
                genotype_changes <- evolution_genotype_changes[[genotype]]
                for (genotype_change in genotype_changes) {
                    if (genotype_change[1] == "whole-genome-duplication") {
                        WGD_status <- "wgd"
                    }
                }
                genotype <- evolution_origin[genotype]
            }
            wgd_status_sims[iteration] <- WGD_status
        }
    }
    if (use_rbindlist == TRUE) {
        CNbins_sims <- rbindlist(CNbins_list_sims, use.names = FALSE, fill = FALSE, idcol = NULL)
        class(CNbins_sims) <- "data.frame"
    }
    #---Transform CN table into copynumber_sims object
    if (use_rbindlist == TRUE) {
        copynumber_sims <- createCNmatrix(CNbins_sims,
            field = plotcol, wholegenome = FALSE,
            fillnaplot = fillna, centromere = FALSE
        )
    } else {
        copynumber_sims <- CNbins_list_sims[[1]]
        copynumber_sims <- copynumber_sims[, 1:3]
        copynumber_sims$width <- copynumber_sims$end - copynumber_sims$start + 1
        for (iteration in 1:length(CNbins_list_sims)) {
            CNbins_iteration <- CNbins_list_sims[[iteration]]
            cell_id <- CNbins_iteration$cell_id[1]
            copynumber_sims[[cell_id]] <- CNbins_iteration$copy
        }
    }
    #---Normalize ploidy of each sample to 2
    if (ploidy_normalization == TRUE) {
        copynumber_sims <- normalize_cell_ploidy(copynumber_sims, state_mode, round)
    }
    #---Get genome coordinates
    if (get_coordinates == TRUE) {
        copynumber_coordinates <- copynumber_sims[, 1:4]
    }
    #---Find gain_sims/loss_sims maps
    n_samples <- ncol(copynumber_sims) - 4
    f1_sims <- rowSums(copynumber_sims[, 5:ncol(copynumber_sims)] > cutoff) / n_samples
    f2_sims <- -rowSums(copynumber_sims[, 5:ncol(copynumber_sims)] < cutoff) / n_samples
    attr(f1_sims, "names") <- NULL
    attr(f2_sims, "names") <- NULL
    #------------------------------Output gain/loss consensus from PCAWG
    output <- list()
    output$delta_gain <- f1_sims
    output$delta_loss <- f2_sims
    if (get_coordinates == TRUE) {
        output$copynumber_coordinates <- copynumber_coordinates
    }
    if (get_CN == TRUE) {
        output$copynumber_sims <- copynumber_sims
    }
    if (get_WGD_status == TRUE) {
        output$wgd_status_sims <- wgd_status_sims
    }
    return(output)
}

calc_state_mode <- function(states) {
    state_levels <- unique(states)

    state_levels <- state_levels[!(state_levels %in% c(0))]

    state_mode <- state_levels[
        which.max(tabulate(match(states, state_levels)))
    ]
    if (!is.finite(state_mode)) {
        state_mode <- 2
    }
    return(state_mode)
}

normalize_cell_ploidy <- function(copynumber, state_mode, round = TRUE) {
    cell_ids <- colnames(copynumber)
    cell_ids <- cell_ids[!(cell_ids %in% c("chr", "start", "end", "width"))]
    for (cell_id in cell_ids) {
        if (is.null(state_mode)) {
            state_mode_cell <- calc_state_mode(copynumber[[cell_id]])
        } else {
            state_mode_cell <- state_mode
        }
        copynumber[[cell_id]] <- copynumber[[cell_id]] / (state_mode_cell / 2)
        if (round == TRUE) {
            copynumber[[cell_id]] <- as.integer(ceiling(copynumber[[cell_id]]))
        }
        copynumber[[cell_id]][copynumber[[cell_id]] > 11] <- 11
    }
    return(copynumber)
}

densityPlot_MODIFIED <- function(object,
                                 obs,
                                 training,
                                 add = TRUE,
                                 main = "Posterior density",
                                 color_prior,
                                 chosen_para = NULL,
                                 color_posterior,
                                 protocol,
                                 color_vline,
                                 log = "",
                                 xlim = NULL,
                                 ylim = NULL,
                                 xlab = NULL,
                                 ylab = NULL,
                                 paral = FALSE,
                                 fontsize = 50,
                                 ncores = if (paral) max(detectCores() - 1, 1) else 1, ...) {
    df_plot <- densityPlot_df(
        object,
        obs,
        training,
        add,
        main,
        color_prior,
        chosen_para,
        color_posterior,
        protocol,
        color_vline,
        log,
        xlim,
        ylim,
        xlab,
        ylab,
        paral,
        ncores
    )

    p_plot <- ggplot(df_plot) +
        geom_area(aes(x = x, y = y_prior), color = color_prior, fill = color_prior, alpha = 0.3) +
        geom_area(aes(x = x, y = y_posterior), color = color_posterior, fill = color_posterior, alpha = 0.3) +
        # geom_density(aes(x = dist_raw, kernel = "gaussian", weight = weight_prior), color = color_prior, fill = color_prior, alpha = 0.3) +
        # geom_density(aes(x = dist_raw, kernel = "gaussian", weight = weight_posterior), color = color_posterior, fill = color_posterior, alpha = 0.3) +
        xlab("") +
        ylab("") +
        ggtitle(main) +
        theme(panel.background = element_rect(fill = "white", colour = "grey50")) +
        theme(text = element_text(size = fontsize)) +
        scale_x_continuous(expand = c(0, 0)) +
        scale_y_continuous(expand = c(0, 0))

    print(chosen_para)
    if (is.null(chosen_para) == FALSE) {
        p_plot <- p_plot +
            geom_vline(aes(xintercept = chosen_para), color = color_vline, size = 1, linetype = "dotted")
    }

    return(p_plot)
}

densityPlot_df <- function(object,
                           obs,
                           training,
                           add = TRUE,
                           main = "Posterior density",
                           color_prior,
                           chosen_para = NULL,
                           color_posterior,
                           protocol = "",
                           color_vline,
                           log = "",
                           xlim = NULL,
                           ylim = NULL,
                           xlab = NULL,
                           ylab = NULL,
                           paral = FALSE,
                           ncores = if (paral) max(detectCores() - 1, 1) else 1, ...) {
    findweights <- getFromNamespace("findweights", "abcrf")
    ### Checking arguments
    if (!inherits(object, "regAbcrf")) {
        stop("object not of class regAbcrf")
    }

    if (!inherits(training, "data.frame")) {
        stop("training needs to be a data.frame object")
    }

    if (!inherits(obs, "data.frame")) {
        stop("obs needs to be a data.frame object")
    }
    if (nrow(obs) == 0L || is.null(nrow(obs))) {
        stop("no data in obs")
    }
    if (nrow(training) == 0L || is.null(nrow(training))) {
        stop("no simulation in the training reference table (response, sumstat)")
    }

    if ((!is.logical(add)) || (length(add) != 1L)) {
        stop("add should be TRUE or FALSE")
    }
    if ((!is.logical(paral)) || (length(paral) != 1L)) {
        stop("paral should be TRUE or FALSE")
    }
    if (is.na(ncores)) {
        warning("Unable to automatically detect the number of CPU cores, \n1 CPU core will be used or please specify ncores.")
        ncores <- 1
    }

    if (!is.character(log)) {
        stop("log needs to be a character string")
    }
    x <- obs
    if (!is.null(x)) {
        if (is.vector(x)) {
            x <- matrix(x, ncol = 1)
        }
        if (nrow(x) == 0) {
            stop("obs has 0 rows")
        }
        if (any(is.na(x))) {
            stop("missing values in obs")
        }
    }

    # resp and sumsta recover

    mf <- match.call(expand.dots = FALSE)
    mf <- mf[1]
    mf$formula <- object$formula


    mf$data <- training


    training <- mf$data

    mf[[1L]] <- as.name("model.frame")
    mf <- eval(mf, parent.frame())
    mt <- attr(mf, "terms")
    resp <- model.response(mf)

    obj <- object$model.rf
    inbag <- matrix(unlist(obj$inbag.counts, use.names = FALSE), ncol = obj$num.trees, byrow = FALSE)

    obj[["origNodes"]] <- predict(obj, training, predict.all = TRUE, num.threads = ncores)$predictions
    obj[["origObs"]] <- model.response(mf)

    #####################

    origObs <- obj$origObs
    origNodes <- obj$origNodes

    nodes <- predict(obj, x, predict.all = TRUE, num.threads = ncores)$predictions
    if (is.null(dim(nodes))) nodes <- matrix(nodes, nrow = 1)
    ntree <- obj$num.trees
    nobs <- object$model.rf$num.samples
    nnew <- nrow(x)

    weights <- findweights(origNodes, nodes, inbag, as.integer(nobs), as.integer(nnew), as.integer(ntree)) # cpp function call
    weights.std <- weights / ntree

    priorDensity <- density(resp)

    if (add) {
        rangex <- range(priorDensity$x)
        rangey <- range(priorDensity$y)

        for (i in 1:nnew) {
            postDensity <- density(resp, weights = weights.std[, i], ...)
            rangex <- range(rangex, postDensity$x)
            rangey <- range(rangey, postDensity$y)
        }

        # plot(priorDensity$x, priorDensity$y,
        #     type = "l", main = main, log = log,
        #     xlim = if (is.null(xlim)) rangex else xlim,
        #     ylim = if (is.null(ylim)) rangey else ylim,
        #     xlab = xlab, ylab = ylab, col = "grey"
        # )
        # for (i in 1:nnew) {
        #     postDensity <- density(resp, weights = weights.std[, i], ...)
        #     points(postDensity$x, postDensity$y, type = "l")
        # }
    } else {
        for (i in 1:nnew) {
            postDensity <- density(resp, weights = weights.std[, i], ...)

            # plot(postDensity$x, postDensity$y,
            #     type = "l", main = main, log = log,
            #     xlim = if (is.null(xlim)) range(postDensity$x, priorDensity$x) else xlim,
            #     ylim = if (is.null(ylim)) range(postDensity$y, priorDensity$y) else ylim,
            #     xlab = xlab, ylab = ylab
            # )
            # points(priorDensity$x, priorDensity$y, type = "l", col = "grey")
            # if (nnew > 1 && i < nnew) readline("Press <ENTER> to Continue")
        }
    }



    if (protocol == "TSG") {
        resp <- 1 / resp
    }
    dist_prior <- density(resp, weights = rep(1 / length(resp), length(resp)))
    dist_posterior <- density(resp, weights = weights.std[, i])

    df_plot_prior <- data.frame(x = dist_prior$x, y = dist_prior$y)
    df_plot_posterior <- data.frame(x = dist_posterior$x, y = dist_posterior$y)

    df_plot <- data.frame(x = dist_prior$x, y_prior = dist_prior$y, y_posterior = dist_posterior$y)

    # df_plot <- data.frame(dist_raw = resp)
    # df_plot$weight_prior <- 1 / nrow(df_plot)
    # df_plot$weight_posterior <- weights.std[, i]

    return(df_plot)
}
