#' @export
plot_cn_per_clone <- function(model = "",
                              n_simulations = 0,
                              CN_data = "TRUTH",
                              width = 1000,
                              height = 1000) {
    for (i in 1:n_simulations) {
        #------------------------------------------Input simulation file
        filename <- paste(model, "_simulation_", i, ".rda", sep = "")
        load(filename)
        #--------------------------------------------Extract CN profiles
        if (CN_data == "TRUTH") {
            #---Extract CN profiles from GROUND TRUTH
            sample_genotype_profiles <- simulation$sample$cn_profiles_long
            #   Remove internal nodes from true CN profiles
            vec_delete <- c(
                grep("Internal-node-", sample_genotype_profiles$cell_id, value = FALSE),
                grep("Initial-clone-", sample_genotype_profiles$cell_id, value = FALSE)
            )
            sample_genotype_profiles <- sample_genotype_profiles[-vec_delete, ]
        }
        #--------------------------------Plot CN profiles for each clone
        filename <- paste(model, "_sim", i, "_CN_total_per_clone.jpeg", sep = "")
        jpeg(file = filename, width = width, height = height)
        p <- plotCNprofile_MODIFIED(sample_genotype_profiles, legend.position = "right")
        print(p)
        dev.off()
    }
}

#' @export
plotCNprofile_MODIFIED <- function(CNbins,
                                   cellid = NULL,
                                   chrfilt = NULL,
                                   pointsize = 1,
                                   alphaval = 0.6,
                                   maxCN = 10,
                                   cellidx = 1,
                                   statecol = "state",
                                   returnlist = FALSE,
                                   raster = FALSE,
                                   y_axis_trans = "identity",
                                   xaxis_order = "genome_position",
                                   legend.position = "bottom",
                                   annotateregions = NULL,
                                   SV = NULL,
                                   svalpha = 0.5,
                                   svwidth = 1.0,
                                   adj = 0.03,
                                   genes = NULL,
                                   tickwidth = 50,
                                   chrstart = NULL,
                                   chrend = NULL,
                                   shape = 16,
                                   ...) {
    if (!xaxis_order %in% c("bin", "genome_position")) {
        stop("xaxis_order must be either 'bin' or 'genome_position'")
    }

    if (is.null(cellid)) {
        cellid <- unique(CNbins$cell_id)[min(cellidx, length(unique(CNbins$cell_id)))]
    }

    if (y_axis_trans == "squashy") {
        ybreaks <- c(0, 2, 5, 10, maxCN)
    } else {
        ybreaks <- seq(0, maxCN, 2)
    }

    if (length(chrfilt) == 1) {
        xlab <- paste0("Chromosome ", chrfilt)
    } else {
        xlab <- "Chromosome"
    }

    statecolpal <- scCNstate_cols()

    message(paste0("Making CN profile and BAF plot for cell - ", cellid))

    if (!is.null(chrfilt)) {
        message(paste0("Filtering for chromosomes: ", paste0(chrfilt, collapse = ",")))
        CNbins <- dplyr::filter(CNbins, chr %in% chrfilt)
    }

    pl <- CNbins %>%
        dplyr::filter(cell_id == cellid) %>%
        plottinglist(., xaxis_order = xaxis_order, maxCN = maxCN, tickwidth = tickwidth, chrstart = chrstart, chrend = chrend)

    if (raster == TRUE) {
        if (!requireNamespace("ggrastr", quietly = TRUE)) {
            stop("Package \"ggrastr\" needed for this function to work. Please install it.",
                call. = FALSE
            )
        }
        gCN <- pl$CNbins %>%
            dplyr::mutate(state = ifelse(state >= 11, "11+", paste0(state))) %>%
            dplyr::mutate(state = factor(paste0(state), levels = c(paste0(seq(0, 10, 1)), "11+"))) %>%
            ggplot2::ggplot(ggplot2::aes(x = idx, y = copy)) +
            ggplot2::geom_vline(xintercept = pl$chrbreaks, col = "grey90", alpha = 0.75) +
            ggrastr::geom_point_rast(ggplot2::aes_string(col = statecol), size = pointsize, alpha = alphaval, shape = shape) +
            ggplot2::scale_color_manual(
                name = "Copy number",
                breaks = names(statecolpal),
                labels = names(statecolpal),
                values = statecolpal,
                drop = FALSE
            ) +
            ggplot2::theme(
                axis.title.y = ggplot2::element_blank(),
                axis.text.y = ggplot2::element_blank(),
                axis.ticks.y = ggplot2::element_blank(),
                legend.position = "none"
            ) +
            ggplot2::scale_x_continuous(breaks = pl$chrticks, labels = pl$chrlabels, expand = c(0, 0), limits = c(pl$minidx, pl$maxidx), guide = ggplot2::guide_axis(check.overlap = TRUE)) +
            ggplot2::scale_y_continuous(breaks = ybreaks, limits = c(0, maxCN), trans = y_axis_trans) +
            ggplot2::xlab(xlab) +
            ggplot2::ylab("Copy Number") +
            cowplot::theme_cowplot(...) +
            ggplot2::guides(colour = ggplot2::guide_legend(
                ncol = 6, byrow = TRUE,
                override.aes = list(alpha = 1, size = 3, shape = 15)
            ))
        # +
        # ggplot2::theme(legend.title = ggplot2::element_blank(), legend.position = legend.position)
    } else {
        gCN <- pl$CNbins %>%
            dplyr::mutate(state = ifelse(state >= 11, "11+", paste0(state))) %>%
            dplyr::mutate(state = factor(paste0(state), levels = c(paste0(seq(0, 10, 1)), "11+"))) %>%
            ggplot2::ggplot(ggplot2::aes(x = idx, y = copy)) +
            ggplot2::geom_vline(xintercept = pl$chrbreaks, col = "grey90", alpha = 0.75) +
            ggplot2::geom_point(ggplot2::aes_string(col = statecol), size = pointsize, alpha = alphaval, shape = 16) +
            ggplot2::scale_color_manual(
                name = "Allele Specific CN",
                breaks = names(statecolpal),
                labels = names(statecolpal),
                values = statecolpal,
                drop = FALSE
            ) +
            ggplot2::theme(
                axis.title.y = ggplot2::element_blank(),
                axis.text.y = ggplot2::element_blank(),
                axis.ticks.y = ggplot2::element_blank(),
                legend.position = "none"
            ) +
            ggplot2::scale_x_continuous(breaks = pl$chrticks, labels = pl$chrlabels, expand = c(0, 0), limits = c(pl$minidx, pl$maxidx), guide = ggplot2::guide_axis(check.overlap = TRUE)) +
            ggplot2::scale_y_continuous(breaks = ybreaks, limits = c(0, maxCN), trans = y_axis_trans) +
            ggplot2::xlab(xlab) +
            ggplot2::ylab("Copy Number") +
            cowplot::theme_cowplot(...) +
            ggplot2::guides(colour = ggplot2::guide_legend(
                ncol = 6, byrow = TRUE,
                override.aes = list(alpha = 1, size = 3, shape = 15)
            ))
        # +
        # ggplot2::theme(legend.title = ggplot2::element_blank(), legend.position = legend.position)
    }

    if (!is.null(genes)) {
        gene_idx <- get_gene_idx(genes, chr = chrfilt)
        npoints <- dim(pl$CNbins)[1]
        gCN <- gCN +
            ggplot2::geom_vline(data = gene_idx, ggplot2::aes(xintercept = idx), lty = 2, size = 0.3) +
            ggrepel::geom_text_repel(data = gene_idx, ggplot2::aes(x = idx - npoints * adj, y = maxCN, label = ensembl_gene_symbol), col = "black", alpha = 0.75)
    }

    if (!is.null(annotateregions)) {
        datidx <- dplyr::inner_join(annotateregions, pl$bins %>% dplyr::select(chr, start, idx)) %>% dplyr::distinct(.)
        gCN <- gCN +
            ggplot2::geom_vline(data = datidx, ggplot2::aes(xintercept = idx), lty = 2, size = 0.3, alpha = 0.5)
    }

    if (!is.null(SV)) {
        svpl <- plottinglistSV(SV, chrfilt = chrfilt)
        bezdf <- get_bezier_df(svpl, pl, maxCN)
        bezdf <- bezdf %>%
            dplyr::filter((position_1 != position_2) | rearrangement_type == "foldback")
        gCN <- gCN +
            ggforce::geom_bezier(ggplot2::aes(x = idx, y = copy, group = id),
                alpha = 0.8,
                size = svwidth,
                col = as.vector(SV_colors["Foldback"]),
                data = bezdf %>% dplyr::filter(rearrangement_type == "foldback")
            ) +
            ggforce::geom_bezier(ggplot2::aes(x = idx, y = copy, group = id),
                alpha = svalpha,
                size = svwidth,
                col = "grey30",
                data = bezdf %>% dplyr::filter(rearrangement_type != "foldback")
            )
    }

    if (returnlist == TRUE) {
        gCN <- list(CN = gCN, plist = pl)
    }


    return(gCN)
}
