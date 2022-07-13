#' @export
plot_gainloss_PCAWG <- function(copynumber_PCAWG) {
    #----------------------------------------------Rework the data frame
    #   Conform data frame into format required by signals
    CNbins <- copynumber_PCAWG[c("chromosome", "start", "end", "total_cn", "minor_cn", "major_cn", "donor_unique_id")]

    # CNbins <- CNbins[which(CNbins$donor_unique_id == CNbins$donor_unique_id[1]), ]

    CNbins$state <- CNbins$total_cn
    CNbins <- CNbins[, c(1, 2, 3, 4, 8, 5, 6, 7)]
    colnames(CNbins) <- c("chr", "start", "end", "copy", "state", "Min", "Maj", "cell_id")
    #

    plotcol <- "state"
    fillna <- TRUE
    cutoff <- NULL
    maxf <- 1
    plotfrequency <- TRUE
    SV <- NULL


    copynumber <- createCNmatrix(CNbins,
        field = plotcol, wholegenome = TRUE,
        fillnaplot = fillna, centromere = FALSE
    )

    message("Normalizing ploidy for each cell to 2")
    copynumber <- normalize_cell_ploidy(copynumber)

    top_annotation <- make_top_annotation_gain(copynumber,
        cutoff = cutoff, maxf = maxf,
        plotfrequency = plotfrequency, plotcol = plotcol, SV = SV
    )

    # ComplexHeatmap::draw(top_annotation)

    print(CNbins)
    print(copynumber)

    return(top_annotation)
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

normalize_cell_ploidy <- function(copynumber) {
    cell_ids <- colnames(copynumber)
    cell_ids <- cell_ids[!(cell_ids %in% c("chr", "start", "end", "width"))]

    for (cell_id in cell_ids) {
        state_mode <- calc_state_mode(copynumber[[cell_id]])
        copynumber[[cell_id]] <- as.integer(ceiling(
            copynumber[[cell_id]] / (state_mode / 2)
        ))
        copynumber[[cell_id]][copynumber[[cell_id]] > 11] <- 11
    }
    return(copynumber)
}

make_top_annotation_gain <- function(copynumber,
                                     plotcol = "state",
                                     plotfrequency = FALSE,
                                     cutoff = NULL,
                                     maxf = NULL,
                                     SV = NULL) {
    ncells <- nrow(copynumber)

    print(ncells)

    f1 <- colSums(copynumber > cutoff, na.rm = TRUE) / ncells
    f2 <- -colSums(copynumber < cutoff, na.rm = TRUE) / ncells
    if (is.null(maxf)) {
        maxf <- ceiling(max(max(f1, max(abs(f2)))) / 0.1) * 0.1
        if (maxf < 0.01) {
            maxf <- 0.01
        }
    }
    ha2 <- ComplexHeatmap::columnAnnotation(
        dist2 = ComplexHeatmap::anno_barplot(
            f1,
            bar_width = 1,
            gp = grid::gpar(col = "#E34A33", fill = "#E34A33"),
            axis_param = list(
                at = c(round(maxf / 2, 2), maxf),
                labels = c(paste0(round(maxf / 2, 2)), paste0(maxf))
            ),
            ylim = c(0, maxf),
            border = FALSE,
        ),
        dist3 = ComplexHeatmap::anno_barplot(
            f2,
            bar_width = 1,
            gp = grid::gpar(col = "#3182BD", fill = "#3182BD"),
            axis_param = list(
                at = c(0.0, -round(maxf / 2, 2), -maxf),
                labels = c("0", paste0(round(maxf / 2, 2)), paste0(maxf))
            ),
            ylim = c(-maxf, 0),
            border = FALSE,
        ),
        show_annotation_name = FALSE,
        height = grid::unit(1.4, "cm")
    )

    ComplexHeatmap::draw(ha2)

    return(ha2)
}
