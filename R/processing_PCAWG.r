#' @export
processing_PCAWG <- function(model_variables,
                             folder_name,
                             copynumber_PCAWG,
                             focal_threshold = 0.5,
                             height = 1000,
                             width = 2000) {
    #---------------------Check and correct model variables if necessary
    model_variables <- CHECK_model_variables(model_variables)
    #-------------------------------Prepare ingredients for the analysis
    CN_bin_length <- as.double(model_variables$general_variables$Value[which(model_variables$general_variables$Variable == "size_CN_block_DNA")])
    list_donors <- unique(copynumber_PCAWG$donor_unique_id)
    list_arms <- model_variables$chromosome_arm_library$Arm_ID
    driver_library <- model_variables$driver_library
    cn_info <- model_variables$cn_info
    cn_info$Centromere_bp <- cn_info$Centromere_location * CN_bin_length
    #   Find genome coordinate
    model_variables$general_variables$Value[which(model_variables$general_variables$Variable == "prob_CN_whole_genome_duplication")] <- 0
    model_variables$general_variables$Value[which(model_variables$general_variables$Variable == "prob_CN_missegregation")] <- 0
    model_variables$general_variables$Value[which(model_variables$general_variables$Variable == "prob_CN_chrom_arm_missegregation")] <- 0
    model_variables$general_variables$Value[which(model_variables$general_variables$Variable == "prob_CN_focal_amplification")] <- 0
    model_variables$general_variables$Value[which(model_variables$general_variables$Variable == "prob_CN_focal_deletion")] <- 0
    model_variables$general_variables$Value[which(model_variables$general_variables$Variable == "prob_CN_cnloh_interstitial")] <- 0
    model_variables$general_variables$Value[which(model_variables$general_variables$Variable == "prob_CN_cnloh_terminal")] <- 0
    model_variables$general_variables$Value[which(model_variables$general_variables$Variable == "rate_driver")] <- 0
    copynumber_sims <- simulator_full_program(
        model = model_variables, model_prefix = "TEST", n_simulations = 1, stage_final = 2,
        save_simulation = FALSE, report_progress = FALSE,
        output_variables = c("all_sample_genotype", "sample_cell_ID", "sample_genotype_unique", "sample_genotype_unique_profile")
    )
    simulation <- copynumber_sims[[1]]
    sample_genotype_unique_profile <- simulation$sample$sample_genotype_unique_profile
    CNbins_iteration <- sample_genotype_unique_profile[[1]]
    copynumber_coordinates <- CNbins_iteration[, 1:3]
    copynumber_coordinates$width <- copynumber_coordinates$end - copynumber_coordinates$start + 1
    #----------------------------------------Get bin-level CN from PCAWG
    #   Open workspace folder
    dir.create(folder_name)
    #   Convert PCAWG data into bin-level CN of equal lengths
    bin_PCAWG <- gainloss_PCAWG(
        copynumber_PCAWG,
        copynumber_coordinates,
        ploidy_normalization = TRUE,
        use_rbindlist = TRUE,
    )
    bin_PCAWG <- bin_PCAWG$copynumber_data
    #   Save CN of each arm of each donor in individual files
    ls_subsample <- c()
    ls_arm_ID <- c()
    ls_chrom <- c()
    for (donor in 1:length(list_donors)) {
        donor_ID <- list_donors[donor]
        donor_CN <- bin_PCAWG[[donor_ID]]
        donor_ID <- sub("::", "-", donor_ID)
        for (arm in 1:length(list_arms)) {
            arm_ID <- list_arms[arm]
            chrom <- substr(arm_ID, 1, nchar(arm_ID) - 1)
            arm <- substr(arm_ID, nchar(arm_ID), nchar(arm_ID))
            centromere <- cn_info$Centromere_bp[which(cn_info$Chromosome == chrom)]
            #   Find appropriate CN substring
            if (arm == "p") {
                arm_CN <- donor_CN[which(bin_PCAWG$chr == chrom & bin_PCAWG$end <= centromere)]
            } else if (arm == "q") {
                arm_CN <- donor_CN[which(bin_PCAWG$chr == chrom & bin_PCAWG$start > centromere)]
            }
            arm_CN_diploid <- rep(2, length(arm_CN))
            #   Save CN of this donor/arm
            ID <- paste(donor_ID, "-", arm_ID, sep = "")
            filename <- file.path(folder_name, paste(ID, ".csv", sep = ""))
            arm_df <- data.frame(From = arm_CN_diploid, To = arm_CN)
            write.csv(arm_df, filename, row.names = FALSE)
            ls_subsample <- c(ls_subsample, ID)
            ls_arm_ID <- c(ls_arm_ID, arm_ID)
            ls_chrom <- c(ls_chrom, chrom)
        }
    }
    #   Save list of donor ID's
    filename <- file.path(folder_name, "ID_list.csv")
    df_subsample <- data.frame(ID = ls_subsample)
    write.csv(df_subsample, filename, row.names = FALSE)
    #--------------------------------Get list of focal events from PCAWG
    ####################################################################
    ################ REWRITE THIS PART TO RUN FROM R!!! ################
    ####################################################################
    #   IN TERMINAL: RUN
    #   python
    #   from WCND import *
    #   WCND_arm('OV-AU')
    ####################################################################
    #-------------------------------------------Get data of focal events
    #   Find bin length of each chromosome arm
    list_arms_length <- rep(0, length(list_arms))
    for (i in 1:length(list_arms_length)) {
        chrom <- substr(list_arms[i], 1, nchar(list_arms[i]) - 1)
        pos <- which(cn_info$Chromosome == chrom)
        arm <- substr(list_arms[i], nchar(list_arms[i]), nchar(list_arms[i]))
        if (arm == "p") {
            list_arms_length[i] <- cn_info$Centromere_location[pos]
        } else if (arm == "q") {
            list_arms_length[i] <- cn_info$Bin_count[pos] - cn_info$Centromere_location[pos]
        }
    }
    df_chrom_arm <- data.frame(arm = list_arms, bin_count = list_arms_length)
    #   Record lengths of focal events and which drivers are hit
    driver_library$gain <- 0
    driver_library$loss <- 0
    ls_foc_amp_len <- c()
    ls_foc_del_len <- c()
    ls_foc_amp_len_ratio <- c()
    ls_foc_del_len_ratio <- c()
    ls_foc_amp_arm_len <- c()
    ls_foc_del_arm_len <- c()
    for (i in 1:length(ls_subsample)) {
        ID <- ls_subsample[i]
        chrom <- ls_chrom[i]
        arm <- ls_arm_ID[i]
        arm_len <- df_chrom_arm$bin_count[which(df_chrom_arm$arm == arm)]
        filename <- file.path(folder_name, paste(ID, "_focal_events.csv", sep = ""))
        ls_focal_events <- read.csv(file = filename)
        if (nrow(ls_focal_events) == 0) {
            next
        }
        #   Compute focal event lengths
        ls_focal_events$len <- ls_focal_events$end - ls_focal_events$start + 1
        #   Remove focal events that are too large
        vec_del <- which(ls_focal_events$len / arm_len > focal_threshold)
        if (length(vec_del) > 0) {
            ls_focal_events <- ls_focal_events[-vec_del, ]
        }
        #   Save focal event lengths
        ls_foc_amp_len <- c(ls_foc_amp_len, ls_focal_events$len[which(ls_focal_events$sign == "+")])
        ls_foc_del_len <- c(ls_foc_del_len, ls_focal_events$len[which(ls_focal_events$sign == "-")])

        ls_foc_amp_len_ratio <- c(ls_foc_amp_len_ratio, (ls_focal_events$len[which(ls_focal_events$sign == "+")]) / arm_len)
        ls_foc_del_len_ratio <- c(ls_foc_del_len_ratio, (ls_focal_events$len[which(ls_focal_events$sign == "-")]) / arm_len)

        ls_foc_amp_arm_len <- c(ls_foc_amp_arm_len, rep(arm_len, length(which(ls_focal_events$sign == "+"))))
        ls_foc_del_arm_len <- c(ls_foc_del_arm_len, rep(arm_len, length(which(ls_focal_events$sign == "-"))))
        #   Check if any driver is hit
        gene_foc_amp <- c()
        gene_foc_del <- c()
        for (foc in 1:nrow(ls_focal_events)) {
            start <- ls_focal_events$start[foc] + 1
            end <- ls_focal_events$end[foc] + 1
            sign <- ls_focal_events$sign[foc]
            loc <- which((driver_library$Chromosome == chrom) & (driver_library$Bin >= start) & (driver_library$Bin <= end))
            if (length(loc) > 0) {
                if (sign == "+") {
                    gene_foc_amp <- c(gene_foc_amp, loc)
                } else if (sign == "-") {
                    gene_foc_del <- c(gene_foc_del, loc)
                }
            }
        }
        driver_library$gain[unique(gene_foc_amp)] <- driver_library$gain[unique(gene_foc_amp)] + 1
        driver_library$loss[unique(gene_foc_del)] <- driver_library$loss[unique(gene_foc_del)] + 1
    }
    vec_del <- which(driver_library$gain == 0 & driver_library$loss == 0)
    if (length(vec_del) > 0) {
        driver_library <- driver_library[-vec_del, ]
    }
    driver_library$gain_ratio <- driver_library$gain / length(list_donors)
    driver_library$loss_ratio <- driver_library$loss / length(list_donors)
    driver_library$Classification <- driver_library$Gene_role
    #------------------------------------Fitting for focal event lengths
    #   Fit focal event lengths for gains/losses
    amp_estimate <- fitdist(ls_foc_amp_len_ratio, "beta")$estimate
    prob_amp_1 <- amp_estimate[1]
    prob_amp_2 <- amp_estimate[2]
    del_estimate <- fitdist(ls_foc_del_len_ratio, "beta")$estimate
    prob_del_1 <- del_estimate[1]
    prob_del_2 <- del_estimate[2]
    #   Compute data for comparison with observed focal event lengths
    step <- 0.01
    fit_amp_x <- seq(0, 1, by = step)
    fit_amp_y <- dbeta(fit_amp_x, prob_amp_1, prob_amp_2)
    df_fit_amp <- data.frame(x = fit_amp_x, y = fit_amp_y)
    fit_del_x <- seq(0, 1, by = step)
    fit_del_y <- dbeta(fit_del_x, prob_del_1, prob_del_2)
    df_fit_del <- data.frame(x = fit_del_x, y = fit_del_y)
    #-----------------------------------------------------------Plotting
    #   Plot distribution of driver genes being focally amplified/deleted
    filename <- paste(folder_name, "_driver_genes.png", sep = "")
    jpeg(filename, width = height, height = height)
    p <- ggplot(driver_library) +
        aes(x = gain_ratio, y = loss_ratio, color = Classification, label = Gene_ID) +
        geom_point(size = 10) +
        geom_label_repel(size = 10) +
        geom_abline() +
        theme(panel.background = element_rect(fill = "white", colour = "grey50")) +
        theme(text = element_text(size = 20), legend.position = "bottom") +
        scale_x_continuous(limits = c(0, 1), expand = c(0, 0)) +
        scale_y_continuous(limits = c(0, 1), expand = c(0, 0)) +
        ylab("Amplifications") +
        xlab("Deletions")
    print(p)
    dev.off()
    #   Plot distribution of focal length vs focal length ratio in chromosome arm
    filename <- paste(folder_name, "_focal_length_vs_focal_ratio.png", sep = "")
    jpeg(filename, width = height, height = height)
    p <- ggplot(
        data.frame(
            focal_length = c(ls_foc_amp_len, ls_foc_del_len),
            focal_ratio = c(ls_foc_amp_arm_len, ls_foc_del_arm_len),
            type = c(rep("AMP", length(ls_foc_amp_len)), rep("DEL", length(ls_foc_del_len)))
        ),
        aes(x = focal_length, y = focal_ratio, color = type)
    ) +
        geom_point() +
        # geom_abline() +
        geom_abline(intercept = 0, slope = 1 / focal_threshold, colour = "red") +
        theme(panel.background = element_rect(fill = "white", colour = "grey50")) +
        theme(text = element_text(size = 20)) +
        # scale_x_continuous(expand = c(1, 1)) +
        # scale_y_continuous(expand = c(1, 1)) +
        ylab("Chromosome arm length") +
        xlab("Focal event length")
    print(p)
    dev.off()
    #   Plot distributions of focal event lengths
    filename <- paste(folder_name, "_focal_length_distribution.png", sep = "")
    jpeg(filename, width = width, height = height)
    p1 <- ggplot(
        data.frame(
            length = ls_foc_amp_len_ratio
        ), aes(x = length)
    ) +
        geom_density(alpha = 0.3, fill = "#E34A33", color = "#E34A33") +
        geom_line(data = df_fit_amp, aes(x = x, y = y), color = "#E34A33", inherit.aes = FALSE, size = 2) +
        theme(panel.background = element_rect(fill = "white", colour = "grey50")) +
        theme(text = element_text(size = 20)) +
        scale_x_continuous(limits = c(min(ls_foc_amp_len_ratio), 1), expand = c(0, 0)) +
        scale_y_continuous(expand = c(0, 0)) +
        ylab("Amplification") +
        xlab("Ratio of chromosome arm")
    p2 <- ggplot(
        data.frame(
            length = ls_foc_del_len_ratio
        ),
        aes(x = length)
    ) +
        geom_density(alpha = 0.3, fill = "#3182BD", color = "#3182BD") +
        geom_line(data = df_fit_del, aes(x = x, y = y), color = "#3182BD", inherit.aes = FALSE, size = 2) +
        theme(panel.background = element_rect(fill = "white", colour = "grey50")) +
        theme(text = element_text(size = 20)) +
        scale_x_continuous(limits = c(min(ls_foc_del_len_ratio), 1), expand = c(0, 0)) +
        scale_y_continuous(expand = c(0, 0)) +
        ylab("Deletion") +
        xlab("Ratio of chromosome arm")
    p <- grid.arrange(p1, p2, nrow = 1)
    print(p)
    dev.off()
    #-------------------------------Print suggestions for better fitting
    cat("===========================================================\n")
    cat("============== SUGGESTION FOR BETTER FITTING ==============\n")
    cat("===========================================================\n")
    cat("FOCAL AMPLIFICATION/DELETION:\n\n")
    cat(paste("Focal amplification: focal length / arm length ~ Beta(shape1=", prob_amp_1, ",shape2=", prob_amp_2, ")\n\n", sep = ""))
    cat(paste("Focal deletion: focal length / arm length ~ Beta(shape1=", prob_del_1, ",shape2=", prob_del_2, ")\n", sep = ""))
    cat("===========================================================\n")
    cat("DRIVER GENE LIBRARY:\n\n")
    for (gene in 1:nrow(driver_library)) {
        cat(paste(driver_library$Gene_ID[gene], "\n", sep = ""))
    }
    cat("===========================================================\n")
    #---------------------Output better parameters for fitting the model
    driver_library <- driver_library[, !(names(driver_library) %in% c("gain", "loss", "gain_ratio", "loss_ratio", "Classification"))]

    output <- list()
    output$prob_CN_focal_amplification_length_shape_1 <- prob_amp_1
    output$prob_CN_focal_amplification_length_shape_2 <- prob_amp_2
    output$prob_CN_focal_deletion_length_shape_1 <- prob_del_1
    output$prob_CN_focal_deletion_length_shape_2 <- prob_del_2
    output$driver_library <- driver_library

    return(output)
}
