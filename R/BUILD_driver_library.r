#' Compute and add the selection rates for driver alleles based on the choice of selection model
#'
#' @description
#' `BUILD_driver_library` returns an updated version of `model_variables` that includes information about the chosen selection model and the resulting selection rates for different driver alleles based on whether they are classified as TSGs or Oncogenes.
#' 
#' @param model_variables The object returned by `BUILD_general_variables`.
#' @param table_arm_selection_rates A dataframe consisting of columns 'Arm_ID', 'Chromosome', 'Bin_start', 'Bin_end', 's_rate'. 'Arm_ID' refers to the name of each chromosome arm eg. 1p, 1q, etc. 'Chromosome' refers to the chromosome number. 'Bin_start' refers to the ???. 'Bin_end' refers to the ???. 's_rate' refers to the selective value of the chromosome arm as a float.
#' @param table_gene_selection_rates
#' @param vec_driver_genes A list of driver genes. Only needed if `selection_model` is 'ancient'.
#' @param vec_driver_role A list corresponding to `vec_driver_genes` classifying each driver gene as a TSG or Oncogene. Only needed if `selection_model` is 'ancient'.
#' @param vec_chromosome ???
#' @param vec_bin ???
#' @param vec_driver_s A list corresponding to `vec_driver_genes` specifying the selection rate of the driver gene. Only needed if `selection_model` is 'ancient'.
#' @param vec_id ???
#' @param vec_start ???
#' @param vec_end ???
#' @param vec_arm_s ???
#'
#' @examples
#' 
#' model_variables_base <- BUILD_driver_library(
#'    model_variables = model_variables_base,
#'    table_arm_selection_rates = data.frame(Arm_ID = arm_id, Chromosome = arm_chromosome, Bin_start = arm_start, Bin_end = arm_end, s_rate = arm_s))
#'
#' @export
BUILD_driver_library <- function(model_variables = list(),
                                 table_arm_selection_rates = list(),
                                 table_gene_selection_rates = list(),
                                 vec_driver_genes = c(),
                                 vec_driver_role = c(),
                                 vec_chromosome = -1,
                                 vec_bin = -1,
                                 vec_driver_s = c(),
                                 vec_id = c(),
                                 vec_start = c(),
                                 vec_end = c(),
                                 vec_arm_s = c()) {
    #------------------------------------Input choice of selection model
    TABLE_SELECTION <- model_variables$selection_model
    selection_model <- TABLE_SELECTION$Value[which(TABLE_SELECTION$Variable == "selection_model")]
    #---------------------Compute selection rates for each driver allele
    #-----------------------------according to choice of selection model
    if ((selection_model == "chrom-arm-selection") | (selection_model == "WGD-chrom-arm-selection")) {
        #---Build the arm driver library
        TABLE_CANCER_ARMS <- table_arm_selection_rates
        #--------------------------------Output the model variable files
        model_variables$chromosome_arm_library <- TABLE_CANCER_ARMS
        model_variables$selection_model <- TABLE_SELECTION
    } else if (selection_model == "driver-gene-selection") {
        #---Build the gene driver library
        TABLE_CANCER_GENES <- table_gene_selection_rates
        TABLE_CANCER_GENES$s_rate_WT <- 0
        TABLE_CANCER_GENES$s_rate_MUT <- 0
        #   Compute selection rates for TSGs
        list_TSG <- which(TABLE_CANCER_GENES$Gene_role == "TSG")
        if (length(list_TSG) > 0) {
            for (driver in 1:length(list_TSG)) {
                row <- list_TSG[driver]
                driver_sel_rate <- TABLE_CANCER_GENES$s_rate[row]
                TABLE_CANCER_GENES$s_rate_WT[row] <- 1 / driver_sel_rate
                TABLE_CANCER_GENES$s_rate_MUT[row] <- 1
            }
        }
        #   Compute selection rates for ONCOGENEs
        list_ONCOGENE <- which(TABLE_CANCER_GENES$Gene_role == "ONCOGENE")
        if (length(list_ONCOGENE) > 0) {
            for (driver in 1:length(list_ONCOGENE)) {
                row <- list_ONCOGENE[driver]
                driver_sel_rate <- TABLE_CANCER_GENES$s_rate[row]
                TABLE_CANCER_GENES$s_rate_WT[row] <- driver_sel_rate
                TABLE_CANCER_GENES$s_rate_MUT[row] <- driver_sel_rate^2
            }
        }
        #--------------------------------Output the model variable files
        model_variables$driver_library <- TABLE_CANCER_GENES
        model_variables$selection_model <- TABLE_SELECTION
    } else if (selection_model == "chrom-arm-and-driver-gene-selection") {
        #---Build the arm driver library
        TABLE_CANCER_ARMS <- table_arm_selection_rates
        #---Build the gene driver library
        TABLE_CANCER_GENES <- table_gene_selection_rates
        if (nrow(TABLE_CANCER_GENES) > 0) {
            TABLE_CANCER_GENES$s_rate_WT <- 0
            TABLE_CANCER_GENES$s_rate_MUT <- 0
            #   Compute selection rates for TSGs
            list_TSG <- which(TABLE_CANCER_GENES$Gene_role == "TSG")
            if (length(list_TSG) > 0) {
                for (driver in 1:length(list_TSG)) {
                    row <- list_TSG[driver]
                    driver_sel_rate <- TABLE_CANCER_GENES$s_rate[row]
                    TABLE_CANCER_GENES$s_rate_WT[row] <- 1 / driver_sel_rate
                    TABLE_CANCER_GENES$s_rate_MUT[row] <- 1
                }
            }
            #   Compute selection rates for ONCOGENEs
            list_ONCOGENE <- which(TABLE_CANCER_GENES$Gene_role == "ONCOGENE")
            if (length(list_ONCOGENE) > 0) {
                for (driver in 1:length(list_ONCOGENE)) {
                    row <- list_ONCOGENE[driver]
                    driver_sel_rate <- TABLE_CANCER_GENES$s_rate[row]
                    TABLE_CANCER_GENES$s_rate_WT[row] <- driver_sel_rate
                    TABLE_CANCER_GENES$s_rate_MUT[row] <- driver_sel_rate^2
                }
            }
        }
        #--------------------------------Output the model variable files
        model_variables$chromosome_arm_library <- TABLE_CANCER_ARMS
        model_variables$driver_library <- TABLE_CANCER_GENES
        model_variables$selection_model <- TABLE_SELECTION
    } else if (selection_model == "ancient") {
        #---Input the Cancer Gene Census
        DATA_cancer_gene_census <- read.csv(file = system.file("extdata", "cancer_gene_census.csv", package = "CancerSimulator"))
        #---Build the driver library
        columns <- c("Gene_ID", "Gene_role", "s_rate")
        TABLE_CANCER_GENES <- data.frame(vec_driver_genes, vec_driver_role, vec_driver_s)
        colnames(TABLE_CANCER_GENES) <- columns
        if (any(vec_chromosome == -1)) {
            TABLE_CANCER_GENES$Chromosome <- -1
        } else {
            TABLE_CANCER_GENES$Chromosome <- vec_chromosome
        }
        if (any(vec_bin == -1)) {
            TABLE_CANCER_GENES$Bin <- -1
        } else {
            TABLE_CANCER_GENES$Bin <- vec_bin
        }
        #---Supplement driver library with gene locations
        if ((any(TABLE_CANCER_GENES$Chromosome == -1) == TRUE) | (any(TABLE_CANCER_GENES$Bin == -1) == TRUE)) {
            #       Get the CN bin length
            size_CN_block_DNA <- as.numeric(model_variables$general_variables$Value[model_variables$general_variables$Variable == "size_CN_block_DNA"])
            #       Find the chromosome and bin for each driver
            for (gene in 1:length(vec_driver_genes)) {
                Gene_ID <- vec_driver_genes[gene]
                loc <- which(DATA_cancer_gene_census$Gene.Symbol == Gene_ID)
                Gene_address <- DATA_cancer_gene_census$Genome.Location[loc]
                Gene_chromosome <- as.numeric(sub(":.*", "", Gene_address))
                Gene_bin <- round(as.numeric(sub("-.*", "", sub(".*:", "", Gene_address))) / size_CN_block_DNA)
                TABLE_CANCER_GENES$Chromosome[gene] <- Gene_chromosome
                TABLE_CANCER_GENES$Bin[gene] <- Gene_bin
            }
        }
        #---Supplement driver library with allele selection rates
        #   Count the number of TSGs and ONCOGENEs
        count_TSG <- sum(TABLE_CANCER_GENES$Gene_role == "TSG")
        count_ONCOGENE <- sum(TABLE_CANCER_GENES$Gene_role == "ONCOGENE")
        #---Compute selection rates for TSGs
        list_TSG <- which(TABLE_CANCER_GENES$Gene_role == "TSG")
        s_normalization <- 1
        if (length(list_TSG) > 0) {
            for (driver in 1:length(list_TSG)) {
                row <- list_TSG[driver]
                #       Get its selection strength
                driver_sel_rate <- vec_driver_s[row]
                #       Compute its selection rate for WT and MUT alleles
                TABLE_CANCER_GENES$s_rate_WT[row] <- 1 / (1 + driver_sel_rate)
                TABLE_CANCER_GENES$s_rate_MUT[row] <- 1
                #       Update normalizer for selection rate
                s_normalization <- s_normalization * (1 + driver_sel_rate)
            }
        }
        # #---Compute selection rates for ONCOGENEs
        # list_ONCOGENE <- which(TABLE_CANCER_GENES$Gene_role == "ONCOGENE")
        # if (length(list_ONCOGENE) > 0) {
        #     for (driver in 1:length(list_ONCOGENE)) {
        #         row <- list_ONCOGENE[driver]
        #         #       Get its selection strength
        #         driver_sel_rate <- vec_driver_s[row]
        #         #       Compute its selection rate for WT and MUT alleles
        #         TABLE_CANCER_GENES$s_rate_WT[row] <- (1 + driver_sel_rate)
        #         TABLE_CANCER_GENES$s_rate_MUT[row] <- (1 + driver_sel_rate)^2
        #     }
        # }
        #---Compute selection rates for ONCOGENEs
        list_ONCOGENE <- which(TABLE_CANCER_GENES$Gene_role == "ONCOGENE")
        if (length(list_ONCOGENE) > 0) {
            for (driver in 1:length(list_ONCOGENE)) {
                row <- list_ONCOGENE[driver]
                #       Get its selection strength
                driver_sel_rate <- vec_driver_s[row]
                #       Compute its selection rate for WT and MUT alleles
                TABLE_CANCER_GENES$s_rate_WT[row] <- 1
                TABLE_CANCER_GENES$s_rate_MUT[row] <- (1 + driver_sel_rate)
            }
        }
        #---Update selection model with normalization constant
        N_row <- nrow(TABLE_SELECTION)
        N_row <- N_row + 1
        TABLE_SELECTION[N_row, ] <- c("s_normalization", s_normalization, "", "Normalization constant")
        #--------------------------------Output the model variable files
        model_variables$driver_library <- TABLE_CANCER_GENES
        model_variables$selection_model <- TABLE_SELECTION
    } else if (selection_model == "old") {
        #---Input the Cancer Gene Census
        DATA_cancer_gene_census <- read.csv(file = system.file("extdata", "cancer_gene_census.csv", package = "CancerSimulator"))
        #---Build the driver library
        columns <- c("Gene_ID", "Gene_role", "s_rate")
        TABLE_CANCER_GENES <- data.frame(vec_driver_genes, vec_driver_role, vec_driver_s)
        colnames(TABLE_CANCER_GENES) <- columns
        if (any(vec_chromosome == -1)) {
            TABLE_CANCER_GENES$Chromosome <- -1
        } else {
            TABLE_CANCER_GENES$Chromosome <- vec_chromosome
        }
        if (any(vec_bin == -1)) {
            TABLE_CANCER_GENES$Bin <- -1
        } else {
            TABLE_CANCER_GENES$Bin <- vec_bin
        }
        #---Supplement driver library with gene locations
        if ((any(TABLE_CANCER_GENES$Chromosome == -1) == TRUE) | (any(TABLE_CANCER_GENES$Bin == -1) == TRUE)) {
            #       Get the CN bin length
            size_CN_block_DNA <- as.numeric(model_variables$general_variables$Value[model_variables$general_variables$Variable == "size_CN_block_DNA"])
            #       Find the chromosome and bin for each driver
            for (gene in 1:length(vec_driver_genes)) {
                Gene_ID <- vec_driver_genes[gene]
                loc <- which(DATA_cancer_gene_census$Gene.Symbol == Gene_ID)
                Gene_address <- DATA_cancer_gene_census$Genome.Location[loc]
                Gene_chromosome <- as.numeric(sub(":.*", "", Gene_address))
                Gene_bin <- round(as.numeric(sub("-.*", "", sub(".*:", "", Gene_address))) / size_CN_block_DNA)
                TABLE_CANCER_GENES$Chromosome[gene] <- Gene_chromosome
                TABLE_CANCER_GENES$Bin[gene] <- Gene_bin
            }
        }
        #---Supplement driver library with allele selection rates
        # #   Count the number of TSGs and ONCOGENEs
        # count_TSG <- sum(TABLE_CANCER_GENES$Gene_role == "TSG")
        # count_ONCOGENE <- sum(TABLE_CANCER_GENES$Gene_role == "ONCOGENE")
        #---Compute selection rates for TSGs
        list_TSG <- which(TABLE_CANCER_GENES$Gene_role == "TSG")
        s_normalization <- 1
        if (length(list_TSG) > 0) {
            for (driver in 1:length(list_TSG)) {
                row <- list_TSG[driver]
                #       Get its selection strength
                driver_sel_rate <- vec_driver_s[row]
                #       Compute its selection rate for WT and MUT alleles
                TABLE_CANCER_GENES$s_rate_WT[row] <- 1 / (1 + driver_sel_rate)
                TABLE_CANCER_GENES$s_rate_MUT[row] <- 1
            }
        }
        # #---Compute selection rates for TSGs
        # list_TSG <- which(TABLE_CANCER_GENES$Gene_role == "TSG")
        # s_normalization <- 1
        # if (length(list_TSG) > 0) {
        #     for (driver in 1:length(list_TSG)) {
        #         row <- list_TSG[driver]
        #         #       Get its selection strength
        #         driver_sel_rate <- vec_driver_s[row]
        #         #       Compute its selection rate for WT and MUT alleles
        #         TABLE_CANCER_GENES$s_rate_WT[row] <- 1 / (1 + driver_sel_rate)
        #         TABLE_CANCER_GENES$s_rate_MUT[row] <- 1
        #         #       Update normalizer for selection rate
        #         s_normalization <- s_normalization * (1 + driver_sel_rate)
        #     }
        # }
        #---Compute selection rates for ONCOGENEs
        list_ONCOGENE <- which(TABLE_CANCER_GENES$Gene_role == "ONCOGENE")
        if (length(list_ONCOGENE) > 0) {
            for (driver in 1:length(list_ONCOGENE)) {
                row <- list_ONCOGENE[driver]
                #       Get its selection strength
                driver_sel_rate <- vec_driver_s[row]
                #       Compute its selection rate for WT and MUT alleles
                TABLE_CANCER_GENES$s_rate_WT[row] <- (1 + driver_sel_rate)
                TABLE_CANCER_GENES$s_rate_MUT[row] <- (1 + driver_sel_rate)^2
            }
        }
        # #---Compute selection rates for ONCOGENEs
        # list_ONCOGENE <- which(TABLE_CANCER_GENES$Gene_role == "ONCOGENE")
        # if (length(list_ONCOGENE) > 0) {
        #     for (driver in 1:length(list_ONCOGENE)) {
        #         row <- list_ONCOGENE[driver]
        #         #       Get its selection strength
        #         driver_sel_rate <- vec_driver_s[row]
        #         #       Compute its selection rate for WT and MUT alleles
        #         TABLE_CANCER_GENES$s_rate_WT[row] <- 1
        #         TABLE_CANCER_GENES$s_rate_MUT[row] <- (1 + driver_sel_rate)
        #     }
        # }
        # #---Update selection model with normalization constant
        # N_row <- nrow(TABLE_SELECTION)
        # N_row <- N_row + 1
        # TABLE_SELECTION[N_row, ] <- c("s_normalization", s_normalization, "", "Normalization constant")
        #--------------------------------Output the model variable files
        model_variables$driver_library <- TABLE_CANCER_GENES
        model_variables$selection_model <- TABLE_SELECTION
    }
    return(model_variables)
}
