BUILD_driver_library <- function(model_variables = list(),
                                 vec_driver_genes = c(),
                                 vec_driver_role = c(),
                                 vec_chromosome = -1,
                                 vec_bin = -1,
                                 vec_driver_s = c()) {
    #-------------------------------------------Input the Cancer Gene Census
    DATA_cancer_gene_census <- read.csv("../data/cancer_gene_census.csv")
    #-----------------------------------------------Build the driver library
    columns <- c("Gene_ID", "Gene_role", "Selective_strength")
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
    #--------------------------Supplement driver library with gene locations
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
    #------------------Supplement driver library with allele selection rates
    #   Count the number of TSGs and ONCOGENEs
    count_TSG <- sum(TABLE_CANCER_GENES$Gene_role == "TSG")
    count_ONCOGENE <- sum(TABLE_CANCER_GENES$Gene_role == "ONCOGENE")
    #---Compute selection rates for TSGs
    list_TSG <- which(TABLE_CANCER_GENES$Gene_role == "TSG")
    s_normalization <- 1
    if (length(list_TSG)>0){
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


    TABLE_VARIABLES <- model_variables$general_variables
    N_row <- nrow(TABLE_VARIABLES)
    N_row <- N_row + 1
    TABLE_VARIABLES[N_row, ] <- c("s_normalization", s_normalization, "", "Normalization constant")
    model_variables$general_variables <- TABLE_VARIABLES


    #---Compute selection rates for ONCOGENEs
    list_ONCOGENE <- which(TABLE_CANCER_GENES$Gene_role == "ONCOGENE")
    if (length(list_ONCOGENE)>0){
        for (driver in 1:length(list_ONCOGENE)) {
            row <- list_ONCOGENE[driver]
            #       Get its selection strength
            driver_sel_rate <- vec_driver_s[row]
            #       Compute its selection rate for WT and MUT alleles
            TABLE_CANCER_GENES$s_rate_WT[row] <- 1
            TABLE_CANCER_GENES$s_rate_MUT[row] <- (1 + driver_sel_rate)
        }
    }
    #----------------------------------------Output the model variable files
    model_variables$driver_library <- TABLE_CANCER_GENES
    return(model_variables)
}
