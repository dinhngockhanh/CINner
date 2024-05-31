cell_lifespan <- 30
T_0 <- list(0, "year")
T_end <- list(80, "year")
Table_sample <- data.frame(Sample_ID = c("SA01"), Cell_count = c(Inf), Age_sample = c(80))
T_tau_step <- cell_lifespan / 2
CN_bin_length <- 500000
selection_model <- "chrom-arm-selection"
#------------------------------------------------------CNA PROBABILITIES
prob_CN_whole_genome_duplication <- 0e-4
prob_CN_missegregation <- 5e-5
prob_CN_chrom_arm_missegregation <- 0e-5
prob_CN_focal_amplification <- 0e-5
prob_CN_focal_deletion <- 0e-5
prob_CN_cnloh_interstitial <- 0e-5
prob_CN_cnloh_terminal <- 0e-5
model_CN_focal_amplification_length <- "beta"
model_CN_focal_deletion_length <- "beta"
prob_CN_focal_amplification_length_shape_1 <- 0.758304780825031
prob_CN_focal_amplification_length_shape_2 <- 5.33873409782625
prob_CN_focal_deletion_length_shape_1 <- 0.814054548726361
prob_CN_focal_deletion_length_shape_2 <- 6.16614890284825
prob_CN_cnloh_interstitial_length <- 0.005
prob_CN_cnloh_terminal_length <- 0.005
rate_driver <- 0
rate_passenger <- 1e-11
#---------------------------------------------------VIABILITY THRESHOLDS
bound_driver <- Inf
bound_average_ploidy <- Inf
bound_homozygosity <- 0
bound_maximum_CN <- Inf
bound_maximum_CN_normalized <- 4
#-----------------------------------------------------------------------
vec_time <- T_0[[1]]:T_end[[1]]
L <- 10000
t_0 <- 20
k <- 0.3
vec_cell_count <- L / (1 + exp(-k * (vec_time - t_0)))
table_population_dynamics <- cbind(vec_time, vec_cell_count)
gc <- read.csv(file = system.file("extdata", "gc_map_500kb.csv", package = "CancerSimulator"))
gc_slope <- 1.2
gc_int <- 0
sigma1 <- 0.1
num_reads <- 2e6
model_variables <- BUILD_general_variables(
    cell_lifespan = cell_lifespan,
    T_0 = T_0, T_end = T_end, T_tau_step = T_tau_step,
    Table_sample = Table_sample,
    CN_bin_length = CN_bin_length,
    prob_CN_whole_genome_duplication = prob_CN_whole_genome_duplication,
    prob_CN_missegregation = prob_CN_missegregation,
    prob_CN_chrom_arm_missegregation = prob_CN_chrom_arm_missegregation,
    prob_CN_focal_amplification = prob_CN_focal_amplification,
    prob_CN_focal_deletion = prob_CN_focal_deletion,
    prob_CN_cnloh_interstitial = prob_CN_cnloh_interstitial,
    prob_CN_cnloh_terminal = prob_CN_cnloh_terminal,
    model_CN_focal_amplification_length = model_CN_focal_amplification_length,
    model_CN_focal_deletion_length = model_CN_focal_deletion_length,
    prob_CN_focal_amplification_length_shape_1 = prob_CN_focal_amplification_length_shape_1,
    prob_CN_focal_amplification_length_shape_2 = prob_CN_focal_amplification_length_shape_2,
    prob_CN_focal_deletion_length_shape_1 = prob_CN_focal_deletion_length_shape_1,
    prob_CN_focal_deletion_length_shape_2 = prob_CN_focal_deletion_length_shape_2,
    prob_CN_cnloh_interstitial_length = prob_CN_cnloh_interstitial_length,
    prob_CN_cnloh_terminal_length = prob_CN_cnloh_terminal_length,
    rate_driver = rate_driver,
    rate_passenger = rate_passenger,
    selection_model = selection_model,
    bound_driver = bound_driver,
    bound_maximum_CN = bound_maximum_CN,
    bound_average_ploidy = bound_average_ploidy,
    bound_homozygosity = bound_homozygosity,
    table_population_dynamics = table_population_dynamics,
    gc = gc,
    gc_slope = gc_slope,
    gc_int = gc_int,
    sigma1 = sigma1,
    num_reads = num_reads
)
arm_id <- c(paste(model_variables$cn_info$Chromosome, "p", sep = ""), paste(model_variables$cn_info$Chromosome, "q", sep = ""))
arm_chromosome <- rep(model_variables$cn_info$Chromosome, 2)
arm_start <- c(rep(1, length(model_variables$cn_info$Chromosome)), model_variables$cn_info$Centromere_location + 1)
arm_end <- c(model_variables$cn_info$Centromere_location, model_variables$cn_info$Bin_count)
arm_s <- rep(1, length(arm_id))

model_variables <- BUILD_driver_library(
    model_variables = model_variables,
    table_arm_selection_rates = data.frame(Arm_ID = arm_id, Chromosome = arm_chromosome, Bin_start = arm_start, Bin_end = arm_end, s_rate = arm_s)
)

cell_count <- 20
CN_matrix <- BUILD_cn_normal_autosomes(model_variables$cn_info)
drivers <- list()
model_variables <- BUILD_initial_population(model_variables = model_variables, cell_count = cell_count, CN_matrix = CN_matrix, drivers = drivers)
library_name <- "ABC-BULK-ARM-CN"
model_variables <- CHECK_model_variables(model_variables)



list_parameters_library <- data.frame(matrix(ncol = 4, nrow = 0))
colnames(list_parameters_library) <- c("Variable", "Type", "Lower_bound", "Upper_bound")
list_parameters_library[nrow(list_parameters_library) + 1, ] <- c("prob_CN_chrom_arm_missegregation", "CNA_probability", 1e-5, 1e-4)
for (i in 1:length(model_variables$chromosome_arm_library$Arm_ID)) {
    list_parameters_library[nrow(list_parameters_library) + 1, ] <- c(
        model_variables$chromosome_arm_library$Arm_ID[i],
        "Arm_selection_rate",
        0.5, 1.5
    )
}



list_targets_library <- model_variables$chromosome_arm_library$Arm_ID



###################################    BUILD LIBRARY OF SIMULATED ARM-CN
start_time <- Sys.time()
library_bulk_arm_CN(
    library_name = library_name,
    model_variables = model_variables,
    list_parameters = list_parameters_library,
    list_targets = list_targets_library,
    ABC_simcount = 2000,
    n_samples = 100,
    R_libPaths = "/burg/iicd/users/knd2127/rpackages"
)
end_time <- Sys.time()
print(end_time - start_time)
########################################################################



#############################################    FIT FOR TCGA PAN-CANCER
copynumber_DATA <- read.csv(file = system.file("extdata", "Davoli_Charm_score.csv", package = "CancerSimulator"))
copynumber_DATA <- copynumber_DATA[, c(1, 2, 3, 14, 15)]
colnames(copynumber_DATA) <- c("Arm", "Del_freq_all", "Amp_freq_all", "Charm.TSG.OG.score", "Charm.TSG.OG.Ess.score")

list_targets <- copynumber_DATA$Arm
list_parameters <- list_parameters_library[c(
    which(list_parameters_library$Type == "CNA_probability"),
    which(list_parameters_library$Type == "Arm_selection_rate" & list_parameters_library$Variable %in% copynumber_DATA$Arm)
), ]

cat(paste0("\n\n\nFITTING FOR PAN-CANCER\n"))
fitting_bulk_arm_CN(
    library_name = library_name,
    model_name = "PAN-CANCER",
    model_variables = model_variables,
    copynumber_DATA = copynumber_DATA,
    list_parameters = list_parameters,
    list_parameters_library = list_parameters_library,
    list_targets = list_targets,
    list_targets_library = list_targets_library,
    library_shuffle = TRUE,
    type_sample_DATA = "average",
    type_cn_DATA = "arm",
    folder_workplace = library_name,
    R_libPaths = "/burg/iicd/users/knd2127/rpackages"
)
########################################################################



#######################################    FIT FOR PCAWG CANCER SUBTYPES
copynumber_PCAWG <- read_excel(system.file("pcawg_specimen_histology_August2016_v9.xlsx", package = "CancerSimulator"))

PCAWG_cancer_types <- unique(copynumber_PCAWG$histology_abbreviation)
PCAWG_cancer_types <- PCAWG_cancer_types[which(!is.na(PCAWG_cancer_types))]

PCAWG_wgd <- read.table(system.file("consensus.20170218.purity.ploidy.txt", package = "CancerSimulator"), header = TRUE)

PCAWG_cancer_type_sample_ids <- vector("list", length(PCAWG_cancer_types))
for (i in 1:length(PCAWG_cancer_types)) {
    sample_ids <- copynumber_PCAWG$tcga_sample_uuid[which(copynumber_PCAWG$histology_abbreviation == PCAWG_cancer_types[i])]
    for (j in 1:length(sample_ids)) {
        sample_id <- sample_ids[j]
        if (!file.exists(system.file(paste0("consensus.20170119.somatic.cna.annotated/", sample_id, ".consensus.20170119.somatic.cna.annotated.txt"), package = "CancerSimulator"))) next
        tmp <- read.table(system.file(paste0("consensus.20170119.somatic.cna.annotated/", sample_id, ".consensus.20170119.somatic.cna.annotated.txt"), package = "CancerSimulator"), header = TRUE)
        if (max(tmp$star, na.rm = TRUE) < 2) next
        PCAWG_cancer_type_sample_ids[[i]] <- c(PCAWG_cancer_type_sample_ids[[i]], sample_id)
    }
}

PCAWG_cancer_type_nonWGD_sample_ids <- vector("list", length(PCAWG_cancer_types))
for (i in 1:length(PCAWG_cancer_types)) {
    sample_ids <- PCAWG_cancer_type_sample_ids[[i]]
    sample_wgd <- PCAWG_wgd[which(PCAWG_wgd$samplename %in% sample_ids), ]
    sample_nonwgd_ids <- sample_wgd$samplename[which(sample_wgd$wgd_uncertain == FALSE & sample_wgd$wgd_status == "no_wgd")]
    PCAWG_cancer_type_nonWGD_sample_ids[[i]] <- sample_nonwgd_ids
}

sample_sizes <- sapply(PCAWG_cancer_type_nonWGD_sample_ids, length)
PCAWG_cancer_types <- PCAWG_cancer_types[-which(sample_sizes < 10)]
PCAWG_cancer_type_sample_ids <- PCAWG_cancer_type_sample_ids[-which(sample_sizes < 10)]
PCAWG_cancer_type_nonWGD_sample_ids <- PCAWG_cancer_type_nonWGD_sample_ids[-which(sample_sizes < 10)]

copynumber_DATA_cancer_types <- vector("list", length(PCAWG_cancer_types))
for (i in 1:length(PCAWG_cancer_types)) {
    cancer_type <- PCAWG_cancer_types[i]
    cancer_type_sample_ids <- PCAWG_cancer_type_nonWGD_sample_ids[[i]]
    copynumber_DATA_ls <- vector("list", length(cancer_type_sample_ids))
    for (j in 1:length(cancer_type_sample_ids)) {
        sample_id <- cancer_type_sample_ids[j]
        tmp <- read.table(system.file(paste0("consensus.20170119.somatic.cna.annotated/", sample_id, ".consensus.20170119.somatic.cna.annotated.txt"), package = "CancerSimulator"), header = TRUE)
        if (any(is.na(tmp$star))) tmp <- tmp[-which(is.na(tmp))]
        tmp <- tmp[which(tmp$star >= 2), ]
        tmp$donor_unique_id <- sample_id
        copynumber_DATA_ls[[j]] <- tmp
    }
    copynumber_DATA_cancer_types[[i]] <- rbindlist(copynumber_DATA_ls)
}

list_targets <- model_variables$chromosome_arm_library$Arm_ID
list_parameters <- list_parameters_library

for (i in 1:6) {
    cancer_type <- PCAWG_cancer_types[i]
    copynumber_DATA <- copynumber_DATA_cancer_types[[i]]
    cat(paste0("\n\n\nFITTING FOR ", cancer_type, "   [", i, "/", length(PCAWG_cancer_types), "]", "\n"))
    fitting_bulk_arm_CN(
        library_name = library_name,
        model_name = cancer_type,
        model_variables = model_variables,
        copynumber_DATA = copynumber_DATA,
        list_parameters = list_parameters,
        list_parameters_library = list_parameters_library,
        list_targets = list_targets,
        list_targets_library = list_targets_library,
        bound_freq = 0.1,
        library_shuffle = TRUE,
        folder_workplace = library_name,
        R_libPaths = "/burg/iicd/users/knd2127/rpackages"
    )
}
########################################################################



################################    STATISTICS FOR PCAWG CANCER SUBTYPES
df_WGD_FGA <- statistics_bulk_arm_WGD_status(
    plotname = "ALL_PCAWG",
    DATA_cancer_types = PCAWG_cancer_types,
    DATA_cancer_type_sample_ids = PCAWG_cancer_type_sample_ids,
    DATA_cancer_type_cn = copynumber_DATA_cancer_types,
    DATA_wgd = PCAWG_wgd,
    model_variables = model_variables
)
########################################################################
