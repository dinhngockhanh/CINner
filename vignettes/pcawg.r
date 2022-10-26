.libPaths("/burg/iicd/users/knd2127/rpackages")
.libPaths()
library(ggplot2)
library(ggtree)
library(signals)
library(dendextend)
library(fishplot)
library(ctc)
library(adephylo)
library(data.table)
library(HMMcopy)
library(plyr)
library(matrixStats)
library(parallel)
library(RColorBrewer)
library(ape)

library(gridExtra)

library(readxl)
library(phyloTop)
library(pbapply)

library(abc)
library(abcrf)

library(reticulate)
library(fitdistrplus)
library(ggrepel)
#-----------------------------------------------------------------------
tmp <- getwd()
setwd("/burg/home/knd2127/R/")
files_sources <- list.files(pattern = "*.r$")
sapply(files_sources, source)
setwd(tmp)

# devtools::install_github("dinhngockhanh/CancerSimulator")
# library(CancerSimulator)
#-----------------------------------------------------------------------



copynumber_PCAWG <- read.csv(file = system.file("extdata", "PCAWG_OV-AU.csv", package = "CancerSimulator"))



cell_lifespan <- 30
T_0 <- list(0, "year")
T_end <- list(80, "year")
Table_sample <- data.frame(Sample_ID = c("SA01"), Cell_count = c(Inf), Age_sample = c(80))
T_tau_step <- cell_lifespan / 2
CN_bin_length <- 500000
selection_model <- "chrom-arm-and-driver-gene-selection"
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
#-----------------------------------------------------------------------
rate_passenger <- 1e-11
bound_driver <- 3
bound_average_ploidy <- 4.5
bound_homozygosity <- 10

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

driver_gene_list <- read_excel("HGSOC_driver_genes.xlsx")
gene_id <- driver_gene_list$Gene_ID
gene_role <- driver_gene_list$Gene_role
gene_chromosome <- driver_gene_list$Chromosome
gene_bin <- round(driver_gene_list$Start / CN_bin_length)
gene_s <- rep(1, length(gene_id))

table_arm_selection_rates <- data.frame(Arm_ID = arm_id, Chromosome = arm_chromosome, Bin_start = arm_start, Bin_end = arm_end, s_rate = arm_s)
table_gene_selection_rates <- data.frame(Gene_ID = gene_id, Gene_role = gene_role, s_rate = gene_s, Chromosome = gene_chromosome, Bin = gene_bin)

model_variables <- BUILD_driver_library(
    model_variables = model_variables,
    table_arm_selection_rates = table_arm_selection_rates,
    table_gene_selection_rates = table_gene_selection_rates
)
cell_count <- 20
CN_matrix <- BUILD_cn_normal_XX(model_variables$cn_info)
drivers <- list()
model_variables <- BUILD_initial_population(model_variables = model_variables, cell_count = cell_count, CN_matrix = CN_matrix, drivers = drivers)
model_name <- "HGSOC"
model_variables <- CHECK_model_variables(model_variables)



list_parameters <- data.frame(matrix(ncol = 3, nrow = 0))
colnames(list_parameters) <- c("Variable", "Lower_bound", "Upper_bound")
list_parameters[nrow(list_parameters) + 1, ] <- c("prob_CN_chrom_arm_missegregation", 1e-5, 1e-4)
for (i in 1:length(model_variables$chromosome_arm_library$Arm_ID)) {
    list_parameters[nrow(list_parameters) + 1, ] <- c(
        model_variables$chromosome_arm_library$Arm_ID[i],
        0.5, 1.5
    )
}



list_targets <- model_variables$chromosome_arm_library$Arm_ID



fitting_arm_PCAWG(
    model_name,
    model_variables,
    copynumber_PCAWG,
    list_parameters,
    list_targets,
    ABC_method = "rf",
    ABC_simcount = 10000
)
