# ======================================================================
# ======================================================================
# ============================================================ TERREMOTO
.libPaths("/moto/iicd/users/knd2127/rpackages")
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
library(hrbrthemes)
library(viridis)
#-----------------------------------------------------------------------
setwd("/moto/home/knd2127/CancerSimulator/R/")
files_sources <- list.files(pattern = "*.r$")
sapply(files_sources, source)
#-----------------------------------------------------------------------
cell_lifespan <- 30
T_0 <- list(0, "year")
T_end <- list(80, "year")
Table_sample <- data.frame(Sample_ID = "SA01", Cell_count = 300, Age_sample = 80)
T_tau_step <- cell_lifespan / 2
CN_bin_length <- 500000
prob_CN_whole_genome_duplication <- 3e-4
prob_CN_missegregation <- 7e-4
prob_CN_chrom_arm_missegregation <- 2e-4
prob_CN_focal_amplification <- 7e-4
prob_CN_focal_deletion <- 2e-4
prob_CN_cnloh_interstitial <- 4e-4
prob_CN_cnloh_terminal <- 6e-4
prob_CN_focal_amplification_length <- 0.5
prob_CN_focal_deletion_length <- 0.5
prob_CN_cnloh_interstitial_length <- 0.5
prob_CN_cnloh_terminal_length <- 0.5
rate_driver <- 8e-16
rate_passenger <- 1e-11
bound_driver <- 3
bound_average_ploidy <- 4.5
bound_homozygosity <- 500
vec_time <- T_0[[1]]:T_end[[1]]
L <- 10000
t_0 <- 20
k <- 0.3
vec_cell_count <- L / (1 + exp(-k * (vec_time - t_0)))
table_population_dynamics <- cbind(vec_time, vec_cell_count)
gc <- read.csv(file = "../data/gc_map_500kb.csv")
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
    prob_CN_focal_amplification_length = prob_CN_focal_amplification_length,
    prob_CN_focal_deletion_length = prob_CN_focal_deletion_length,
    prob_CN_cnloh_interstitial_length = prob_CN_cnloh_interstitial_length,
    prob_CN_cnloh_terminal_length = prob_CN_cnloh_terminal_length,
    rate_driver = rate_driver,
    rate_passenger = rate_passenger,
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
vec_driver_genes <- c(
    "TP53", "BRCA1", "BRCA2", "MSH2", "MSH6", "ARID1A",
    "STK11", "BRAF", "MLH1", "PIK3R1", "ERBB2", "CTNNB1", "AKT1", "PPP2R1A"
)
vec_driver_role <- c(
    "TSG", "TSG", "TSG", "TSG", "TSG", "TSG",
    "TSG", "ONCOGENE", "TSG", "TSG", "ONCOGENE", "ONCOGENE", "ONCOGENE", "TSG"
)
vec_driver_s <- rep(0.007, length(vec_driver_genes))
vec_driver_s[1] <- 0.01
vec_driver_s[2] <- 0.009
vec_driver_s[3] <- 0.009
model_variables <- BUILD_driver_library(model_variables = model_variables, vec_driver_genes = vec_driver_genes, vec_driver_role = vec_driver_role, vec_driver_s = vec_driver_s)
cell_count <- 20
CN_matrix <- BUILD_cn_normal_XX(model_variables$cn_info)
drivers <- list()
model_variables <- BUILD_initial_population(model_variables = model_variables, cell_count = cell_count, CN_matrix = CN_matrix, drivers = drivers)
model_name <- "TEST"
SAVE_model_variables(model_name = model_name, model_variables = model_variables)
n_simulations <- 1000
stage_final <- 3
n_clones_min <- 1
n_clones_max <- Inf

start_time <- Sys.time()
simulator_full_program(
    model = model_name,
    n_simulations = n_simulations,
    stage_final = stage_final,
    n_clones_min = n_clones_min,
    n_clones_max = n_clones_max,
    save_simulation = TRUE,
    internal_nodes_cn_info = FALSE,
    save_newick_tree = TRUE,
    save_cn_profile = FALSE,
    model_readcount = FALSE,
    format_cn_profile = "both",
    apply_HMM = FALSE,
    apply_UMAP_on_HMM = FALSE,
    report_progress = FALSE,
    seed = 10,
    compute_parallel = TRUE
)
end_time <- Sys.time()
print(end_time - start_time)
plot_all(model = model_name, n_simulations = n_simulations)
