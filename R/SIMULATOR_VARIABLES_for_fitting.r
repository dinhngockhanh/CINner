#================================SET ALL PARAMETERS REQUIRED FOR FITTING
SIMULATOR_VARIABLES_for_fitting <- function(model,parameter_set) {
#-------------------------------------Set up the default model variables
    SIMULATOR_VARIABLES_for_simulation(model)
#-------------------Set up the driver rate and CN event rates from input
    rate_driver                         <<- parameter_set$rate_driver

    prob_CN_whole_genome_duplication    <<- parameter_set$prob_CN_whole_genome_duplication
    prob_CN_missegregation              <<- parameter_set$prob_CN_missegregation
    prob_CN_chrom_arm_missegregation    <<- parameter_set$prob_CN_chrom_arm_missegregation
    prob_CN_focal_amplification         <<- parameter_set$prob_CN_focal_amplification
    prob_CN_focal_deletion              <<- parameter_set$prob_CN_focal_deletion
    prob_CN_cnloh_interstitial          <<- parameter_set$prob_CN_cnloh_interstitial
    prob_CN_cnloh_terminal              <<- parameter_set$prob_CN_cnloh_terminal

    prob_CN_focal_amplification_length  <<- parameter_set$prob_CN_focal_amplification_length
    prob_CN_focal_deletion_length       <<- parameter_set$prob_CN_focal_deletion_length
    prob_CN_cnloh_interstitial_length   <<- parameter_set$prob_CN_cnloh_interstitial_length
    prob_CN_cnloh_terminal_length       <<- parameter_set$prob_CN_cnloh_terminal_length
#---------------------------Set up the driver selection rates from input
    vec_selection_rates                 <<- parameter_set$vec_selection_rates
#   Count the number of TSGs and ONCOGENEs
    count_TSG                           <- sum(driver_library$Gene_role=='TSG')
    count_ONCOGENE                      <- sum(driver_library$Gene_role=='ONCOGENE')
#---Compute selection rates for TSGs
    list_TSG                            <- which(driver_library$Gene_role=='TSG')
    s_normalization                     <- 1
    for (driver=1:length(list_TSG)){
        row                             <- list_TSG[driver]
#       Get its selection strength
        driver_sel_rate                 <- vec_selection_rates[row]
#       Compute its selection rate for WT and MUT alleles
        driver_library$s_rate_WT[row]   <- 1/(1+driver_sel_rate)
        driver_library$s_rate_MUT[row]  <- 1
#       Update normalizer for selection rate
        s_normalization                 <- s_normalization*(1+driver_sel_rate)
    }
#---Compute selection rates for ONCOGENEs
    s_normalization                     <- s_normalization^(1/count_ONCOGENE)
    list_ONCOGENE                       <- which(driver_library$Gene_role=='ONCOGENE')
    for (driver=1:length(list_ONCOGENE)){
        row                             <- list_ONCOGENE[driver]
#       Get its selection strength
        driver_sel_rate                 <- vec_selection_rates[row]
#       Compute its selection rate for WT and MUT alleles
        driver_library$s_rate_WT[row]   <- s_normalization
        driver_library$s_rate_MUT[row]  <- s_normalization*(1+driver_sel_rate)
    }
}
