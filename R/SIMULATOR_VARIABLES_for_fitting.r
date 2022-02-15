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
    
print(vec_selection_rates)



}
