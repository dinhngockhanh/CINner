#============================CREATE ONE SIMULATION WITH INPUT PARAMETERS
ABC_one_simulation <- function(model,parameter_set) {
#---------------------------------------------Set up the model variables
    SIMULATOR_VARIABLES_for_fitting(model,parameter_set)
#------------------------------------------Simulate the clonal evolution
    flag_success                                    <- 0
    while (flag_success==0) {
print('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   CLONAL EVOLUTION')
start.time  <- Sys.time()
        output                                      <- SIMULATOR_FULL_PHASE_1_main()
        flag_success                                <- output[[1]]
        package_clonal_evolution                    <- output[[2]]
end.time    <- Sys.time()
time.taken  <- end.time - start.time
print(time.taken)
    }
    package_output                                  <- list()
    package_output[[1]]                             <- package_clonal_evolution
print('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   SAMPLING')
start.time  <- Sys.time()
    package_sample                                  <- SIMULATOR_FULL_PHASE_2_main(package_clonal_evolution)
    package_output[[2]]                             <- package_sample
end.time    <- Sys.time()
time.taken  <- end.time - start.time
print(time.taken)
print('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   SAMPLE PHYLOGENY')
start.time  <- Sys.time()
    package_sample_phylogeny                        <- SIMULATOR_FULL_PHASE_3_main(package_clonal_evolution,package_sample)
    package_output[[3]]                             <- package_sample_phylogeny
end.time    <- Sys.time()
time.taken  <- end.time - start.time
print(time.taken)
print('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   DONE WITH SIMULATION')
#------------------------------------------Output the simulation package
    return(package_output)
}
