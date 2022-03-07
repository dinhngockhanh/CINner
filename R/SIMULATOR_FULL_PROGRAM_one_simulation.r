#==================================================CREATE ONE SIMULATION
SIMULATOR_FULL_PROGRAM_one_simulation <- function(model,stage_final) {
    SIMULATOR_VARIABLES_for_simulation(model)
#------------------------------------------Simulate the clonal evolution
    flag_success                                    <- 0
    while (flag_success==0){
print('')
print('CLONAL EVOLUTION...')
start.time  <- Sys.time()
        output                                      <- SIMULATOR_FULL_PHASE_1_main()
        flag_success                                <- output[[1]]
        package_clonal_evolution                    <- output[[2]]
end.time    <- Sys.time()
time.taken  <- end.time - start.time
print(time.taken)
    }
    if(stage_final>=1){
        package_output                              <- list()
        package_output[[1]]                         <- package_clonal_evolution
    }
#----------------------------------------------------Simulate the sample
    if(stage_final>=2){
print('')
print('SAMPLING...')
start.time  <- Sys.time()
        package_sample                              <- SIMULATOR_FULL_PHASE_2_main(package_clonal_evolution)
        package_output[[2]]                         <- package_sample
end.time    <- Sys.time()
time.taken  <- end.time - start.time
print(time.taken)
    }
#   Save the sample's cell CN profiles
    sample_genotype_profiles                        <- package_sample[[1]]
    filename                                        <- paste(model,'-output-cn-profiles','.rda',sep='')
    # write.csv(sample_genotype_profiles,filename,row.names=FALSE)
    save(sample_genotype_profiles,file=filename)
#-----------------------------------Simulate the phylogeny of the sample
    if(stage_final>=3){
print('')
print('SAMPLE PHYLOGENY...')
start.time  <- Sys.time()
        package_sample_phylogeny                    <- SIMULATOR_FULL_PHASE_3_main(package_clonal_evolution,package_sample)
        package_output[[3]]                         <- package_sample_phylogeny
end.time    <- Sys.time()
time.taken  <- end.time - start.time
print(time.taken)
    }
#   Save the sample's cell phylogeny
    phylogeny_clustering_truth                      <- package_sample_phylogeny[[1]]
    filename                                        <- paste(model,'-output-clustering','.rda',sep='')
    save(phylogeny_clustering_truth,file=filename)



print('')
print('DONE WITH SIMULATION...')

#------------------------------------------Output the simulation package
    return(package_output)
}
