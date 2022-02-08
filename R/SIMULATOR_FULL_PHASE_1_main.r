#==============================================PHASE 1: CLONAL EVOLUTION
#---List of global variables:
    # print('CHECKING:')
    # print(genotype_list_ploidy_chrom)

    # print(genotype_list_ploidy_block)

    # print(genotype_list_driver_count)
    # print(genotype_list_driver_map)

    # print(genotype_list_selection_rate)
    # print(genotype_list_DNA_length)
    # print(genotype_list_prob_new_drivers)

    # print(N_clones)
    # print(clonal_ID_current)
    # print(clonal_population_current)
    # print(clonal_population_next)
    # print(evolution_origin)
    # print(evolution_genotype_changes)
SIMULATOR_FULL_PHASE_1_main <- function() {
#-----------------------------------------Set up the initial CN genotype
    cell_vec_ploidy_chrom                       <- rep(1,N_chromosomes)
    cell_vec_ploidy_chrom[1:N_chromosomes]      <- 2

    cell_mat_ploidy_block                       <- list()
    for(chrom in 1:N_chromosomes){
        ploidy                                  <- cell_vec_ploidy_chrom[chrom]
        no_blocks                               <- vec_CN_block_no[chrom]
        cell_mat_ploidy_block[[chrom]]          <- list()
        for(strand in 1:ploidy){
            cell_mat_ploidy_block[[chrom]][[strand]]    <- rep(1,no_blocks)
        }
    }

    cell_vec_ploidy_allele                      <- list()
    for(chrom in 1:N_chromosomes){
        ploidy                                  <- cell_vec_ploidy_chrom[chrom]
        no_blocks                               <- vec_CN_block_no[chrom]
        cell_vec_ploidy_allele[[chrom]]         <- list()
        for(strand in 1:ploidy){
            cell_vec_ploidy_allele[[chrom]][[strand]]   <- rep(strand,no_blocks)
        }
    }

    genotype_list_ploidy_chrom                  <<- list()
    genotype_list_ploidy_chrom[[1]]             <<- cell_vec_ploidy_chrom

    genotype_list_ploidy_allele                 <<- list()
    genotype_list_ploidy_allele[[1]]            <<- cell_vec_ploidy_allele

    genotype_list_ploidy_block                  <<- list()
    genotype_list_ploidy_block[[1]]             <<- cell_mat_ploidy_block

    assign('genotype_list_ploidy_chrom',genotype_list_ploidy_chrom,envir=.GlobalEnv)
    assign('genotype_list_ploidy_allele',genotype_list_ploidy_chrom,envir=.GlobalEnv)
    assign('genotype_list_ploidy_block',genotype_list_ploidy_block,envir=.GlobalEnv)
#-------------------------------------Set up the initial driver genotype
    cell_mat_drivers                            <- c(0)
    genotype_list_driver_count                  <<- matrix(0,1,1)
    genotype_list_driver_map                    <<- list()
    genotype_list_driver_map[[1]]               <<- cell_mat_drivers
    genotype_list_selection_rate                <<- matrix(0,1,1)
    genotype_list_selection_rate[1]             <<- SIMULATOR_FULL_PHASE_1_selection_rate(0,cell_mat_drivers,cell_vec_ploidy_chrom,cell_mat_ploidy_block)
    assign('genotype_list_driver_count',genotype_list_driver_count,envir=.GlobalEnv)
    assign('genotype_list_driver_map',genotype_list_driver_map,envir=.GlobalEnv)
    assign('genotype_list_selection_rate',genotype_list_selection_rate,envir=.GlobalEnv)
#--------------------------Set up the DNA length of the initial genotype
    DNA_length                                  <- 0
    for (chrom in 1:N_chromosomes) {
        for (strand in 1:cell_vec_ploidy_chrom[chrom]) {
            DNA_length                          <- DNA_length+sum(cell_mat_ploidy_block[[chrom]][[strand]])
        }
    }
    DNA_length                                  <- size_CN_block_DNA*DNA_length
    genotype_list_DNA_length                    <<- list(DNA_length)
    genotype_list_prob_new_drivers              <<- 1-dpois(x=0,lambda=rate_driver*DNA_length)
    assign('genotype_list_DNA_length',genotype_list_DNA_length,envir=.GlobalEnv)
    assign('genotype_list_prob_new_drivers',genotype_list_prob_new_drivers,envir=.GlobalEnv)
#-------------------------------------Set up the clonal evolution record
    N_clones                                    <<- 1
    assign('N_clones',N_clones,envir=.GlobalEnv)
#   Set up the record for the clonal evolution
    clonal_ID_current                           <<- c(N_clones)
    clonal_population_current                   <<- c(1)
    clonal_population_next                      <<- clonal_population_current
    assign('clonal_ID_current',clonal_ID_current,envir=.GlobalEnv)
    assign('clonal_population_current',clonal_population_current,envir=.GlobalEnv)
    assign('clonal_population_next',clonal_population_next,envir=.GlobalEnv)

    evolution_origin                            <<- c(0)
    evolution_genotype_changes                  <<- list(c())
    evolution_traj_time                         <- c(T_start_time)
    evolution_traj_clonal_ID                    <- list(clonal_ID_current)
    evolution_traj_population                   <- list(clonal_population_current)
    evolution_traj_divisions                    <- list()
    assign('evolution_origin',evolution_origin,envir=.GlobalEnv)
    assign('evolution_genotype_changes',evolution_genotype_changes,envir=.GlobalEnv)
#--------------------------------------Set up counts for the simulations
#   Current time
    T_current                       <- T_start_time
    T_goal                          <- T_end_time
#   Count of cells still alive at current time
    N_cells_current                 <- 1
    N_cells_goal                    <- Population_end
#   Count of events
    N_events_current                <- 0
    N_events_goal                   <- Max_events
#-----------------------------------------------------Simulation process
    while ((T_current<T_goal)&&(N_cells_current<N_cells_goal)&&(N_events_current<N_events_goal)) {
#       Find the Poisson propensities of event count for all clones
        rate_base_lifetime          <- func_event_rate(T_current)
        all_propensity              <- T_tau_step*rate_base_lifetime*clonal_population_current
#       Find the probability of division for all clones
        clonal_portion              <- genotype_list_selection_rate[clonal_ID_current]
        all_prob_division           <- func_expected_population(T_current)/(func_expected_population(T_current)+N_cells_current) * sum(clonal_population_current)*clonal_portion/sum(clonal_portion*clonal_population_current)
#       Find next time step and initiate next clonal population vector
        T_next                      <- T_current+T_tau_step
        clonal_population_next      <<- clonal_population_current
        if ((N_cells_current<=0) || (is.null(N_cells_current))) {
            flag_success            <<- 0 # Let flag_success as global variable
            break
            }
        else { if (T_next>=T_end_time) {
            flag_success            <<- 1 # Let flag_success as global variable
            T_current               <- T_end_time
            break
            }
        }
#       Report on progress
        # if(floor(T_current/365)<floor(T_next/365)){
        #     cat('Event = ',N_events_current,'/',N_events_goal,';   Time = ',round(T_current/365),'/',T_goal/365,';   Population = ',N_cells_current,'/',N_cells_goal,';   Number of clones = ',N_clones,'\n')
        # }
#       Initialize the matrix of divisions for this step
        mat_divisions               <- c()
#       Find all existing clones
        all_existing_clones         <- clonal_ID_current
#       For every existing clones...
        if (length(all_existing_clones) >=1) {for (i in 1:length(all_existing_clones)) {
#           Find clone ID
            position_to_react       <- i
            clone_to_react          <- all_existing_clones[i]
            genotype_to_react       <- clone_to_react
#           Find current clonal population
            clone_population        <- clonal_population_current[i]
#           Find probability of division
            prob_division           <- all_prob_division[i]
#           Find probability of new genotype
            DNA_length              <- genotype_list_DNA_length[[clone_to_react]]
            prob_new_drivers        <- genotype_list_prob_new_drivers[clone_to_react]
            prob_new_genotype       <- 1-(1-prob_CN_whole_genome_duplication)*(1-prob_CN_missegregation)*(1-prob_CN_chrom_arm_missegregation)*(1-prob_CN_focal_amplification)*(1-prob_CN_focal_deletion)*(1-prob_new_drivers)
#           Find number of events
            prop                    <- all_propensity[i]
            count_new_events        <- Inf
            while ((count_new_events>clone_population)||(count_new_events<0)) {
                if (prop>1000) {
                    count_new_events <- round(rnorm(n=1,mean=prop,sd=sqrt(prop)))
                    }
                else {
                count_new_events    <- rpois(n=1,lambda=prop)
                }
            }
            N_events_current        <- N_events_current+count_new_events
            count_event_types       <- rmultinom(n=1,size=count_new_events,prob=c((1-prob_division), prob_division*(1-prob_new_genotype), prob_division*prob_new_genotype))
            count_deaths            <- count_event_types[1]
            count_div_old           <- count_event_types[2]
            count_div_new_tmp       <- count_event_types[3]
            count_div_new           <- count_div_new_tmp
#-----------Perform death events
            clonal_population_next[i]               <<- clonal_population_next[i]-count_deaths
#-----------Perform division events with no new genotype
            clonal_population_next[i]               <<- clonal_population_next[i]+count_div_old
#-----------Perform division events with new genotype
            if (count_div_new_tmp>=1){ for (j in 1:count_div_new_tmp) {
#               Find what events lead to the new genotype
                flag_whole_genome_duplication       <- 0
                flag_missegregation                 <- 0
                flag_chrom_arm_missegregation       <- 0
                flag_amplification                  <- 0
                flag_drivers                        <- 0
                while (max(c(flag_whole_genome_duplication, flag_missegregation, flag_chrom_arm_missegregation, flag_amplification, flag_drivers))==0) {
                    flag_drivers                    <- as.numeric(runif(1)<(prob_new_drivers/prob_new_genotype))
                    flag_whole_genome_duplication   <- as.numeric(runif(1)<(prob_CN_whole_genome_duplication/prob_new_genotype))
                    flag_missegregation             <- as.numeric(runif(1)<(prob_CN_missegregation/prob_new_genotype))
                    flag_chrom_arm_missegregation   <- as.numeric(runif(1)<(prob_CN_chrom_arm_missegregation/prob_new_genotype))
                    flag_amplification              <- as.numeric(runif(1)<(prob_CN_focal_amplification/prob_new_genotype))
                    flag_deletion                   <- as.numeric(runif(1)<(prob_CN_focal_deletion/prob_new_genotype))
                }
#               Initiate the two new genotypes
                output              <- SIMULATOR_FULL_PHASE_1_genotype_initiation(genotype_to_react)
                genotype_daughter_1 <- output[[1]]
                genotype_daughter_2 <- output[[2]]
                position_daughter_1 <- output[[3]]
                position_daughter_2 <- output[[4]]
#               Simulate new driver event
                if (flag_drivers==1) {
                    SIMULATOR_FULL_PHASE_1_drivers(genotype_to_react,genotype_daughter_1,genotype_daughter_2)
                }
#               Simulate whole genome duplication event
                if (flag_whole_genome_duplication==1) {
                    SIMULATOR_FULL_PHASE_1_CN_whole_genome_duplication(genotype_to_react,genotype_daughter_1,genotype_daughter_2)
                }
#               Simulate missegregation event
                if (flag_missegregation==1) {
                    SIMULATOR_FULL_PHASE_1_CN_missegregation(genotype_to_react,genotype_daughter_1,genotype_daughter_2)
                }
#               Simulate chromosome-arm missegregation event
                if (flag_chrom_arm_missegregation==1) {
                    SIMULATOR_FULL_PHASE_1_CN_chrom_arm_missegregation(genotype_to_react,genotype_daughter_1,genotype_daughter_2)
                }
#               Simulate focal amplification event
                if (flag_amplification==1) {
                    if (sample.int(2,size=1)==1) {
                        SIMULATOR_FULL_PHASE_1_CN_focal_amplification(genotype_to_react,genotype_daughter_1)
                        }
                    else {
                        SIMULATOR_FULL_PHASE_1_CN_focal_amplification(genotype_to_react,genotype_daughter_2)
                    }
                }
#               Simulate focal deletion event
                if (flag_deletion==1) {
                    if (sample.int(2,size=1)==1) {
                        SIMULATOR_FULL_PHASE_1_CN_focal_deletion(genotype_to_react,genotype_daughter_1)
                        }
                    else {
                        SIMULATOR_FULL_PHASE_1_CN_focal_deletion(genotype_to_react,genotype_daughter_2)
                    }
                }
#               Update DNA length and selection rates of daughter cells
                SIMULATOR_FULL_PHASE_1_genotype_update(genotype_daughter_1,genotype_daughter_2)
#               Check if either daughter cell did not create a new clone
                output                                      <- SIMULATOR_FULL_PHASE_1_genotype_cleaning(genotype_to_react,genotype_daughter_1,genotype_daughter_2,position_to_react,position_daughter_1,position_daughter_2)
                genotype_to_react                           <- output[[1]]
                genotype_daughter_1                         <- output[[2]]
                genotype_daughter_2                         <- output[[3]]
                position_to_react                           <- output[[4]]
                position_daughter_1                         <- output[[5]]
                position_daughter_2                         <- output[[6]]
#               Adjust the event count accordingly, and add event to
#               matrix of divisions
                if ((genotype_daughter_1==genotype_to_react)&&(genotype_daughter_2==genotype_to_react)) {
                   count_div_new                            <- count_div_new-1
                   count_div_old                            <- count_div_old+1
                   }
                else {
                    mat_divisions                           <- rbind(mat_divisions, c(1, genotype_to_react, genotype_daughter_1, genotype_daughter_2))
                }
#               Update the clonal population according to what happens
                clonal_population_next[position_to_react]   <<- clonal_population_next[position_to_react]-1
                clonal_population_next[position_daughter_1] <<- clonal_population_next[position_daughter_1]+1
                clonal_population_next[position_daughter_2] <<- clonal_population_next[position_daughter_2]+1
            }}
#           Add the divisions with old genotype to the matrix of divisions
            if (count_div_old>0) {
                mat_divisions                               <- rbind(mat_divisions, c(count_div_old, genotype_to_react, genotype_to_react, genotype_to_react))
            }
        }}
#       Clean clonal populations
        SIMULATOR_FULL_PHASE_1_clonal_population_cleaning()
#       Update clonal populations
        clonal_population_current           <<- clonal_population_next
#       Update time
        T_current                           <- T_next
#       Update count of cells
        N_cells_current                     <- sum(clonal_population_current)
#       Update record of clonal evolution over time
        evolution_traj_time                                              <- c(evolution_traj_time, T_current)
        evolution_traj_clonal_ID[[length(evolution_traj_clonal_ID)+1]]   <- clonal_ID_current
        evolution_traj_population[[length(evolution_traj_population)+1]] <- clonal_population_current
        evolution_traj_divisions[[length(evolution_traj_divisions)+1]]   <- mat_divisions
    }
#---------------------------------Output package of data from simulation
    if (is.null(N_cells_current)) {
        flag_success                        <<- 0 # Let flag_success as global variable
    }
#   Let package_clonal_evolution as a vector of 15 lists or by list() and then define each element
    package_clonal_evolution                <- list()
    package_clonal_evolution[[1]]           <- T_current
    package_clonal_evolution[[2]]           <- N_cells_current
    package_clonal_evolution[[3]]           <- N_events_current
    package_clonal_evolution[[4]]           <- N_clones
    package_clonal_evolution[[5]]           <- genotype_list_ploidy_chrom
    package_clonal_evolution[[6]]           <- genotype_list_ploidy_block



    package_clonal_evolution[[7]]           <- genotype_list_ploidy_allele



    package_clonal_evolution[[8]]           <- genotype_list_driver_count
    package_clonal_evolution[[9]]           <- genotype_list_driver_map
    package_clonal_evolution[[10]]          <- genotype_list_selection_rate
    package_clonal_evolution[[11]]          <- evolution_origin
    package_clonal_evolution[[12]]          <- evolution_genotype_changes
    package_clonal_evolution[[13]]          <- evolution_traj_time
    package_clonal_evolution[[14]]          <- evolution_traj_divisions
    package_clonal_evolution[[15]]          <- evolution_traj_clonal_ID
    package_clonal_evolution[[16]]          <- evolution_traj_population


final_clonal_ID     <- evolution_traj_clonal_ID[[length(evolution_traj_clonal_ID)]]
for(i in 1:length(final_clonal_ID)){
    clonal_ID       <- final_clonal_ID[i]
    print('-----------------------------------------------------------')
    print(genotype_list_ploidy_chrom[[clonal_ID]])
    # print(genotype_list_ploidy_block[[clonal_ID]])
    # print(genotype_list_ploidy_chrom[[clonal_ID]])
}


    output                                  <- list()
    output[[1]]                             <- flag_success
    output[[2]]                             <- package_clonal_evolution
    return(output)
}
