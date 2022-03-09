#==============================================PHASE 2: SAMPLE PHYLOGENY
SIMULATOR_FULL_PHASE_2_main <- function(package_clonal_evolution) {
#---------------------------------------------Input the clonal evolution
    T_current                       <- package_clonal_evolution[[1]]
    N_cells_current                 <- package_clonal_evolution[[2]]
    N_clones                        <- package_clonal_evolution[[4]]
    genotype_list_ploidy_chrom      <- package_clonal_evolution[[5]]
    genotype_list_ploidy_block      <- package_clonal_evolution[[6]]
    genotype_list_ploidy_allele     <- package_clonal_evolution[[7]]

    genotype_list_driver_map        <- package_clonal_evolution[[9]]

    evolution_traj_time             <- package_clonal_evolution[[13]]
    evolution_traj_divisions        <- package_clonal_evolution[[14]]
    evolution_traj_clonal_ID        <- package_clonal_evolution[[15]]
    evolution_traj_population       <- package_clonal_evolution[[16]]

    for (row in 1:nrow(Table_sampling)){
        if (Table_sampling$Cell_count[row]==Inf){
            loc                     <- which.min(abs(evolution_traj_time-Table_sampling$T_sample[row]))
            N_cells_total           <- sum(evolution_traj_population[[loc]])
            Table_sampling$Cell_count[row]  <- N_cells_total
        }
    }

#-------------------------------Find a random sample of final population
    all_sample_genotype             <- c()
    all_sample_ID                   <- c()
    all_sample_sampled_time         <- c()
    for (sample in 1:nrow(Table_sampling)){
        N_sample                    <- Table_sampling$Cell_count[sample]
        ID_sample                   <- Table_sampling$Sample_ID[sample]
        T_sample                    <- Table_sampling$T_sample[sample]

        loc                         <- which.min(abs(evolution_traj_time-T_sample))
        vec_clonal_ID               <- evolution_traj_clonal_ID[[loc]]
        vec_clonal_population       <- evolution_traj_population[[loc]]
        vec_population              <- c()
        for (i in 1:length(vec_clonal_ID)) {
            clone                   <- vec_clonal_ID[i]
            clonal_population       <- vec_clonal_population[i]
            vec_population          <- c(vec_population,rep(clone,1,clonal_population))
        }

        sample_genotype             <- sample(x=vec_population,size=N_sample,replace=FALSE)
        all_sample_genotype         <- c(all_sample_genotype,sample_genotype)
        all_sample_ID               <- c(all_sample_ID,rep(ID_sample,N_sample))
        all_sample_sampled_time     <- c(all_sample_sampled_time,rep(T_sample,N_sample))

        print(paste('DETECTED ',length(unique(sample_genotype)),' CLONES IN SAMPLE ',ID_sample,sep=''))
    }
    print(paste('DETECTED ',length(unique(all_sample_genotype)),' CLONES IN ALL SAMPLES',sep=''))
#---------------------------------Create CN object for the sampled cells
#---Find the CN profiles for each clone found in the sample
    sample_genotype_unique          <- unique(all_sample_genotype)
    sample_genotype_unique_profile  <- list()
    for(i_clone in 1:length(sample_genotype_unique)){
#       Extract CN information for the clone from clonal evolution data
        clone_ID                    <- sample_genotype_unique[i_clone]
        ploidy_chrom                <- genotype_list_ploidy_chrom[[clone_ID]]
        ploidy_block                <- genotype_list_ploidy_block[[clone_ID]]
        ploidy_allele               <- genotype_list_ploidy_allele[[clone_ID]]
#       Build the CN profile in SIGNALS style for the clone
        vec_clone_chr               <- c()
        vec_clone_start             <- c()
        vec_clone_end               <- c()
        vec_clone_copy              <- c()
        vec_clone_state             <- c()
        vec_clone_Min               <- c()
        vec_clone_Maj               <- c()
        for (chrom in 1:N_chromosomes){
            chrom_block_count       <- vec_CN_block_no[chrom]
            chrom_ploidy            <- ploidy_chrom[chrom]
#           Find location information of each chromosome block
            vec_chr                 <- rep(as.character(chrom),1,chrom_block_count)
            vec_start               <- seq(0,size_CN_block_DNA*(chrom_block_count-1),by=size_CN_block_DNA)+1
            vec_end                 <- seq(size_CN_block_DNA,size_CN_block_DNA*chrom_block_count,by=size_CN_block_DNA)
#           Find CN counts for each allele of each chromosome block
            vec_Allele_1            <- rep(0,1,chrom_block_count)
            vec_Allele_2            <- rep(0,1,chrom_block_count)
            for (strand in 1:chrom_ploidy){
                mat_allele                  <- ploidy_allele[[chrom]][[strand]]
                for (CN_row in 1:nrow(mat_allele)){
                    list_1                  <- which(mat_allele[CN_row,]==1)
                    vec_Allele_1[list_1]    <- vec_Allele_1[list_1]+1
                    list_2                  <- which(mat_allele[CN_row,]==2)
                    vec_Allele_2[list_2]    <- vec_Allele_2[list_2]+1
                }
            }
#           Find Major/Minor CN counts of each chromosome block
            if (mean(vec_Allele_1)<=mean(vec_Allele_2)){
                vec_Min             <- vec_Allele_1
                vec_Maj             <- vec_Allele_2
            }
            else{
                vec_Min             <- vec_Allele_2
                vec_Maj             <- vec_Allele_1
            }
#           Find total CN count of each chromosome block
            vec_copy                <- vec_Min+vec_Maj
            vec_state               <- vec_copy
#           Update the CN information of the clone
            vec_clone_chr           <- c(vec_clone_chr,vec_chr)
            vec_clone_start         <- c(vec_clone_start,vec_start)
            vec_clone_end           <- c(vec_clone_end,vec_end)
            vec_clone_copy          <- c(vec_clone_copy,vec_copy)
            vec_clone_state         <- c(vec_clone_state,vec_state)
            vec_clone_Min           <- c(vec_clone_Min,vec_Min)
            vec_clone_Maj           <- c(vec_clone_Maj,vec_Maj)
        }
#       Store the CN profile for the clone
        genotype_unique_profile                     <- data.frame(vec_clone_chr,vec_clone_start,vec_clone_end,vec_clone_copy,vec_clone_state,vec_clone_Min,vec_clone_Maj)
        names(genotype_unique_profile)              <- c("chr","start","end","copy","state","Min","Maj")
        sample_genotype_unique_profile[[i_clone]]   <- genotype_unique_profile
    }
#---Find the CN profiles for each cell in the sample
    sample_cell_ID                  <- c()
    sample_clone_ID                 <- all_sample_genotype
    sample_time                     <- all_sample_sampled_time

    for (i_cell in 1:length(all_sample_genotype)){
        sample_ID                   <- all_sample_ID[i_cell]
        clone_ID                    <- sample_clone_ID[i_cell]
        i_clone                     <- which(sample_genotype_unique==clone_ID)[1]
#       Find the CN profile for this cell
        cell_genotype_profile       <- sample_genotype_unique_profile[[i_clone]]
#       Add column for cell ID
        cell_ID                     <- paste(sample_ID,'-Library-',as.character(i_cell),'-',as.character(i_cell),sep='')
        # cell_ID                     <- paste('Sample-Library-',as.character(i_cell),'-',as.character(i_cell),sep='')
        sample_cell_ID[i_cell]      <- cell_ID
        cell_genotype_profile       <- cbind(cell_genotype_profile,rep(cell_ID,nrow(cell_genotype_profile)))
        names(cell_genotype_profile)<- c("chr","start","end","copy","state","Min","Maj","cell_id")
#       Update table of CN profiles for all cells in the sample
        if (i_cell==1){
            sample_genotype_profile <- cell_genotype_profile
        }
        else{
            sample_genotype_profile <- rbind(sample_genotype_profile,cell_genotype_profile)
        }
    }
#--------------------------Give each clone a character index (A,B,C,...)
    sample_clone_ID_numeric                         <- sample_clone_ID
    sample_clone_ID_unique_numeric                  <- unique(sample_clone_ID_numeric)
    sample_clone_ID_unique_letters                  <- LETTERS[as.numeric(1:length(sample_clone_ID_unique_numeric))]
    table_clone_ID_vs_letters                       <- data.frame(sample_clone_ID_unique_numeric,sample_clone_ID_unique_letters)
    colnames(table_clone_ID_vs_letters)             <- c('Clone_ID_number','Clone_ID_letter')
    sample_clone_ID_letters                         <- c()
    for (i_clone in 1:length(sample_clone_ID_unique_numeric)){
        clone_ID_numeric                            <- sample_clone_ID_unique_numeric[i_clone]
        clone_ID_letters                            <- sample_clone_ID_unique_letters[i_clone]
        vec_cell_ID                                 <- which(sample_clone_ID_numeric==clone_ID_numeric)
        sample_clone_ID_letters[vec_cell_ID]        <- clone_ID_letters
    }


    for (row in 1:nrow(table_clone_ID_vs_letters)){
        clone_ID_numeric                            <- table_clone_ID_vs_letters$Clone_ID_number[row]
        clone_ID_letters                            <- table_clone_ID_vs_letters$Clone_ID_letter[row]
        clone_driver_list                           <- unique(genotype_list_driver_map[[clone_ID_numeric]][,1])
        text_report                                 <- paste('Clone ',clone_ID_letters,': ',sep='')
        if (length(clone_driver_list)==0){
            text_report                             <- paste(text_report,'no drivers',sep='')
        }else{
            for (i in 1:length(clone_driver_list)){
                if (i>=2){
                    text_report                     <- paste(text_report,', ',sep='')
                }
                text_report                         <- paste(text_report,driver_library$Gene_ID[clone_driver_list[i]],sep='')
            }
            print(text_report)
        }
    }


#---------------------------------Output package of data from simulation
    output                                  <- list()
    output[[1]]                             <- sample_genotype_profile
    output[[2]]                             <- sample_cell_ID
    output[[3]]                             <- sample_clone_ID
    output[[4]]                             <- sample_clone_ID_letters
    output[[5]]                             <- table_clone_ID_vs_letters
    output[[6]]                             <- sample_time
    return(output)
}
