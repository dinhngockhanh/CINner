#==================================COMPUTE THE SELECTION RATE OF A CLONE
SIMULATOR_FULL_PHASE_1_selection_rate <- function(driver_count,driver_map,ploidy_chrom,ploidy_block,ploidy_allele) {
#-------------------------Cell is not viable if losing whole chromosomes
#------------------------------------or exceeding maximum average ploidy
    vec_CN_all                      <- c()
    for (chrom in 1:N_chromosomes){
        no_strands                  <- ploidy_chrom[chrom]
        if (no_strands<=0){
            clone_selection_rate    <- 0
            return(clone_selection_rate)
        }
        vec_CN                      <- rep(0,vec_CN_block_no[chrom])
        for (strand in 1:no_strands){
            vec_CN                  <- vec_CN+ploidy_block[[chrom]][[strand]]
        }
        if (max(vec_CN)==0){
            clone_selection_rate    <- 0
            return(clone_selection_rate)
        }
        vec_CN_all                  <- c(vec_CN_all,vec_CN)
    }
    if (mean(vec_CN_all)>bound_average_ploidy){
        clone_selection_rate        <- 0
        return(clone_selection_rate)
    }
#---------Cell is not viable is exceeding maximum length of homozygosity
    L_homozygosity                  <- 0
    for (chrom in 1:N_chromosomes){
        no_strands                  <- ploidy_chrom[chrom]
        vec_allele_1                <- rep(0,vec_CN_block_no[chrom])
        vec_allele_2                <- rep(0,vec_CN_block_no[chrom])

        for (strand in 1:no_strands){
            N_rows                  <- nrow(ploidy_allele[[chrom]][[strand]])
            for (row in 1:N_rows){
                vec_loc_1               <- which(ploidy_allele[[chrom]][[strand]][row,]==1)
                vec_allele_1[vec_loc_1] <- 1

                vec_loc_2               <- which(ploidy_allele[[chrom]][[strand]][row,]==2)
                vec_allele_2[vec_loc_2] <- 1
            }
        }

        vec_homozygosity            <- which(vec_allele_1==0 | vec_allele_2==0)
        L_homozygosity              <- L_homozygosity+length(vec_homozygosity)

        if (L_homozygosity>bound_homozygosity){
            clone_selection_rate        <- 0
            return(clone_selection_rate)
        }
    }


#-------------------Cell is not viable if exceeding maximum driver count
    if (driver_count>0){
        driver_count_unique         <- unique(driver_map[,1])
        if (length(driver_count_unique)>bound_driver){
            clone_selection_rate    <- 0
            return(clone_selection_rate)
        }
    }
#----------If driver library is empty, then viable cells have sel rate 1
    if (nrow(driver_library)==0){
        clone_selection_rate        <- 1
        return(clone_selection_rate)
    }
#------------------------------Compute selection rates for viable clones
    driver_library_copy                             <- driver_library
    driver_library_copy$Copy_WT                     <- 0
    driver_library_copy$Copy_MUT                    <- 0
#   Find WT and MUT allele counts for each driver
    for (i_driver in 1:nrow(driver_library_copy)){
        chrom                                       <- driver_library_copy$Chromosome[i_driver]
        block                                       <- driver_library_copy$Bin[i_driver]
        no_strands                                  <- ploidy_chrom[chrom]
        driver_copy                                 <- 0
        # if (no_strands>0){
            for (strand in 1:no_strands){
                driver_copy                         <- driver_copy+ploidy_block[[chrom]][[strand]][block]
            }
        # }
        driver_library_copy$Copy_WT[i_driver]       <- driver_copy
    }
    if (driver_count>=1){
        for (i_driver in 1:driver_count){
            driver_ID                               <- driver_map[i_driver,1]
            driver_library_copy$Copy_MUT[driver_ID] <- driver_library_copy$Copy_MUT[driver_ID]+1
            driver_library_copy$Copy_WT[driver_ID]  <- driver_library_copy$Copy_WT[driver_ID]-1
        }
    }
#   Compute selection rate
    clone_selection_rate                            <- prod(driver_library_copy$s_rate_WT^driver_library_copy$Copy_WT)*prod(driver_library_copy$s_rate_MUT^driver_library_copy$Copy_MUT)
    return(clone_selection_rate)
}
