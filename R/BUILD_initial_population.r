BUILD_initial_population <- function(MODEL_VARIABLES    = list(),
                                     cell_count         = 1,
                                     CN_arm             = cbind(c(rep('A-A',length=23),''),c(rep('B-B',length=23),'')),
                                     CN_focal           = list(),
                                     drivers            = list()){
#------------------------------Update initial clones - other information
    DRIVERS                                                 <- ''
    if (length(drivers)>0){
        for (i in 1:length(drivers)){
            strand                                          <- drivers[[i]][[1]]
            unit                                            <- drivers[[i]][[2]]
            driver_ID                                       <- drivers[[i]][[3]]
            if (DRIVERS!=''){
                DRIVERS                                     <- paste(DRIVERS,';',sep='')
            }
            DRIVERS                                         <- paste(DRIVERS,driver_ID,'_strand',strand,'_unit',unit,sep='')
        }
    }
    if (is.null(MODEL_VARIABLES$initial_others)){
        columns                                             <- c('Clone','Cell_count','Drivers')

        TABLE_INITIAL_OTHERS                                <- data.frame(1,cell_count,DRIVERS)
        colnames(TABLE_INITIAL_OTHERS)                      <- columns
        I_clone                                             <- 1
    }else{
        TABLE_INITIAL_OTHERS                                <- MODEL_VARIABLES$initial_others
        I_clone                                             <- nrow(TABLE_INITIAL_OTHERS)+1

        TABLE_INITIAL_OTHERS[I_clone,]                      <- c(I_clone,cell_count,DRIVERS)
    }
#------------------------------------Update initial clones - CN profiles
    TABLE_CHROMOSOME_CN_INFO                                <- MODEL_VARIABLES$cn_info
#   Initialize table of alleles if necessary
    if (is.null(MODEL_VARIABLES$initial_cn)){
        vec_chrom                                           <- c()
        vec_bin                                             <- c()
        for (chrom in 1:nrow(TABLE_CHROMOSOME_CN_INFO)){
            bin_count                                       <- TABLE_CHROMOSOME_CN_INFO$Bin_count[chrom]
            vec_chrom                                       <- c(vec_chrom,rep(TABLE_CHROMOSOME_CN_INFO$Chromosome[chrom],bin_count))
            vec_bin                                         <- c(vec_bin,(1:bin_count))
        }
        columns                                             <- c('Chromosome','Bin')
        TABLE_INITIAL_CN                                    <- data.frame(vec_chrom,vec_bin)
        colnames(TABLE_INITIAL_CN)                          <- columns
    }else{
        TABLE_INITIAL_CN                                    <- MODEL_VARIABLES$initial_cn
    }
#---Build the CN profile for the new clone - start from arm level
    max_no_strands                                          <- ncol(CN_arm)
    for (strand in 1:max_no_strands){
        TABLE_INITIAL_CN[paste('Clone_',I_clone,'_strand_',strand,sep='')]  <- 'NA'

    }
    for (i_chrom in 1:nrow(TABLE_CHROMOSOME_CN_INFO)){
        chrom                                               <- TABLE_CHROMOSOME_CN_INFO$Chromosome[i_chrom]
        centromere                                          <- TABLE_CHROMOSOME_CN_INFO$Centromere_location[i_chrom]
        for (strand in 1:max_no_strands){
            col                                             <- which(colnames(TABLE_INITIAL_CN)==paste('Clone_',I_clone,'_strand_',strand,sep=''))
            strand_allele                                   <- CN_arm[i_chrom,strand]
            if (strand_allele==''){
                next
            }

            strand_allele_left                              <- sub('-.*','',strand_allele)
            if (strand_allele_left!=''){
                rows                                        <- which((TABLE_INITIAL_CN$Chromosome==chrom)&(TABLE_INITIAL_CN$Bin<=centromere))
                TABLE_INITIAL_CN[rows,col]                  <- strand_allele_left
            }

            strand_allele_right                             <- sub('.*-','',strand_allele)
            if (strand_allele_right!=''){
                rows                                        <- which((TABLE_INITIAL_CN$Chromosome==chrom)&(TABLE_INITIAL_CN$Bin>centromere))
                TABLE_INITIAL_CN[rows,col]                  <- strand_allele_right
            }
        }
    }
#---Build the CN profile for the new clone - continue with local regions
    if (length(CN_focal)>0){
        for (i_focal in 1:length(CN_focal)){
            focal_event                                     <- CN_focal[[i_focal]]
            chrom                                           <- focal_event[[1]]
            strand                                          <- focal_event[[2]]
            bin_start                                       <- focal_event[[3]]
            bin_end                                         <- focal_event[[4]]
            allele                                          <- focal_event[[5]]
            if (allele==''){
                allele                                      <- 'NA'
            }
            rows                                            <- which((TABLE_INITIAL_CN$Chromosome==chrom)&(TABLE_INITIAL_CN$Bin>=bin_start)&(TABLE_INITIAL_CN$Bin<=bin_end))
            col                                             <- which(colnames(TABLE_INITIAL_CN)==paste('Clone_',I_clone,'_strand_',strand,sep=''))
            TABLE_INITIAL_CN[rows,col]                      <- allele
        }
    }
#----------------------------------------Output the model variable files
    MODEL_VARIABLES$initial_others                          <- TABLE_INITIAL_OTHERS
    MODEL_VARIABLES$initial_cn                              <- TABLE_INITIAL_CN
    return(MODEL_VARIABLES)
}
