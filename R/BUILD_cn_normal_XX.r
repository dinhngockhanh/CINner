BUILD_cn_normal_XX <- function(TABLE_CHROMOSOME_CN_INFO){
#----------------------------------------------------Remove chromosome Y
    TABLE_CHROMOSOME_CN_INFO    <- TABLE_CHROMOSOME_CN_INFO[-which(TABLE_CHROMOSOME_CN_INFO$Chromosome=='Y'),]
#--------------------------------------------------------Build CN matrix
    vec_chrom                   <- c()
    vec_bin_start               <- c()
    vec_bin_end                 <- c()
    for (i_chrom in 1:nrow(TABLE_CHROMOSOME_CN_INFO)){
        chrom                   <- TABLE_CHROMOSOME_CN_INFO$Chromosome[i_chrom]
        bin_count               <- TABLE_CHROMOSOME_CN_INFO$Bin_count[i_chrom]


        vec_chrom               <- c(vec_chrom,chrom)
        vec_bin_start           <- c(vec_bin_start,1)
        vec_bin_end             <- c(vec_bin_end,bin_count)
    }
    columns                     <- c('Chromosome','Strand','Bin_start','Bin_end','Allele')

    CN_strand_1                 <- data.frame(vec_chrom,1,vec_bin_start,vec_bin_end,'A')
    colnames(CN_strand_1)       <- columns

    CN_strand_2                 <- data.frame(vec_chrom,2,vec_bin_start,vec_bin_end,'B')
    colnames(CN_strand_2)       <- columns

    CN_profile                  <- rbind(CN_strand_1,CN_strand_2)
    return(CN_profile)
}
