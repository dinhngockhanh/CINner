BUILD_cn_normal_XX <- function(TABLE_CHROMOSOME_CN_INFO){
#----------------------------------------------------Remove chromosome Y
    TABLE_CHROMOSOME_CN_INFO    <- TABLE_CHROMOSOME_CN_INFO[-which(TABLE_CHROMOSOME_CN_INFO$Chromosome=='Y'),]
#--------------------------------------------------------Build CN matrix
    vec_chrom                   <- c()
    vec_bin                     <- c()
    for (i_chrom in 1:nrow(TABLE_CHROMOSOME_CN_INFO)){
        chrom                   <- TABLE_CHROMOSOME_CN_INFO$Chromosome[i_chrom]
        bin_count               <- TABLE_CHROMOSOME_CN_INFO$Bin_count[i_chrom]
        vec_chrom               <- c(vec_chrom,rep(chrom,bin_count))
        vec_bin                 <- c(vec_bin,(1:bin_count))
    }
    columns                     <- c('Chromosome','Strand','Bin','Allele')

    CN_strand_1                 <- data.frame(vec_chrom,1,vec_bin,'A')
    colnames(CN_strand_1)       <- columns

    CN_strand_2                 <- data.frame(vec_chrom,2,vec_bin,'B')
    colnames(CN_strand_2)       <- columns

    CN_profile                  <- rbind(CN_strand_1,CN_strand_2)
    return(CN_profile)
}
