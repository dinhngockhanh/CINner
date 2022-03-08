BUILD_cn_normal_XX <- function(TABLE_CHROMOSOME_CN_INFO){
#----------------------------------------------------Remove chromosome Y
    TABLE_CHROMOSOME_CN_INFO    <- TABLE_CHROMOSOME_CN_INFO[-which(TABLE_CHROMOSOME_CN_INFO$Chromosome=='Y'),]

    print(TABLE_CHROMOSOME_CN_INFO)


}
