#' Dataframe of information about autosome bins
#'
#' @description
#' `BUILD_cn_normal_autosomes` returns a dataframe of information about each autosome, including bin sizes and associated strand/allele.
#' 
#' @param TABLE_CHROMOSOME_CN_INFO A dataframe containing columns "Chromosome", "Bin_count", "Centromere_location". "Chromosome" refers to the chromosome number. "Bin_count" refers to the number of bins in each chromosome. "Centromere_location" refers to the approximate base pair position of the centromere of the chromosome.
#' 
#' @examples
#' 
#' CN_bin_length <- 500000
#' vec_chromosome_name <- c(1:22, "X", "Y")
#' if (!CN_arm_level) {
#'    vec_chromosome_bp <- c(
#'        248956422, 242193529, 198295559, 190214555,
#'        181538259, 170805979, 159345973, 145138636,
#'        138394717, 133797422, 135086622, 133275309,
#'        114364328, 107043718, 101991189, 90338345,
#'        83257441, 80373285, 58617616, 64444167,
#'        46709983, 50818468, 156040895, 57227415
#'    )
#'    vec_centromere_bp <- c(
#'        125, 93.3, 91, 50.4,
#'        48.4, 61, 59.9, 45.6,
#'        49, 40.2, 53.7, 35.8,
#'        17.9, 17.6, 19, 36.6,
#'        24, 17.2, 26.5, 27.5,
#'        13.2, 14.7, 60.6, 10.4
#'    ) * 10^6
#'    vec_bin_count <- ceiling(vec_chromosome_bp / CN_bin_length)
#'    vec_centromere_location <- round(vec_centromere_bp / CN_bin_length)
#' } else {
#'    vec_bin_count <- rep(2, length(vec_chromosome_name))
#'    vec_centromere_location <- rep(1, length(vec_chromosome_name))
#' }
#'
#' TABLE_CHROMOSOME_CN_INFO <- data.frame(vec_chromosome_name, vec_bin_count, vec_centromere_location)
#' columns <- c("Chromosome", "Bin_count", "Centromere_location")
#' colnames(TABLE_CHROMOSOME_CN_INFO) <- columns
#' 
#' CN_matrix <- BUILD_cn_normal_autosomes(TABLE_CHROMOSOME_CN_INFO)
#'
#' @export
BUILD_cn_normal_autosomes <- function(TABLE_CHROMOSOME_CN_INFO) {
    #---------------------------------------------Remove sex chromosomes
    if (length(which(TABLE_CHROMOSOME_CN_INFO$Chromosome %in% c("X", "Y"))) > 0) TABLE_CHROMOSOME_CN_INFO <- TABLE_CHROMOSOME_CN_INFO[-which(TABLE_CHROMOSOME_CN_INFO$Chromosome %in% c("X", "Y")), ]
    #--------------------------------------------------------Build CN matrix
    vec_chrom <- c()
    vec_bin_start <- c()
    vec_bin_end <- c()
    for (i_chrom in 1:nrow(TABLE_CHROMOSOME_CN_INFO)) {
        chrom <- TABLE_CHROMOSOME_CN_INFO$Chromosome[i_chrom]
        bin_count <- TABLE_CHROMOSOME_CN_INFO$Bin_count[i_chrom]


        vec_chrom <- c(vec_chrom, chrom)
        vec_bin_start <- c(vec_bin_start, 1)
        vec_bin_end <- c(vec_bin_end, bin_count)
    }
    columns <- c("Chromosome", "Strand", "Bin_start", "Bin_end", "Allele")

    CN_strand_1 <- data.frame(vec_chrom, 1, vec_bin_start, vec_bin_end, "A")
    colnames(CN_strand_1) <- columns

    CN_strand_2 <- data.frame(vec_chrom, 2, vec_bin_start, vec_bin_end, "B")
    colnames(CN_strand_2) <- columns

    CN_profile <- rbind(CN_strand_1, CN_strand_2)
    return(CN_profile)
}

#' Dataframe of information about chromosome bins
#'
#' @description
#' `BUILD_cn_normal_autosomes` returns a dataframe of information about each chromosome excluding the Y chromosome, including bin sizes and associated strand/allele.
#' 
#' @param TABLE_CHROMOSOME_CN_INFO A dataframe containing columns "Chromosome", "Bin_count", "Centromere_location". "Chromosome" refers to the chromosome number. "Bin_count" refers to the number of bins in each chromosome. "Centromere_location" refers to the approximate base pair position of the centromere of the chromosome.
#' 
#' @examples
#' 
#' CN_bin_length <- 500000
#' vec_chromosome_name <- c(1:22, "X", "Y")
#' if (!CN_arm_level) {
#'    vec_chromosome_bp <- c(
#'        248956422, 242193529, 198295559, 190214555,
#'        181538259, 170805979, 159345973, 145138636,
#'        138394717, 133797422, 135086622, 133275309,
#'        114364328, 107043718, 101991189, 90338345,
#'        83257441, 80373285, 58617616, 64444167,
#'        46709983, 50818468, 156040895, 57227415
#'    )
#'    vec_centromere_bp <- c(
#'        125, 93.3, 91, 50.4,
#'        48.4, 61, 59.9, 45.6,
#'        49, 40.2, 53.7, 35.8,
#'        17.9, 17.6, 19, 36.6,
#'        24, 17.2, 26.5, 27.5,
#'        13.2, 14.7, 60.6, 10.4
#'    ) * 10^6
#'    vec_bin_count <- ceiling(vec_chromosome_bp / CN_bin_length)
#'    vec_centromere_location <- round(vec_centromere_bp / CN_bin_length)
#' } else {
#'    vec_bin_count <- rep(2, length(vec_chromosome_name))
#'    vec_centromere_location <- rep(1, length(vec_chromosome_name))
#' }
#'
#' TABLE_CHROMOSOME_CN_INFO <- data.frame(vec_chromosome_name, vec_bin_count, vec_centromere_location)
#' columns <- c("Chromosome", "Bin_count", "Centromere_location")
#' colnames(TABLE_CHROMOSOME_CN_INFO) <- columns
#' 
#' CN_matrix <- BUILD_cn_normal_XX(TABLE_CHROMOSOME_CN_INFO)
#'
#' @export
BUILD_cn_normal_XX <- function(TABLE_CHROMOSOME_CN_INFO) {
    #---------------------------------------------Remove sex chromosomes
    if (length(which(TABLE_CHROMOSOME_CN_INFO$Chromosome %in% c("Y"))) > 0) TABLE_CHROMOSOME_CN_INFO <- TABLE_CHROMOSOME_CN_INFO[-which(TABLE_CHROMOSOME_CN_INFO$Chromosome %in% c("Y")), ]
    # if (length(which(TABLE_CHROMOSOME_CN_INFO$Chromosome %in% c("X", "Y"))) > 0) TABLE_CHROMOSOME_CN_INFO <- TABLE_CHROMOSOME_CN_INFO[-which(TABLE_CHROMOSOME_CN_INFO$Chromosome %in% c("X", "Y")), ]
    #--------------------------------------------------------Build CN matrix
    vec_chrom <- c()
    vec_bin_start <- c()
    vec_bin_end <- c()
    for (i_chrom in 1:nrow(TABLE_CHROMOSOME_CN_INFO)) {
        chrom <- TABLE_CHROMOSOME_CN_INFO$Chromosome[i_chrom]
        bin_count <- TABLE_CHROMOSOME_CN_INFO$Bin_count[i_chrom]


        vec_chrom <- c(vec_chrom, chrom)
        vec_bin_start <- c(vec_bin_start, 1)
        vec_bin_end <- c(vec_bin_end, bin_count)
    }
    columns <- c("Chromosome", "Strand", "Bin_start", "Bin_end", "Allele")

    CN_strand_1 <- data.frame(vec_chrom, 1, vec_bin_start, vec_bin_end, "A")
    colnames(CN_strand_1) <- columns

    CN_strand_2 <- data.frame(vec_chrom, 2, vec_bin_start, vec_bin_end, "B")
    colnames(CN_strand_2) <- columns

    CN_profile <- rbind(CN_strand_1, CN_strand_2)
    return(CN_profile)
}
