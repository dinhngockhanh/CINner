#' @export
p0_write_cn_as_wig <- function(filename, noisy_cn_profiles_long, cell_ID) {
    library(readr)
    library(data.table)
    #-----------------------------------Extract readcounts for this cell
    noisy_cn_profiles_long <- noisy_cn_profiles_long[noisy_cn_profiles_long$cell_id == cell_ID, ]
    #--------------------------------------Create dataframe for WIG file
    # df_wig <- data.frame(matrix(0, nrow = 0, ncol = 1))
    # colnames(df_wig) <- "track type=wiggle_0"
    # list_chr <- unique(noisy_cn_profiles_long$chr)
    # for (chrom in 1:length(list_chr)) {
    #     chrom_ID <- list_chr[chrom]
    #     df_wig[nrow(df_wig) + 1, 1] <- paste("fixedStep chrom=", chrom_ID, " start=", 1, " step=", as.integer(size_CN_block_DNA), " span=", as.integer(size_CN_block_DNA), sep = "")
    #     vec_loc <- which(noisy_cn_profiles_long$chr == chrom_ID)
    #     for (loc in 1:length(vec_loc)) {
    #         df_wig[nrow(df_wig) + 1, 1] <- noisy_cn_profiles_long$reads[vec_loc[loc]]
    #     }
    # }




    list_chr <- unique(noisy_cn_profiles_long$chr)
    df_wig_list <- vector("list", length = length(list_chr))
    for (chrom in 1:length(list_chr)) {
        chrom_ID <- list_chr[chrom]
        vec_loc <- which(noisy_cn_profiles_long$chr == chrom_ID)
        df_wig_list[[chrom]] <- data.frame(c(paste0("fixedStep chrom=", chrom_ID, " start=", 1, " step=", as.integer(size_CN_block_DNA), " span=", as.integer(size_CN_block_DNA)), noisy_cn_profiles_long$reads[vec_loc]))
        colnames(df_wig_list[[chrom]]) <- "track type=wiggle_0"
    }
    df_wig <- rbindlist(df_wig_list)
    #----------------------------------------Print dataframe to WIG file
    write_csv(df_wig, file = filename)
}
