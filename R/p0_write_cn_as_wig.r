#' @export
p0_write_cn_as_wig <- function(filename, noisy_cn_profiles_long, cell_ID) {
    #-----------------------------------Extract readcounts for this cell
    noisy_cn_profiles_long <- noisy_cn_profiles_long[noisy_cn_profiles_long$cell_id == cell_ID, ]
    #--------------------------------------------Initialize the WIG file
    write("track type=wiggle_0", file = filename, append = FALSE, sep = "\n")
    #---------------------------Append the WIG file with each chromosome
    list_chr <- unique(noisy_cn_profiles_long$chr)
    for (chrom in 1:length(list_chr)) {
        chrom_ID <- list_chr[chrom]
        write(paste("fixedStep chrom=", chrom_ID, " start=", 1, " step=", as.integer(size_CN_block_DNA), " span=", as.integer(size_CN_block_DNA), sep = ""), file = filename, append = TRUE, sep = "\n")
        vec_loc <- which(noisy_cn_profiles_long$chr == chrom_ID)
        for (loc in 1:length(vec_loc)) {
            write(noisy_cn_profiles_long$reads[vec_loc[loc]], file = filename, append = TRUE, sep = "\n")
        }
    }
}
