#' @export
p0_write_gc_map_as_wig <- function(filename_gc, filename_map) {
    #-----------------------Initialize data frame for GC and mappability
    vec_chr <- c()
    vec_start <- c()
    vec_end <- c()
    for (chrom in 1:N_chromosomes) {
        N_bins <- vec_CN_block_no[chrom]
        vec_chr <- c(vec_chr, rep(vec_chromosome_id[chrom], length = N_bins))
        vec_start <- c(vec_start, (size_CN_block_DNA * (0:(N_bins - 1)) + 1))
        vec_end <- c(vec_end, size_CN_block_DNA * (1:N_bins))
    }
    df_gc <- data.frame(chr = vec_chr, start = vec_start, end = vec_end)
    #--------------------------------Input values for GC and mappability
    df_gc$gc <- -1
    df_gc$map <- -1
    for (row in 1:nrow(table_gc)) {
        val_chr <- table_gc$chr[row]
        val_start <- table_gc$start[row]
        val_end <- table_gc$end[row]
        val_gc <- table_gc$gc[row]
        val_map <- table_gc$map[row]
        loc <- which(df_gc$chr == val_chr & df_gc$start == val_start & df_gc$end == val_end)
        df_gc$gc[loc] <- val_gc
        df_gc$map[loc] <- val_map
    }
    #-------------------------------Save GC and mappability as WIG files
    write("track type=wiggle_0", file = filename_gc, append = FALSE, sep = "\n")
    write("track type=wiggle_0", file = filename_map, append = FALSE, sep = "\n")
    list_chr <- unique(df_gc$chr)
    for (chrom in 1:length(list_chr)) {
        chrom_ID <- list_chr[chrom]
        write(paste("fixedStep chrom=", chrom_ID, " start=", 1, " step=", as.integer(size_CN_block_DNA), " span=", as.integer(size_CN_block_DNA), sep = ""), file = filename_gc, append = TRUE, sep = "\n")
        write(paste("fixedStep chrom=", chrom_ID, " start=", 1, " step=", as.integer(size_CN_block_DNA), " span=", as.integer(size_CN_block_DNA), sep = ""), file = filename_map, append = TRUE, sep = "\n")
        vec_loc <- which(df_gc$chr == chrom_ID)
        for (loc in 1:length(vec_loc)) {
            write(df_gc$gc[vec_loc[loc]], file = filename_gc, append = TRUE, sep = "\n")
            write(df_gc$map[vec_loc[loc]], file = filename_map, append = TRUE, sep = "\n")
        }
    }
}
