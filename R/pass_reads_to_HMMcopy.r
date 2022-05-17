# library("HMMcopy")
# library("GenomicRanges")
# library("rtracklayer")
#
# # Preprocessing steps (we need map.wig and gc.wig files).
# # Hopefully the simulated data can make sense with hg19; otherwise we need a simulated reference genome to replace hg19
# # 0. Get hg19 reference genome or the corresponding simulated reference genome
# # 1. Create BigWig for mappability using http://bioinformatics-ca.github.io/bioinformatics_for_cancer_genomics_mod5_lab_data_prep_2015/
# # 2. Convert BigWig to Wig file
# # 3. Create gc out of reference genome using https://github.com/shahcompbio/hmmcopy_utils
#
#
# ################################################# import, convert and order data (below is an example of dlp data)
#
# # the file is from https://zenodo.org/record/3445364#.YRZEWjYzbRZ
# # take raw column from the file
# dlp_counts <- read.csv("./ov2295_cell_cn.csv")
#
# grange_counts <- makeGRangesFromDataFrame(dlp_counts[, 4:7], keep.extra.columns = TRUE)
#
# # the wig files should be available from HMMcopy repository
# gcdata <- wigToRangedData("./gc.wig")
# gcdata$chr <- factor(gcdata$chr, levels = c(
#     "1", "2", "3", "4", "5", "6", "7", "8", "9", "10",
#     "11", "12", "13", "14", "15", "16", "17", "18", "19",
#     "20", "21", "22", "X", "Y"
# ))
# gc_ordered <- gcdata[order(gcdata$chr), ]
#
# mapdata <- wigToRangedData("./map.wig")
# mapdata$chr <- factor(mapdata$chr, levels = c(
#     "1", "2", "3", "4", "5", "6", "7", "8", "9", "10",
#     "11", "12", "13", "14", "15", "16", "17", "18", "19",
#     "20", "21", "22", "X", "Y"
# ))
# map_ordered <- mapdata[order(mapdata$chr), ]
#
#
# ################################################# Some tests
# length(which(unique(round(dlp_counts$copy)) < 0))
# # result: [1] 54
# length(which(round(dlp_counts$copy) < 0))
# # result: [1] 102
# length((unique(round(dlp_counts$copy))))
# # result: [1] 946
# sum(is.na(dlp_counts$copy))
# # result: [1] 2066024
#
#
# ################################################# For each cell perform HMMcopy
#
# cell <- data.frame(matrix(NA, nrow = 6206, ncol = 6)) # 6206 is the number of bins in the data
# c <- list()
# cell_num <- 1
#
# new_cell <- grange_counts[1:(1 + 6205)]
# new_cell$gc <- gc_ordered$value
# new_cell$map <- map_ordered$value
#
# for (r in seq(1, length(grange_counts), 6206)) {
#     new_c <- grange_counts[r:(r + 6205)]
#     new_c$gc <- gc_ordered$value
#     new_c$map <- map_ordered$value
#     corr_copy <- correctReadcount(new_c)
#     segmented_copy <- HMMsegment(corr_copy)
#     c[[cell_num]] <- segmented_copy
#     cell_num <- cell_num + 1
# }
#
# write.csv(c, file = paste0("./hmmcopy_result_segmented_copy.csv"))
