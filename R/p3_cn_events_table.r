#' @export
p3_cn_events_table <- function(simulation) {
    # hclust_internal_genotypes <- simulation$sample_phylogeny$package_cell_phylogeny_hclust_extra$hclust_internal_genotypes
    hclust_CN_events <- simulation$sample_phylogeny$package_cell_phylogeny_hclust_extra$hclust_CN_events

    sample_cell_ID <- simulation$sample$sample_cell_ID
    sample_internal_node_ID <- simulation$sample$sample_internal_node_ID

    evolution_origin <- simulation$clonal_evolution$evolution_origin
    evolution_genotype_changes <- simulation$clonal_evolution$evolution_genotype_changes
    #---------------------------------------Add column for mother row ID
    hclust_CN_events$Mother_ID <- ""
    for (row in 1:nrow(hclust_CN_events)) {
        hclust_mother <- hclust_CN_events$Mother_hclust[row]
        if (hclust_mother > 0) {
            hclust_CN_events$Mother_ID[row] <- sample_internal_node_ID[hclust_mother]
        } else {
            if (hclust_mother < 0) {
                hclust_CN_events$Mother_ID[row] <- sample_cell_ID[-hclust_mother]
            } else {
                hclust_CN_events$Mother_ID[row] <- paste("Initial-clone-", hclust_CN_events$Mother_genotype[row], sep = "")
            }
        }
    }
    #-------------------------------------Add column for daughter row ID
    hclust_CN_events$Daughter_ID <- ""
    for (row in 1:nrow(hclust_CN_events)) {
        hclust_daughter <- hclust_CN_events$Daughter_hclust[row]
        if (hclust_daughter > 0) {
            hclust_CN_events$Daughter_ID[row] <- sample_internal_node_ID[hclust_daughter]
        } else {
            if (hclust_daughter < 0) {
                hclust_CN_events$Daughter_ID[row] <- sample_cell_ID[-hclust_daughter]
            } else {
                hclust_CN_events$Daughter_ID[row] <- paste("Initial-clone-", hclust_CN_events$Daughter_genotype[row], sep = "")
            }
        }
    }
    #------------------------------Remove columns for row hclust indices
    hclust_CN_events$Mother_hclust <- NULL
    hclust_CN_events$Daughter_hclust <- NULL
    #-------------------------------------------Add column for CN events
    #--------------------------from mother genotype to daughter genotype
    hclust_CN_events$CN_events <- ""
    for (row in 1:nrow(hclust_CN_events)) {
        genotype_start <- hclust_CN_events$Mother_genotype[row]
        genotype_end <- hclust_CN_events$Daughter_genotype[row]
        edge_CN_events <- ""
        genotype_current <- genotype_end
        while (genotype_current != genotype_start) {
            genotype_CN_events <- evolution_genotype_changes[[genotype_current]]
            for (event in length(genotype_CN_events):1) {
                event_info <- genotype_CN_events[[event]]
                event_type <- event_info[1]
                if (edge_CN_events != "") {
                    edge_CN_events <- paste(";", edge_CN_events, sep = "")
                }
                if (event_type == "new-driver") {
                    edge_CN_events <- paste("Driver-mutation:", driver_library$Gene_ID[event_info[2]], edge_CN_events, sep = "")
                }
                if (event_type == "whole-genome-duplication") {
                    edge_CN_events <- paste("whole-genome-duplication", edge_CN_events, sep = "")
                }
                if (event_type == "missegregation") {
                    if (event_info[4] == 1) {
                        edge_CN_events <- paste("missegregation-gain:chrom", event_info[2], edge_CN_events, sep = "")
                    } else {
                        edge_CN_events <- paste("missegregation-loss:chrom", event_info[2], edge_CN_events, sep = "")
                    }
                }
                if (event_type == "focal-deletion") {
                    edge_CN_events <- paste("focal-deletion:chrom", event_info[2], "blocks", event_info[4], "-", event_info[5], edge_CN_events, sep = "")
                }
                if (event_type == "focal-amplification") {
                    edge_CN_events <- paste("focal-amplification:chrom", event_info[2], "blocks", event_info[4], "-", event_info[5], edge_CN_events, sep = "")
                }
                if (event_type == "cnloh-terminal") {
                    edge_CN_events <- paste("terminal-CNLOH:chrom", event_info[2], "blocks", event_info[5], "-", event_info[6], edge_CN_events, sep = "")
                }
                if (event_type == "cnloh-interstitial") {
                    edge_CN_events <- paste("interstitial-CNLOH:chrom", event_info[2], "blocks", event_info[5], "-", event_info[6], edge_CN_events, sep = "")
                }
                if (event_type == "chromosome-arm-missegregation") {
                    if (event_info[5] == 1) {
                        edge_CN_events <- paste("arm-missegregation-gain:chrom", event_info[2], "arm", event_info[4], edge_CN_events, sep = "")
                    } else {
                        edge_CN_events <- paste("arm-missegregation-loss:chrom", event_info[2], "arm", event_info[4], edge_CN_events, sep = "")
                    }
                }
            }
            genotype_current <- evolution_origin[genotype_current]
        }
        hclust_CN_events$CN_events[row] <- edge_CN_events
    }
    #-----------------------------------Remove columns for row genotypes
    hclust_CN_events$Mother_genotype <- NULL
    hclust_CN_events$Daughter_genotype <- NULL
    #------------------------------------------Output table of CN events
    # print(hclust_CN_events)
    simulation$sample_phylogeny$package_cell_phylogeny_hclust_extra$hclust_CN_events <- hclust_CN_events
    return(simulation)
}
