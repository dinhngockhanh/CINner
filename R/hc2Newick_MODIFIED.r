hc2Newick_MODIFIED <- function(hc, flat = TRUE) {
    dist <- 0
    if (is.null(hc$labels)) {
        labels <- seq(along = hc$order)
    } else {
        labels <- hc$labels
    }

    putparenthesis <- function(i) {
        ## recursive function
        ## i: index of row in hc$merge
        j <- hc$merge[i, 1]
        k <- hc$merge[i, 2]

        if (j < 0) {
            left <- labels[-j]
            if (k > 0) {
                dist <- hc$height[i] - hc$height[k]
            } else {
                dist <- hc$height[i]
            }
        } else {
            left <- putparenthesis(j)
        }

        if (k < 0) {
            right <- labels[-k]
            if (j > 0) {
                dist <- hc$height[i] - hc$height[j]
            } else {
                dist <- hc$height[i]
            }
        } else {
            right <- putparenthesis(k)
        }
        if (flat) {
            merged_node_label <- paste("Internal-node-", i, sep = "")
            return(paste("(", left, ":", dist / 2, ",", right, ":", dist / 2, ")", merged_node_label, sep = ""))
            # return(paste("(", left, ":", dist / 2, ",", right, ":", dist / 2, ")", sep = ""))
        } else {
            return(list(left = left, right = right, dist = dist))
        }
    }

    n <- putparenthesis(nrow(hc$merge))
    if (flat) {
        n <- paste(n, ";", sep = "")
    }

    return(n)
}
