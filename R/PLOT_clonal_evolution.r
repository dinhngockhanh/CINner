#=====================================PLOT CLONAL EVOLUTION AS FISH PLOT
PLOT_clonal_evolution <- function(package_simulation,vec_time_plot,unit){
    if (unit=='year'){
        vec_time_plot               <- 365*vec_time_plot
    }
#---------------------------------------------Input the clonal evolution
    package_clonal_evolution                        <- package_simulation[[1]]
    evolution_origin                                <- package_clonal_evolution[[11]]
    evolution_traj_time                             <- package_clonal_evolution[[13]]
    evolution_traj_clonal_ID                        <- package_clonal_evolution[[15]]
    evolution_traj_population                       <- package_clonal_evolution[[16]]
#---------------------------------------------Input the clonal phylogeny
    package_sample_phylogeny                        <- package_simulation[[3]]
    package_clone_phylogeny                         <- package_sample_phylogeny[[4]]
    clone_phylogeny_labels                          <- package_clone_phylogeny[[1]]
    clone_phylogeny_origin                          <- package_clone_phylogeny[[3]]
    clone_phylogeny_genotype                        <- package_clone_phylogeny[[4]]
    clone_hclust_nodes                              <- package_clone_phylogeny[[7]]
    clone_hclust_merge                              <- package_clone_phylogeny[[8]]

    N_clones                                        <- length(clone_phylogeny_labels)
#----------------------------Build the genotype list for each clone node
    clone_phylogeny_all_genotypes                                       <- vector('list',length=length(clone_phylogeny_genotype))
#   Initialize genotype lists for clone leaves
    for (clone in length(clone_phylogeny_genotype):(length(clone_phylogeny_genotype)-N_clones+1)){
        all_genotypes                                                   <- clone_phylogeny_genotype[clone]
        while (all_genotypes[1]!=0){
            ancestor_genotype                                           <- evolution_origin[all_genotypes[1]]
            all_genotypes                                               <- c(ancestor_genotype,all_genotypes)
        }
        clone_phylogeny_all_genotypes[[clone]]                          <- all_genotypes
    }
#   Update genotype lists with information from clone hclust
    for (clone_hclust_mother_node in 1:nrow(clone_hclust_merge)){
#       Find mother clone node's index in phylogeny our style
        clone_phylogeny_mother_node                                     <- which(is.element(clone_hclust_nodes,clone_hclust_mother_node))
#       Find daughter clone nodes' indices in phylogeny our style
        vec_clone_hclust_daughter_nodes                                 <- clone_hclust_merge[clone_hclust_mother_node,]
        vec_clone_phylogeny_daughter_nodes                              <- which(is.element(clone_hclust_nodes,vec_clone_hclust_daughter_nodes))
#       Update genotype lists for the 3 clone nodes
        clone_phylogeny_daughter_node_1                                 <- vec_clone_phylogeny_daughter_nodes[[1]]
        clone_phylogeny_daughter_node_2                                 <- vec_clone_phylogeny_daughter_nodes[[2]]
        mother_genotypes                                                <- intersect(clone_phylogeny_all_genotypes[[clone_phylogeny_daughter_node_1]],clone_phylogeny_all_genotypes[[clone_phylogeny_daughter_node_2]])
        clone_phylogeny_all_genotypes[[clone_phylogeny_mother_node]]    <- mother_genotypes
        clone_phylogeny_all_genotypes[[clone_phylogeny_daughter_node_1]]<- setdiff(clone_phylogeny_all_genotypes[[clone_phylogeny_daughter_node_1]],mother_genotypes)
        clone_phylogeny_all_genotypes[[clone_phylogeny_daughter_node_2]]<- setdiff(clone_phylogeny_all_genotypes[[clone_phylogeny_daughter_node_2]],mother_genotypes)
    }
#---------------------------------Find clonal populations as time series
    table_clonal_populations                                            <- matrix(0,nrow=length(clone_phylogeny_all_genotypes),ncol=length(vec_time_plot))
    vec_total_populations                                               <- rep(0,length=length(vec_time_plot))
    for (col in 1:length(vec_time_plot)){
        time                                                            <- vec_time_plot[col]
        loc                                                             <- which.min(abs(evolution_traj_time-time))
        vec_clonal_ID                                                   <- evolution_traj_clonal_ID[[loc]]
        vec_clonal_population                                           <- evolution_traj_population[[loc]]
        total_clonal_population                                         <- sum(vec_clonal_population)
        for (row in 1:length(clone_phylogeny_all_genotypes)){
            vec_loc                                                     <- which(is.element(vec_clonal_ID,clone_phylogeny_all_genotypes[[row]]))
            if (length(vec_loc)==0){
                next
            }
            table_clonal_populations[row,col]                           <- sum(vec_clonal_population[vec_loc])
        }
        vec_total_populations[col]                                      <- total_clonal_population
    }
#--------------------------------------------------Find clonal parentage
    vec_clonal_parentage                                                <- clone_phylogeny_origin
#----------------------------------------Add a clone for other genotypes
    table_clonal_populations                                            <- rbind(rep(0,length=length(vec_time_plot)),table_clonal_populations)
    for (col in 1:length(vec_time_plot)){
        table_clonal_populations[1,col]                                 <- vec_total_populations[col]-sum(table_clonal_populations[,col])
    }
    vec_clonal_parentage                                                <- c(0,(vec_clonal_parentage+1))
#----------------------------------------------Remove unnecessary clones
#-----------------------------------------i.e. clones that are always 0%
#   Find unnecessary clones
    vec_unnecessary_clones                                              <- c()
    for (clone in 1:nrow(table_clonal_populations)){
        if (all(table_clonal_populations[clone,]==0)){
            vec_unnecessary_clones                                      <- c(vec_unnecessary_clones,clone)
        }
    }
#   Rewire clones to new mother clones
    for (clone_daughter in 1:length(vec_clonal_parentage)){
        clone_mother                                                    <- vec_clonal_parentage[clone_daughter]
        while ((clone_mother!=0) & (is.element(clone_mother,vec_unnecessary_clones))){
            clone_mother                                                <- vec_clonal_parentage[clone_mother]
        }
        vec_clonal_parentage[clone_daughter]                            <- clone_mother
    }
#   Remove unnecessary clones from existence
    table_clonal_populations                                            <- table_clonal_populations[-vec_unnecessary_clones,]
    vec_clonal_parentage                                                <- vec_clonal_parentage[-vec_unnecessary_clones]
#-----------------Scale the clonal populations to match total population
    max_total_population                                                <- max(vec_total_populations)
    for (col in 1:length(vec_time_plot)){
#       Total population size can only go up to 90% to have leeway with
#       numerical errors
        table_clonal_populations[,col]                                  <- 49*table_clonal_populations[,col]/max_total_population
    }
#---------------Conform clonal populations to nested format of fish plot
    table_clonal_populations_tmp                                        <- table_clonal_populations
    for (clone_daughter in 1:length(vec_clonal_parentage)){
        clone_mother                                                    <- vec_clonal_parentage[clone_daughter]
        while (clone_mother>0){
#           Add in 0.001 to offset potential numerical errors to make
#           sure mother clone always has more than sum of daughter clones
            table_clonal_populations[clone_mother,]                     <- table_clonal_populations[clone_mother,]+table_clonal_populations_tmp[clone_daughter,]+0.001
            clone_mother                                                <- vec_clonal_parentage[clone_mother]
        }
    }
#----------Final check to make sure total population doesn't exceed 100%
    for (row in 1:nrow(table_clonal_populations)){
        vec_row                                 <- table_clonal_populations[row,]
        vec_fix                                 <- which(vec_row>100)
        vec_row[vec_fix]                        <- 100
        table_clonal_populations[row,]          <- vec_row
    }
#----------------------------------------------Plot the clonal evolution
    if (unit=='year'){
        vec_time_plot               <- vec_time_plot/365
    }

    fish    <- createFishObject(table_clonal_populations,vec_clonal_parentage,timepoints=vec_time_plot)

    fish    <- layoutClones(fish)

    # fishPlot(fish,shape="spline",title.btm="Sample1",cex.title=0.5)

    fishPlot(fish,shape="spline")

    # vlines=c(0,150),vlab=c("day 0","day 150"))


}
