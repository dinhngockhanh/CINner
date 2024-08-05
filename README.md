#   Simulating chromosomal instability during tumorigenesis

CINner is a tool to simulate cancer evolution as a selective model driven by mutations and Copy Number Aberrations (CNAs).
It is designed to reproduce data from single-cell DNA sequencing such as Direct Library Preperation+ (DLP+).

## Installation

CINner requires some external libraries:

```{r}
BiocManager::install("ggtree")
BiocManager::install("HMMcopy")
devtools::install_github("shahcompbio/signals")
devtools::install_github("chrisamiller/fishplot")
```

You can then install CINner with the following command:

```{r}
devtools::install_github("dinhngockhanh/CINner")
```

## Input data

`CINner` takes as input several `csv` files that dictate the evolution model. The user can either create them from scratch or by using a suite of R functions. The prefix `*` in the filenames must match the model name to be given later to the simulator.

### Model variables

The file `*-input-variables.csv` contains the model variables, including the beginning and final time for the simulations, cell lifespan, size of CN bin, rate of driver mutations and probabilities of CN events.

### Sampling information

The file `*-input-sampling.csv` contains the time where each sample is taken, cell count in each sample, and a name for the sample.

### Chromosome bin count and centromere location

The file `*-input-copy-number-blocks.csv` contains the length of each chromosome as number of bins, as well as the location of its centromere also as bin index.

### Population dynamics

CINner is designed so that the total cell population follows a given dynamics. This is given in the file `*-input-population-dynamics.csv`, which lists the time points and desired population sizes at those time points. Within two time points, the simulator interpolates to find the desired population size.

### Driver gene library.

The file `*-input-cancer-genes.csv` contains information of driver genes to be included in the model. Each gene should be given a role (TSG/ONCOGENE) and the location on the genome (chromosome & bin). The selective strength s for the gene is encoded in `Selective_strength`; larger s means the mutated driver is more likely to be fixated and its fixation takes less time. There are two selection rates for each gene, `s_rate_WT` for the WT allele and `s_rate_MUT` for the mutated allele. These rates are computed from `Selective_strength` according to the following formula:
- For TSG's, s_rate_WT=1/(1+s) and s_rate_MUT=1
- For oncogenes, s_rate_WT=alpha and s_rate_MUT=alpha*(1+s)
where alpha is such that the product of all s_rate_WT is 1.

### Initial population

The information about the initial cell population is given in two files. The file `*-input-initial-others.csv` lists the clones in the population, their cell counts, and mutated drivers. The file `*-input-initial-cn-profiles.csv` specifies the CN profile of each clone. For each clone and for each chromosome, we need to specify the allele identity of each strand.

## Run the simulator

There are several stages for how far the simulator will do, depending on what the user requires from the simulations. The variable `stage_final` dictates this:
- `stage_final`=1 only produces the clonal evolution.
- `stage_final`=2 samples cells as specified above, and outputs their CN profiles as a `rda` file.
- `stage_final`=3 further builds the phylogeny for these cells, and outputs this as a `rda` file. The phylogeny is built based on information from the clonal evolution.

`N_clones_min` and `N_clones_max` are optional variables in case we want to limit the number of clones in the simulation.

To produce one simulation for model `*`, we can call
```r
simulation <- SIMULATOR_FULL_PROGRAM_one_simulation(model=*,stage_final=stage_final,N_clones_min=N_clones_min,N_clones_max=N_clones_max)
```

See the vignettes for further analyses that can be performed from `simulation`.
