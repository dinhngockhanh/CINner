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

## Run CINner

See the vignettes for further analyses that can be performed with `CINner`.

For more detailed examples and usage, please refer to `vignettes/CINner_introduction.html`.