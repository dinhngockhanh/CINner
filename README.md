#   Simulating chromosomal instability during tumorigenesis

CINner [1] is a tool to simulate cancer evolution as a selective model driven by mutations and Copy Number Aberrations (CNAs).
It is designed to reproduce data from single-cell DNA sequencing such as Direct Library Preparation+ (DLP+).

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

For more detailed examples and usage, please refer to `vignettes/CINner_introduction.html`
and other files in `vignettes`, which recreate results in [1].

##  References

1.  Dinh KN, Vázquez-García I, Chan A, Malhotra R, Weiner A, McPherson AW, Tavaré S. CINner: modeling and simulation of chromosomal instability in cancer at single-cell resolution. bioRxiv. 2024 Apr 3.