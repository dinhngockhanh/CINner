#   Simulation algorithm for chromosomal instability during cancer development

##  Installation

CINner [[1]](https://www.biorxiv.org/content/10.1101/2024.04.03.587939v1) is an algorithm to simulate cancer evolution
as a selective model driven by mutations and Copy Number Aberrations (CNAs).
It is designed to reproduce data from single-cell DNA sequencing
such as Direct Library Preparation+ (DLP+).

CINner requires some external libraries:

```{r}
BiocManager::install("ggtree")
BiocManager::install("HMMcopy")
devtools::install_github("shahcompbio/signals")
devtools::install_github("chrisamiller/fishplot")
```

The CINner library can then be installed with

```R
devtools::install_github("dinhngockhanh/CINner")
```

## Run CINner

For more detailed examples and usage, please refer to `vignettes/CINner_introduction.html`
and other files in `vignettes`, which recreate results in [1].

##  References

1.  Dinh KN, Vázquez-García I, Chan A, Malhotra R, Weiner A, McPherson AW, Tavaré S.
CINner: modeling and simulation of chromosomal instability in cancer at single-cell resolution.
bioRxiv. 2024 Apr 3.