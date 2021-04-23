# SPmethylseq <img src="https://raw.githubusercontent.com/systemPipeR/systemPipeR.github.io/main/static/images/SPR-Workflows.png" align="right" height="139" />

<!-- badges: start -->
[![R-CMD-check](https://github.com/systemPipeR/SPmethylseq/actions/workflows/R_CMD.yml/badge.svg)](https://github.com/systemPipeR/SPmethylseq/actions/workflows/R_CMD.yml)
[![Lifecycle: experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://www.tidyverse.org/lifecycle/#experimental)
[![Project Status: Active – The project has reached a stable, usable state and is being actively developed.](https://www.repostatus.org/badges/latest/active.svg)](https://www.repostatus.org/#active)
<!-- badges: end -->

### :construction: Under Development!

> This pipeline is currently under development and does not have a stable release yet.

### Installation

To install the package, please use the _`BiocManager::install`_ command:
```
if (!requireNamespace("BiocManager", quietly=TRUE))
    install.packages("BiocManager")
BiocManager::install("systemPipeR/SPmethylseq", build_vignettes=TRUE, dependencies=TRUE)
```
To obtain the *systemPipeR* and *systemPipeRdata*, please run as follow:
```
if (!requireNamespace("BiocManager", quietly=TRUE))
    install.packages("BiocManager")
BiocManager::install("systemPipeR")
BiocManager::install("systemPipeRdata")
```

