---
title: "riboseq project: cutomized script from deltaTE"
author: "Daeho Joe"
---

## Aim:

* Demonstrate the processing using ([dealtaTE](https://github.com/SGDDNB/translational_regulation)).

* Adaptor sequences were removed and aligned to contanimants using Bowtie2 following manufacturer's protocol.
  See [LACEseq bioIT guidelines](https://immaginabiotech.com/storage/2022/10/24/378bdb72845722f9db650562abadf0981d57f5f3.pdf)

* Cleaned output were aligned and counted by procotol as described previously.
  See [deltaTE Support Protocol](https://currentprotocols.onlinelibrary.wiley.com/doi/full/10.1002/cpmb.108#cpmb108-prot-0003)

* 

## Installation prerequisites

Install the required CRAN R packages and load packages by following code:

```r
# CRAN packages
pkg <- c("BiocManager", "ggplot2", "tibble", "dplyr", "grid", "ggrepel", "stringr") 

# Check if the packages were installed if not install 
pkg.new<-pkg[!(pkg%in%installed.packages())] 
if (length(pkg.new)) { 
  install.packages(pkg.new) 
} 

# load packages 
sapply(pkg, require, character.only = TRUE) 
```

Install the required Bioconductor R packages and load packages by following code:

```r
# Bioconductor packages
Biopkg<-c("DESeq2", "scales", "ggsci", "stringr", "DEGreport", "tidyverse",
          "scater", "ReactomePA", "Homo.sapiens", "multtest", "org.Hs.eg.db",
          "AnnotationHub", "ensembldb", "clusterProfiler", "ComplexHeatmap") 

# Check if the packages were installed if not install 
Biopkg.new<-Biopkg[!(Biopkg%in%installed.packages())] 
if (length(Biopkg.new)) { 
  if (!requireNamespace("BiocManager", quietly = TRUE)) 
    install.packages("BiocManager") 
  BiocManager::install(Biopkg.new, type="source") 
} 

# load bioconductor packages 
sapply(Biopkg, require, character.only = TRUE)
```

