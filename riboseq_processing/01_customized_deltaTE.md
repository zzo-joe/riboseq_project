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


## Read in countdata

deltaTE paper: [Chotani et al., 2019](https://doi.org/10.1002/cpmb.108)

Explaination credit: [SGDDNB github](https://github.com/SGDDNB/translational_regulation)

Calculating differential translation genes (DTGs) requires the count matrices from Ribo-seq and RNA-seq. These should be the raw counts obtained from feature counts or any other tool for counting reads, they should not be normalized or batch corrected. It also requires a sample information file which should be in the same order as samples in the count matrices. It should include information on sequencing type, treatment, batch or any other covariate you need to model.


Countdata should look like table below

Ribo-seq count matrix (RPFs): 

 | Gene ID | Sample 1 | Sample 2 | Sample 3 | Sample 4 |
 | --------|----------|----------|----------|----------|
 | Gene 1  | 1290     | 130      | 2	   | 1000     |
 | Gene 2  | 2	     | 10	| 5	   | 1	      |
 | ..	  |	     | 		|	   |	      |	
 | ..	  |	     | 		|	   |	      |	
 | ..	  |	     | 		|	   |	      |	
 | ..	  |	     | 		|	   |	      |	
 | Gene Z  | 200	     | 140	| 15	   | 11	      |


RNA-seq count matrix (mRNA counts): 

 | Gene ID | Sample 5 | Sample 6 | Sample 7 | Sample 8 |
 | --------|----------|----------|----------|----------|
 | Gene 1  | 1290     | 130      | 2	   | 1000     |
 | Gene 2  | 2	     | 10	| 5	   | 1	      |
 | ..	  |	     | 		|	   |	      |	
 | ..	  |	     | 		|	   |	      |	
 | ..	  |	     | 		|	   |	      |	
 | ..	  |	     | 		|	   |	      |	
 | Gene Z  | 200	     | 140	| 15	   | 11	      |

 We obtained countdata table by following code:

```r
## Import riboseq featurecounts results
countdata <- read.table(file.path("~/Documents/riboseq/2023_Pygo2_KD_Mitosis_riboseq_featurecounts.txt"),
                        header = TRUE,
                        check.names = TRUE,
                        row.names = 1)

# 3 by 3 comparison
sample <- c("RPF2a", "RPF2b","RPF2c", 
            "RPF35a", "RPF35b", "RPF35c")
siRNA <- c("RPF", "RPF", "RPF", 
           "RPF35", "RPF35","RPF35")
control <- c("control", "control", "control",
             "siRNA", "siRNA","siRNA")

colnames(countdata) <- sample

ribo <- countdata


## Import rnaseq featurecounts results
countdata <- read.table(file.path("~/Documents/riboseq/2023_S6K_Pygo2_eIF3D_RNAseq_Featurecounts.txt"),
                        header = TRUE,
                        check.names = TRUE,
                        row.names = 1)
# 3 by 3 comparison
sample <- c("RNA21", "RNA22","RNA23", 
            "RNA35a", "RNA35b", "RNA35c")
siRNA <- c("RNA", "RNA","RNA",
           "RNA35", "RNA35","RNA35")
control <- c("control", "control","control",
             "siRNA", "siRNA","siRNA")

colnames(countdata) <- sample

rna <- countdata

# merge ribo and rna count data
merge <- cbind(ribo,rna)
```

Next, we also need to prepare a also requires a sample information file which should follow the same sample order as count matrices. This file outlines the condition, sequencing type and batch for each sample. This script requires the same header names as shown below (case sensitive). 

Sample information file:

 | SampleID | Condition | SeqType |
 | --------|----------|----------|
 | RPF2a  | 1     | RIBO      |
 | RPF2b  | 1     | RIBO      |
 | RPF2c  | 1     | RIBO      |
 | RPF35a | 2     | RIBO      |
 | RPF35b | 2     | RIBO      |
 | RPF35c | 2     | RIBO      | 
 | RNA21  | 1     | RNA      |
 | RNA22  | 1     | RNA      |
 | RNA23  | 1     | RNA      |
 | RNA35a | 2     | RNA      |
 | RNA35b | 2     | RNA      |
 | RNA35c | 2     | RNA      |

 We obtained sample information table by following code:
 
```r
SampleID <- c("RPF2a", "RPF2b","RPF2c", 
              "RPF35a", "RPF35b", "RPF35c",
              "RNA21", "RNA22","RNA23",
              "RNA35a", "RNA35b", "RNA35c")

Condition <- c("CTRL", "CTRL", "CTRL", "KD", "KD", "KD",
               "CTRL", "CTRL", "CTRL", "KD", "KD", "KD")

SeqType <- c("RIBO", "RIBO", "RIBO", "RIBO", "RIBO", "RIBO",
             "RNA", "RNA", "RNA", "RNA", "RNA", "RNA")

sampleinfo <- as.data.frame(cbind(SampleID = SampleID,
                                  Condition = Condition,
                                  SeqType = SeqType))
```

Detecting differential translation regulation using DESeq2 by following code:

```r
ddsMat <- DESeqDataSetFromMatrix(countData = merge,
                                 colData = coldata, design =~ Condition + SeqType + Condition:SeqType)

ddsMat$SeqType = relevel(ddsMat$SeqType,"RNA")
ddsMat <- DESeq(ddsMat)
resultsNames(ddsMat)

##> resultsNames(ddsMat)
## [1] "Intercept"               "Condition_KD_vs_CTRL"   
## [3] "SeqType_RIBO_vs_RNA"     "ConditionKD.SeqTypeRIBO"

```


Choose the term you want to look at from resultsNames(ddsMat). 
ConditionKD.SeqTypeRIBO means Changes in Ribo-seq levels in KD vs CTRL accounting for changes in RNA-seq levels in KD vs CTRL.
This case, changes in riboseq in KD vs. CTRL accounting for changes in RNA seq levels

```r
res <- results(ddsMat, contrast=list("ConditionKD.SeqTypeRIBO"))
summary(res)

length(which(res$padj < 0.05))
write.table(rownames(res)[which(res$padj < 0.05)],
            "DTEGs.txt",
            quote=F, sep="\t", col.names = F, row.names = F)

res_df <- data.frame(res)
```




