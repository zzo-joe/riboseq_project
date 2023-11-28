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

## Preprocessing
We trimmed our raw **riboseq data** fastq files following manufacturer’s guide(https://immaginabiotech.com/product/laceseq-12-rxns, LACEseq bioIT guidelines). Reads were processed using Cutadapt(version 2.8), UMI-tools(version 1.1.4), and were trimmed from the 3’ end of each read and only reads longer than X+9 were retrieved. The processed sequences were first mapped to the contaminants using bowtie2 with parameters “—un”, which is composed of PhiX genom([GenBank accession J02482.1](https://www.ncbi.nlm.nih.gov/nuccore/J02482.1)), predicted tRNA genes ([GtRNAdb](http://gtrnadb.ucsc.edu/genomes/eukaryota/Hsapi38/Hsapi38-seq.html)), lncRNA transcripts ([Gencode](https://www.gencodegenes.org/human/)), rRNA repeat units ([GenBank accession U13369.1](https://www.ncbi.nlm.nih.gov/nuccore/U13369.1)) and all sequences for ribosomal RNA of human (retrieved from Rfam 11.0 of the Wellcome Trust Sanger Institute). Reads not mapped to the genes above were used for following analysis.  

Trimmed riboseq data and raw RNA-seq data were mapped to human genome GRCh38 using STAR (v2.7.10b) following [deltaTE Support Protocol](https://currentprotocols.onlinelibrary.wiley.com/doi/full/10.1002/cpmb.108#cpmb108-prot-0003). Aligned reads were sorted by samtools(v1.10) followed by quantification using featureCounts(v2.0.0).

## Read in countdata

deltaTE paper: [Chotani et al., 2019](https://doi.org/10.1002/cpmb.108)

Explaination credit: [SGDDNB github](https://github.com/SGDDNB/translational_regulation)

Calculating differential translation genes (DTGs) requires the count matrices from Ribo-seq and RNA-seq. These should be the raw counts obtained from feature counts or any other tool for counting reads, they should not be normalized or batch corrected. It also requires a sample information file which should be in the same order as samples in the count matrices. It should include information on sequencing type, treatment, batch or any other covariate you need to model.

We obtained countdata table from featureCounts output by following code:

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


DESeq2 object with batch for Ribo-seq

```r
ind = which(coldata$SeqType == "RIBO")
coldata_ribo = coldata[ind,]
ddsMat_ribo <- DESeq(ddsMat_ribo)
res_ribo <- results(ddsMat_ribo, contrast=c("Condition","KD","CTRL"))
res_ribo <- lfcShrink(ddsMat_ribo, coef=2,res=res_ribo,type="apeglm")
```


DESeq2 object with batch for RNA-seq

```r
ind = which(coldata$SeqType == "RNA")
coldata_rna = coldata[ind,]
ddsMat_rna <- DESeq(ddsMat_rna)
res_rna <- results(ddsMat_rna, contrast=c("Condition","KD","CTRL"))
res_rna <- lfcShrink(ddsMat_rna, coef=2,type="apeglm",res=res_rna)
```

For further analysis, classify genes as follow

```r
system("mkdir gene_lists")
forwarded = rownames(res)[which(res$padj > 0.05 & res_ribo$padj < 0.05 & res_rna$padj < 0.05)]
write.table(forwarded,"gene_lists/forwarded.txt",quote=F,sep="\t",col.names = F,row.names = F)

exclusive = rownames(res)[which(res$padj < 0.05 & res_ribo$padj < 0.05 & res_rna$padj > 0.05)]
write.table(exclusive,"gene_lists/exclusive.txt",quote=F,sep="\t",col.names = F,row.names = F)

both = which(res$padj < 0.05 & res_ribo$padj < 0.05 & res_rna$padj < 0.05)
neither = which(res$padj >= 0.05 & res_ribo$padj >= 0.05 & res_rna$padj >= 0.05)

intensified = rownames(res)[both[which(res[both,2]*res_rna[both,2] > 0)]]
write.table(intensified,"gene_lists/intensified.txt",quote=F,sep="\t",col.names = F,row.names = F)

buffered = rownames(res)[both[which(res[both,2]*res_rna[both,2] < 0)]]
buffered = c(rownames(res)[which(res$padj < 0.05 & res_ribo$padj > 0.05 & res_rna$padj < 0.05)],buffered)
write.table(buffered,"gene_lists/buffered.txt",quote=F,sep="\t",col.names = F,row.names = F)
```

Generated gene class will be plotted with following code

```r
plot(y=res_ribo[,2],x=res_rna[,2],
     xlab = "RNA-seq log2 fold change",
     ylab = "Ribo-seq log2 fold change",
     asp=1,pch=16,
     col=rgb(128/255,128/255,128/255,0.1),
     ylim=c(-6,6),
     xlim=c(-6,6),
     cex=0.4)

abline(a=0,b=1,col="gray")
abline(h=0,v=0,col="gray")

points(y=res_ribo[intensified,2],x=res_rna[intensified,2],pch=16,col="#CBC7B1")
points(y=res_ribo[buffered,2],x=res_rna[buffered,2],pch=16,col="#C39B6B")
points(y=res_ribo[forwarded,2],x=res_rna[forwarded,2],pch=16,col="#8AB0A2")
points(y=res_ribo[exclusive,2],x=res_rna[exclusive,2],pch=16,col="#021B4E")


legend("topleft",
       c("Forwarded","Exclusive","Buffered", "Intensified"),
       fill=c("#8AB0A2", "#021B4E", "#C39B6B", "#CBC7B1"),
       cex=1, border = NA, bty="n")
```

To find out gene profiles, we plotted with following code

```r
par(mfrow=c(2,2))

goi = forwarded[1]
y_u = max(res[goi,2],res_ribo[goi,2], res_rna[goi,2],0)
y_l = min(res[goi,2],res_ribo[goi,2], res_rna[goi,2],0)
plot(c(1,2),c(0,res[goi,2]),type="l",col="#AA4A44",xaxt="n",
     xlab="Conditions",ylim=c(y_l,y_u),ylab="Log2 Fold Change", 
     main="Forwarded gene")
lines(c(1,2),c(0,res_ribo[goi,2]),col="#BFBA9F")
lines(c(1,2),c(0,res_rna[goi,2]),col="#011640")
axis(1,at=c(1,2),labels=c("CTRL","KD"),las=1)
legend("bottomleft",c("mRNA","RPF","TE"), fill=c("#011640","#BFBA9F","#AA4A44"),
       cex=1, border = NA, bty="n")

goi = exclusive[1]
y_u = max(res[goi,2],res_ribo[goi,2], res_rna[goi,2],0)
y_l = min(res[goi,2],res_ribo[goi,2], res_rna[goi,2],0)
plot(c(1,2),c(0,res[goi,2]),type="l",col="#AA4A44",xaxt="n",
     xlab="Conditions",ylim=c(y_l,y_u),ylab="Log2 Fold Change", 
     main="Exclusive gene")
lines(c(1,2),c(0,res_ribo[goi,2]),col="#BFBA9F")
lines(c(1,2),c(0,res_rna[goi,2]),col="#011640")
axis(1,at=c(1,2),labels=c("CTRL","KD"),las=1)
legend("bottomleft",c("mRNA","RPF","TE"), fill=c("#011640","#BFBA9F","#AA4A44"),
       cex=1, border = NA, bty="n")

goi = buffered[1]
y_u = max(res[goi,2],res_ribo[goi,2], res_rna[goi,2],0)
y_l = min(res[goi,2],res_ribo[goi,2], res_rna[goi,2],0)
plot(c(1,2),c(0,res[goi,2]),type="l",col="#AA4A44",xaxt="n",
     xlab="Conditions",ylim=c(y_l,y_u),ylab="Log2 Fold Change", 
     main="Buffered gene")
lines(c(1,2),c(0,res_ribo[goi,2]),col="#BFBA9F")
lines(c(1,2),c(0,res_rna[goi,2]),col="#011640")
axis(1,at=c(1,2),labels=c("CTRL","KD"),las=1)
legend("bottomleft",c("mRNA","RPF","TE"), fill=c("#011640","#BFBA9F","#AA4A44"),
       cex=1, border = NA, bty="n")

goi = intensified[1]
y_u = max(res[goi,2],res_ribo[goi,2], res_rna[goi,2],0)
y_l = min(res[goi,2],res_ribo[goi,2], res_rna[goi,2],0)
plot(c(1,2),c(0,res[goi,2]),type="l",col="#AA4A44",xaxt="n",
     xlab="Conditions",ylim=c(y_l,y_u),ylab="Log2 Fold Change", 
     main="Intensified gene")
lines(c(1,2),c(0,res_ribo[goi,2]),col="#BFBA9F")
lines(c(1,2),c(0,res_rna[goi,2]),col="#011640")
axis(1,at=c(1,2),labels=c("CTRL","KD"),las=1)
legend("bottomleft",c("mRNA","RPF","TE"), fill=c("#011640","#BFBA9F","#AA4A44"),
       cex=1, border = NA, bty="n")

```

