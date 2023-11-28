---
title: "riboseq project: quality control using riboWaltz"
author: "Daeho Joe"
---

## Aim:

* Demonstrate the processing using [riboWaltz](https://github.com/LabTranslationalArchitectomics/riboWaltz).

* Samples were processed using [IMMAGINA's proprietary RiboLace technology](https://immaginabiotech.com/) and processed raw files by following document. [LACEseq bioIT guidelines](https://immaginabiotech.com/storage/2022/10/24/378bdb72845722f9db650562abadf0981d57f5f3.pdf)

## Code we use  
#### Prerequisites
```r
require("riboWaltz")
require("BiocManager")
require("remotes")
require("ggplot2")
require("ggrepel")
require("Biostrings")
require("GenomicAlignments")
require("GenomicFeatures")
require("GenomicRanges")
require("IRanges")
```

#### Acquire input files
> * BAM files were made using STAR with option "--quantMode TranscriptomeSAM"  
> * Reference was retrieved from gencode, used [release 43](https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_43/) for analysis.  
```r
annotation_table <- create_annotation(gtfpath = "~/Documents/reference/gencode.v43.chr_patch_hapl_scaff.annotation.gtf")
# sample name should not include digits
reads_list <- bamtolist(bamfolder = "transcript_bam/",
                        annotation = annotation_table,
                        refseq_sep = " ")
```  

#### Data filtering and visualization
```r
filtered_list <- length_filter(data = filtered_list,
                               length_filter_mode = "custom",
                               length_range = 26:34)
psite_offset <- psite(filtered_list, flanking = 6, extremity = "auto")
reads_psite_list <- psite_info(filtered_list, psite_offset)

# codon coverage
codon_coverage_example <- codon_coverage(reads_psite_list, annotation_table)
# cds coverage
cds_coverage_example <- cds_coverage(reads_psite_list, annotation_table)

## average of read length of all sample
length_dist_rep <- rlength_distr(reads_psite_list,
                                 sample = list("Samp_avg" = c("RPF2a", "RPF2b", "RPF2c", "RPF35a", "RPF35b", "RPF35c")),
                                 cl = 99, multisamples = "average",
                                 colour = "gray70")
length_dist_rep[["plot"]]

frames <- frame_psite(reads_psite_list,
                      sample = sample,  
                      region = "all")
frames[["plot"]]

# Applied for loop to run all the samples in use
for (sp in c("RPF2a","RPF2b", "RPF2c", "RPF35a","RPF35b", "RPF35c")){
  
  ## Distribution of reads lengths per specified sample (A)
  length_dist <- rlength_distr(reads_psite_list, sample = sp)
  p <- length_dist[[paste0("plot_", sp)]]
  ggsave(paste0(sp, "_length_dist.png"), plot = p, width = 8, height = 8, dpi = 300)
  
  ## P-sites per region (B)
  ## computed output of the p-sites falling in the three annoted transcript regions
  ## a bar called "RNAs" displaying the expected read distribution from a random fragmentation of RNA
  region_psite <- region_psite(reads_psite_list,
                               annotation_table,
                               sample = sp)
  p <- region_psite[["plot"]]
  ggsave(paste0(sp, "_region_psite.png"), plot = p, width = 8, height = 8, dpi = 300)
  
  ## Trinucleotied periodicity (C)
  ## showing if, and to which extent, the identified P-sites results in codon periodicity on the CDS
  frames_stratified <- frame_psite_length(reads_psite_list,
                                          sample = sp,
                                          region = "all",
                                          cl = 90)
  
  p <- frames_stratified[["plot"]]
  ggsave(paste0(sp, "_trinucleotied_periodicity.png"), plot = p, width = 8, height = 6, dpi = 300)
  
  ## barplot style
  frames <- frame_psite(reads_psite_list,
                        sample = sp,  
                        region = "all")
  p <- frames[["plot"]]
  ggsave(paste0(sp, "_trinucleotied_periodicity_bar.png"), plot = p, width = 8, height = 8, dpi = 300)
  
  ## Metaheatmaps (D)
  ## consisted of four metaheatmaps
  ## displaying the abundance of the 5' and 3' extremity of reads mapping on 
  ## and around the start and stop codon of annoted CDS
  
  ends_heatmap <- rends_heat(reads_psite_list,
                             annotation_table,
                             sample = sp,
                             cl = 95, utr5l = 25, cdsl = 40, utr3l = 25)
  p <- ends_heatmap[["plot"]]
  ggsave(paste0(sp, "_metaheatmap.png"), plot = p, width = 8, height = 6, dpi = 300)
  
  metaprofile <- metaprofile_psite(reads_psite_list, 
                                   annotation_table, 
                                   sample = sp,
                                   utr5l = 20, cdsl = 40, utr3l = 20,                
                                   plot_title = paste0(sp, "_transcript"))
  p <- metaprofile[[paste0("plot_", sp)]]
  ggsave(paste0(sp, "_metaprofile_psite.png"), plot = p, width = 8, height = 6, dpi = 300)
  
  
  ## Codon usage (E)
  cu_barplot <- codon_usage_psite(reads_psite_list,
                                  annotation_table,
                                  sample = sp,
                                  fastapath = "~/Documents/reference/GRCh38.p13.genome.fa",
                                  fasta_genome = TRUE,
                                  gtfpath = "~/Documents/reference/gencode.v43.chr_patch_hapl_scaff.annotation.gtf",
                                  frequency_normalization = TRUE,
                                  refseq_sep = " ")
  
  p <- cu_barplot[["plot"]]
  ggsave(paste0(sp, "_codon_usage.png"), plot = p, width = 8, height = 8, dpi = 300)
  
}

```
