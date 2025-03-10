---
title: "riboseq_project: MS data processing"
author: "Seunghee Jung, Daeho Joe"
---

## Data collection  
Database searching of all raw data files was performed in Proteome Discoverer 2.5 software(Thermo Fisher Scientific). SEQUEST-HT was used for database searching against Swissprot-Human database. Database searching against the corresponding reversed database was also performed to evaluate the false discovery rate (FDR) of peptide identification. The database searching parameters included precursor ion mass tolerance 10 ppm, fragment ion mass tolerance 0.02 Da, fixed modification for carbamidomethyl cysteine (+57.021 Da / C) and variable modifications for methionine oxidation (+15.995 Da / M) and phosphorylation (+79.966 Da / S, T, Y). We obtained FDR of less than 1% on the peptide level and filtered with the high peptide confidence.

## GO analysis  
#### Installation prerequisites
```r
## list of packages required by this analysis 
pkg <- c("BiocManager", "ggplot2", "dplyr", "tibble", "grid", "gridBase",
         "magrittr", "viridis")

## Check if the packages were installed if not install 
pkg.new<-pkg[!(pkg%in%installed.packages())] 
if (length(pkg.new)) { 
  install.packages(pkg.new) 
} 

## load packages
sapply(pkg, require, character.only = TRUE) 

## list of bioconductor packages required by this analysis 
Biopkg<-c("clusterProfiler", "org.Hs.eg.db", "circlize", "ComplexHeatmap",
          "enrichplot", "scales", "ggsci") 

## check if the packages were installed if not install 
Biopkg.new<-Biopkg[!(Biopkg%in%installed.packages())] 
if (length(Biopkg.new)) { 
  if (!requireNamespace("BiocManager", quietly = TRUE)) 
    install.packages("BiocManager") 
  BiocManager::install(Biopkg.new, type="source") 
} 

## load bioconductor packages 
sapply(Biopkg, require, character.only = TRUE)
```  

#### UniProt ID was converted into gene symbol using OrgDb. Following code is example of its usage.  
```r
enRich <- bitr(names(gene_List), fromType = "UNIPROT", toType = c("SYMBOL"),
               OrgDb = "org.Hs.eg.db")
```

#### GO was retrieved using clusterProfiler with following code, using converted gene symbols.
```r
ego <- enrichGO(gene          = enRich$SYMBOL,
                OrgDb         = org.Hs.eg.db,
                ont           = "BP",
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.01,
                qvalueCutoff  = 0.05,
                keyType       = "SYMBOL")

ego_df <- data.frame(ego)
```

#### Dotplot was generated by following code.  
```r
category_selected <- c("ribonucleoprotein complex biogenesis",
                       "ribonucleoprotein complex subunit organization",
                       "ribosome biogenesis",
                       "cytoplasmic translation",
                       "ribonucleoprotein complex assembly",
                       "ribosomal small subunit biogenesis",
                       "ribosome assembly",
                       "viral process",
                       "translational initiation",
                       "mitochondrial translation",
                       "cytoplasmic translational initiation",
                       "formation of cytoplasmic translation initiation complex")
padjust_df <- ego_df[order(ego_df$p.adjust),]
sorted_categories <- padjust_df$Description[padjust_df$Description %in% category_selected]
category_selected <- factor(category_selected,
                            levels = sorted_categories,
                            order = T)

d1 <- dotplot(ego,
              showCategory = category_selected)
d_data <- d1$data

min.value <- as.numeric(sprintf("%.1e", min(d1$data$p.adjust, na.rm = TRUE)))
max.value <- as.numeric(sprintf("%.1e", max(d1$data$p.adjust, na.rm = TRUE)))

max.ratio <- as.numeric(max(d1$data$GeneRatio, na.rm = TRUE))

d1 +
  theme_classic() +
  theme(plot.title = element_text(face = "bold"),
        axis.title.x = element_text(angle = 0, vjust = -1),
        axis.title.y = element_text(angle = 90, vjust = 1),
        axis.text.x = element_text(size = 10),
        strip.text = element_text(size = rel(1.1), face = "bold"),
        panel.grid.major = element_line(colour = "grey95")) +
  scale_x_continuous(limits = c(0, max.ratio+0.005),
                     breaks = c(0.05, 0.10),
                     labels = c("0.05", "0.10")) +
  scale_y_discrete(
    limits = rev(levels(category_selected)),
    label = stringr::str_wrap(rev(levels(category_selected)), 31)) +
  scale_color_gradientn(name = "p.adjust",
                        colours = c("#021B4E", "#CBC7B1"),
                        limits= c(min.value, max.value),
                        breaks=c(min.value, c(min.value, max.value)/2, max.value)
  )
```

#### Chorddiagram was generated by following code.  
```r
#### prepare data
## select category
category_selected <- c("translational initiation",
                       "cytoplasmic translational initiation",
                       "formation of cytoplasmic translation initiation complex")
category_selected <- factor(category_selected, levels = category_selected, order = T)

## prepare data with selected category
heat_df2 <- heatplot(ego,
                     showCategory = category_selected)$data
heat_df2 <- na.omit(heat_df2)
heat_df2 <- heat_df2 %>%
  left_join(., res_df[,c("SYMBOL", "Sum.PEP.Score")], by = c("Gene" = "SYMBOL")) %>%
  .[order(.$Sum.PEP.Score, decreasing = T),] %>%
  na.omit()
colnames(heat_df2) <- c("categoryID", "Gene", "foldChange")

## generate order
group <- c(unique(heat_df2$categoryID),
           unique(heat_df2$Gene))

## prepare data for the legend
uni.heat_df2 <- heat_df2 %>%
  arrange(desc(foldChange)) %>%
  .[, c("Gene", "foldChange")] %>%
  unique()
values <- sort(uni.heat_df2$foldChange)
normalized_values <- (values - min(values)) / (max(values) - min(values))
color_palette <- colorRampPalette(c("white", "#EEE9CE", "#B64F4A"))(100)
matched_colors <- color_palette[round(normalized_values * 99) + 1]

## color a layer of the chorddiagram
mycolors <- c(colorRampPalette(rev(c("#B64F4A", "#CBC7B1", "#326482")))(length(unique(heat_df2$categoryID))),
              c(rev(matched_colors)))
names(mycolors) <- group




#### start plotting
circos.par(canvas.ylim=c(-1.1, 1.1),
           canvas.xlim=c(-0.5, 1.1),
           track.margin = c(0.01, 0.01),
           track.height = 0.05,
           start.degree = 90)

cd <- chordDiagram(heat_df2[1:2],
                   
                   order = group,
                   
                   grid.col = mycolors,
                   
                   ## plot only grid (no labels, no axis)
                   annotationTrack = "grid", 
                   preAllocateTracks = 2,
                   
                   ## adjust grid width and spacing
                   annotationTrackHeight = c(0.03, 0.01),
                   
                   ## adjust height of all links
                   h.ratio = 0.8,
                   
                   ## transparency
                   transparency = 0.05
                   
)

## present gene labels
circos.track(track.index = 2, panel.fun = function(x, y) {
  circos.text(CELL_META$xcenter,
              CELL_META$ylim[1],
              CELL_META$sector.index[CELL_META$sector.index %in% heat_df2$Gene],
              facing = "reverse.clockwise",
              niceFacing = F,
              adj = c(1, 0.5),
              cex = 0.6)
}, bg.border = NA)

## present highlighted gene labels
gene_names <- unique(heat_df2$Gene)
gene_colors <- character(length(gene_names))
 
for (i in 1:length(gene_names)) {
  if (startsWith(gene_names[i], "EIF")) {
    gene_colors[gene_names[i]] <- "red"
  } else {
    gene_colors[gene_names[i]] <- "grey50"
  }
}

## add legends
min.value <- min(uni.heat_df2$foldChange, na.rm = TRUE)
max.value <- max(uni.heat_df2$foldChange, na.rm = TRUE)
col_fun = colorRamp2(c(max.value, (max.value + min.value)/2, 0), c("#B64F4A", "#EEE9CE", "white"))

legend_GO <- Legend(title = "GO term", labels = unique(heat_df2$categoryID),
                    type = "point", background = unique(cd$col),
                    labels_gp = gpar(font = 1))
legend_FC <- Legend(title = "Sum.PEP.Score",
                    col_fun = col_fun,
                    direction = "horizontal",
                    title_gp = gpar(fontface = "plain",
                                    fontsize = 9.5),
                    legend_width = unit(5, "cm"))
h = dev.size()[2]
legend_list = packLegend(legend_FC, legend_GO)
draw(legend_list, x = x_pos, y = y_pos, just = c("right", "bottom"))

circos.clear()
```
