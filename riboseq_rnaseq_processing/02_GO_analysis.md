---
title: "riboseq project: GO analysis using customized deltaTE"
author: "Daeho Joe, Seunghee Jung"
---

## Aim:

* Demonstrate the processing using ([clusterProfiler](https://github.com/YuLab-SMU/clusterProfiler), [enrichR](https://maayanlab.cloud/Enrichr/)).


## GO Analysis  
We applied GO analysis using both clusterProfiler and enrichR on every gene that foldchange exceeds 1 or less than 1. We found useful terms in category of biological process and visualized them using dotplot and chorddiagram.  
#### Prerequisites
```r
require("BiocManager")
require("remotes")
require("ggplot2")
require("ggrepel")
require("dplyr")
require("tibble")
require("grid")
require("gridBase")
require("magrittr")
require("viridis")
require("clusterProfiler")
require("org.Hs.eg.db")
require("circlize")
require("ComplexHeatmap")
require("enrichplot"
require("scales")
require("ggsci")
```
#### For dotplot, we used follosing code  
```r
# enrichR result
GO_combined <- read_tsv("~/Documents/riboseq/enrichR/results.txt")
colnames(GO_combined) <- c("Term", "GeneRatio", "GeneCount", "padj", "cluster")

GO_combined$Term <- factor(x = GO_combined$Term,
                           levels = c("Regulation Of Transcription By RNA Polymerase II (GO:0006357)",
                                      "Regulation Of DNA-templated Transcription (GO:0006355)",
                                      "Positive Regulation Of Transcription By RNA Polymerase II (GO:0045944)",
                                      "Positive Regulation Of DNA-templated Transcription (GO:0045893)",
                                      "Positive Regulation Of Nucleic Acid-Templated Transcription (GO:1903508)",
                                      "mRNA Splicing, Via Spliceosome (GO:0000398)",
                                      "mRNA Processing (GO:0006397)",
                                      "Translation (GO:0006412)",
                                      "Gene Expression (GO:0010467)",
                                      "Ribosome Biogenesis (GO:0042254)",
                                      "Mitotic Sister Chromatid Segregation (GO:0000070)",
                                      "RNA Processing (GO:0006396)",
                                      "Mitotic Spindle Organization (GO:0007052)",
                                      "Aerobic Electron Transport Chain (GO:0019646)",
                                      "Mitochondrial ATP Synthesis Coupled Electron Transport (GO:0042775)",
                                      "Positive Regulation Of Translation (GO:0045727)",
                                      "Positive Regulation Of Macromolecule Biosynthetic Process (GO:0010557)",
                                      "Neutral Amino Acid Transport (GO:0015804)"),
                           order = TRUE)

GO_combined$cluster <- factor(x = GO_combined$cluster,
                              levels = c("Exc_1Down",
                                         "Exc_1Up",
                                         "Buf_1Down",
                                         "Buf_1Up",
                                         "Int_1Up"),
                              order = TRUE)

ggplot(GO_combined,
            aes_string(x = "cluster",
                       y = "Term")) +
  geom_point(aes(size = as.numeric(GeneCount),
                 color = as.numeric(-log10(padj)))) +
  scale_color_gradientn(name = "-log10(padj)",
                        colours = c("#8AB0A2", "#021B4E", "#021B4E"),
                        #limits= c(1E-1, 1E-4), 
                        breaks = c(5, 15, 25, 35)) +
  scale_size_area(name = c("GeneCount"),
                  breaks = c(20, 40, 60, 80),
                  max_size = 8) +
  theme_classic() +
  theme(plot.title = element_text(face = "bold"),
        axis.title.x = element_text(angle = 0, vjust = -1),
        axis.title.y = element_text(angle = 90, vjust = 1),
        axis.text.x = element_text(size = 10),
        strip.text = element_text(size = rel(1.1), face = "bold"),
        panel.grid.major = element_line(colour = "grey95"))
```

#### Chorddiagram was generatd by following code.
```r
df <- read.csv("TE_GO.csv", check.names = F, header = T)
df <- df[order(df$FoldChange, decreasing = T),]
group <- c(unique(df$categoryID),
           unique(df$Gene))

# Preparing legends
df_lg <- df %>%
  arrange(desc(FoldChange)) %>%
  .[, c("Gene", "FoldChange")] %>%
  unique()

values <- sort(df_lg$FoldChange)
normalized_values <- (values - min(values)) / (max(values) - min(values))
color_palette <- colorRampPalette(c("#EEE9CE", "#C39B6B"))(100)
matched_colors <- color_palette[round(normalized_values * 99) + 1]

# Setting colors
mycolors <- c(colorRampPalette(rev(c("#B64F4A", "#C39B6B", "#CBC7B1", "#8AB0A2", "#326482")))(length(unique(df$categoryID))),
              c(rev(matched_colors)))

names(mycolors) <- group

## plotting
circos.par(## edit  canvas size 
  canvas.ylim=c(-1.1, 1.1),
  canvas.xlim=c(-0.5, 1.1),
  ## adjust bottom and top margin
  ## track.margin = c(0.01, 0.1)
  track.margin = c(0.01, 0.01),
  track.height = 0.05,
  start.degree = 90)


cd <- chordDiagram(df[1:2],
                   order = group,
                   grid.col = mycolors,
                   annotationTrack = "grid", # plot only grid
                   preAllocateTracks = 2, # adjust grid
                   annotationTrackHeight = c(0.03, 0.01), 
                   h.ratio = 0.8, # adjust height of all links
                   transparency = 0.05 # adjust transparency
                   
)

# gene label
circos.track(track.index = 2, panel.fun = function(x, y) {
  circos.text(CELL_META$xcenter,
              CELL_META$ylim[1],
              CELL_META$sector.index[CELL_META$sector.index %in% df$Gene],
              facing = "reverse.clockwise",
              niceFacing = F,
              adj = c(1, 0.5),
              cex = 0.6)
}, bg.border = NA)



# add legend
min.value <- min(df$FoldChange, na.rm = TRUE)
max.value <- max(df$FoldChange, na.rm = TRUE)

if(max.value > 0 & min.value >= 0){
  col_fun = colorRamp2(c(max.value, min.value), c("#C39B6B", "#EEE9CE"))
} else if(max.value > 0 & min.value <= 0){
  col_fun = colorRamp2(c(max.value, 0, min.value), c("#807676", "white", "navy"))
} else if (max.value <= 0 & min.value <= 0) {
  col_fun = colorRamp2(c(max.value, min.value), c("white", "navy"))
}

legend_GO <- Legend(title = "GO term", labels = unique(df$categoryID),
                    type = "point", background = unique(cd$col),
                    labels_gp = gpar(font = 1))

legend_FC <- Legend(title = colnames(df[4]),
                    col_fun = col_fun,
                    direction = "horizontal",
                    title_gp = gpar(fontface = "plain",
                                    fontsize = 9.5),
                    # grid_with = unit(1, "cm"),
                    legend_width = unit(5, "cm"))

h = dev.size()[2]
legend_list = packLegend(legend_FC, legend_GO)
x_pos = unit(4*h, "cm")
y_pos = unit(0.2*h, "cm")

draw(legend_list, x = x_pos, y = y_pos, just = c("right", "bottom"))

circos.clear()
```
