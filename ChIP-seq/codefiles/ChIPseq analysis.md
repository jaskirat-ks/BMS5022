---
title: "ChIPseq analysis"
author: "Jaskirat"
date: "23/09/2021"
output: html_document
---

# References
```{r}
# https://guangchuangyu.github.io/2016/02/covplot-supports-grangeslist/
# https://guangchuangyu.github.io/2014/04/visualization-methods-in-chipseeker/
# http://bioconductor.riken.jp/packages/3.4/bioc/vignettes/ChIPseeker/inst/doc/ChIPseeker.html
# https://bioconductor.org/packages/release/bioc/vignettes/ChIPseeker/inst/doc/ChIPseeker.html
# https://github.com/hbctraining/Intro-to-ChIPseq/blob/master/schedule/2-day.md
```

```{r}
# Load libraries
library(ChIPseeker)
library(TxDb.Scerevisiae.UCSC.sacCer3.sgdGene)
library(clusterProfiler)
library(annotables)
library(GenomicRanges)
library(org.Sc.sgd.db)
library(ggplot2)
```

# ChIP profiling
```{r}
# read in all the bed files 
samplefiles <-list.files(
  "/Users/jaskiratkaursandhu/Documents/bioinformatics/BMS5022/Topic 2/chipseeker R/chr.bed",
  pattern= "*.bed",
  full.names=T) 
samplefiles <- as.list(samplefiles) 
names(samplefiles) <- c("Hog1IP_YEPD",
                        "Hog1IP_KCI",
                        "Sko1IP_YEPD",
                        "Sko1IP_KCI",
                        "Hot1IP_KCI")
print(samplefiles)

# read in the bed files as GRangesList 
peak=GenomicRanges::GRangesList(
  Hog1IP_KCI=readPeakFile(samplefiles[["Hog1IP_KCI"]]),
  Hog1IP_YEPD=readPeakFile(samplefiles[["Hog1IP_YEPD"]]),
  Sko1IP_KCI=readPeakFile(samplefiles[["Sko1IP_KCI"]]),
  Sko1IP_YEPD=readPeakFile(samplefiles[["Sko1IP_YEPD"]]),
  Hot1IP_KCI=readPeakFile(samplefiles[["Hot1IP_KCI"]]))
```

### ChIP peaks coverage plot
```{r, fig.height=5, fig.width=6}
# ChIP peaks coverage plot which shows the peak locations over the whole genome across the GRangesList object (peak). The default setting merges the coverage plots with different colours. 
p <- covplot(peak)
print(p)
```

```{r, fig.height=5, fig.width=6}
# Separate the coverage plots using facet_grid

col <- c(Hog1IP_YEPD= 'red', 
         Hog1IP_KCI= 'green', 
         Sko1IP_YEPD= 'pink',
         Sko1IP_KCI= 'blue',
         Hot1IP_KCI= 'orange')

p + 
  facet_grid(chr ~ .id) + 
  scale_color_manual(values=col) + 
  scale_fill_manual(values=col) +
  theme(plot.title = element_text(size = 20))+
  theme(axis.title.x = element_text(size = 15)) + 
  theme(legend.position = "none")

# download plot as image
ggsave('chip_peaks.png')
```

# ChIP peak data set comparison 

### Average Profile of ChIP peaks binding to TSS region
```{r, fig.height=4, fig.width=6}
# To calculate the profile of ChIP peaks binding to TSS regions, prepare the TSS regions, which are defined as the flanking sequence of the TSS sites. Then, align the peaks that are mapping to these regions, and generate the tagMatrix.
# Note: the tagMatrix is not restricted to TSS regions. The regions can be other types that defined by the user.

# TxDb object contains the transcript-related features of a particular genome. In this case, the yeast genome is used.
txdb <- TxDb.Scerevisiae.UCSC.sacCer3.sgdGene

# Prepare the promotor regions
promoter <- getPromoters(TxDb=txdb, upstream=3000, downstream=3000)

# Calculate the tag matrix
tagMatrixList <- lapply(as.list(samplefiles), getTagMatrix, windows=promoter)

# plotAvgProf() - plots the profile of peaks
# Confidence interval estimated by bootstrap method is also supported for characterizing ChIP binding profiles.
plotAvgProf(tagMatrixList, xlim=c(-3000, 3000), conf=0.95,resample=500, facet="row")

# download plot as image
ggsave('avg-plot.png')
```

### Heatmap of ChIP binding to TSS regions
```{r, fig.height=3, fig.width=6}
# Plot heatmap
tagHeatmap(tagMatrixList, xlim=c(-3000, 3000), color=NULL)

# download plot as image
ggsave('heatmap.png')
```

# Peak Annotation
```{r}
# Perform peak annotations across all files where the TSS (transcription start site) region is defined from -3kb to +3kb.

# annotatePeak() assigns peaks to genomic annotation in the annotation column of the output, which includes whether a peak is in the TSS, Exon, 5'UTR, 3'UTR, Intronic or Intergenic regions. Moreover, the distance from the peak (binding site) to the TSS of the nearest gene is calculated by annotatePeak() and reported in the output. 

# peakAnnoList is a named list of annotated peaks

peakAnnoList <- lapply(samplefiles, annotatePeak, TxDb=txdb, 
                       tssRegion=c(-3000, 3000), verbose=TRUE)
```

```{r}
# view the peak annotations for each file
peakAnnoList
```

# Visualisations of genomic region annotation
```{r}
# Genomic annotation by barplot
plotAnnoBar(peakAnnoList)

# download plot as image
ggsave('feature_dist.png')
```

```{r}
# This plot breaks down the location of each peak relative to the TSS. The x-axis gives the percentage of sites, whilst the colour represents the distance from the TSS. Basically, plotDistToTSS calculates the percentage of binding sites upstream and downstream from the TSS of the nearest genes, and visualize the distribution. 
plotDistToTSS(peakAnnoList,
              title="Distribution of transcription factor-binding loci relative to TSS") +
  theme(plot.title= element_text(hjust = 0.5))

# download plot as image
ggsave('tss.png')
```

```{r}
# Download the peak Annotations for each bed file as separate csv files

Map(function(x, y) write.csv(x@anno, y, row.names = FALSE), peakAnnoList, paste0(names(peakAnnoList), '.csv'))
```

# Overlap of peaks and annotated genes 
```{r, fig.height=3, fig.width=5.5}
genes= lapply(peakAnnoList, function(i) as.data.frame(i)$geneId)
vennplot(genes)

ggsave("vennplot.png")
```