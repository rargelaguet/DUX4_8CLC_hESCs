---
"HNES1_DUX4-OE_scATAC"
---

```{r}
suppressPackageStartupMessages(library(biovizBase))
suppressPackageStartupMessages(library(Signac))
suppressPackageStartupMessages(library(Seurat))
suppressPackageStartupMessages(library(EnsDb.Hsapiens.v86))
suppressPackageStartupMessages(library(BSgenome.Hsapiens.UCSC.hg38))
suppressPackageStartupMessages(library(rhdf5))
suppressPackageStartupMessages(library(hdf5r))
suppressPackageStartupMessages(library(Rcpp))
suppressPackageStartupMessages(library(SeuratDisk))
suppressPackageStartupMessages(library(ggplot2))
# suppressPackageStartupMessages(library(dplyr))
# suppressPackageStartupMessages(library(tidyverse))
# suppressPackageStartupMessages(library(Matrix))
# suppressPackageStartupMessages(library(spatstat))
# suppressPackageStartupMessages(library(patchwork))
# suppressPackageStartupMessages(library(cellranger))
# suppressPackageStartupMessages(library(ggrepel))
set.seed(13)
```

```{r message=FALSE}
# stuff
ZGA_genes <- c("CCNA1","DUXA","KDM4E","KHDC1L","LEUTX","MBD3L2","MBD3L3","PRAMEF1","PRAMEF12",
               "PRAMEF11","RFPL2","RFPL4A","RFPL4B","SLC34A2","TRIM43","TRIM43B","TRIM49","TRIM49B","ZNF296","ZSCAN4")

theme_set(theme_bw(base_size = 10))
theme(plot.title = element_text (hjust=0.5))
```

```{r}
# load the RNA and ATAC data
DUX4_counts <- Read10X_h5("/Users/stowersj/OneDrive - BABRAHAM/DataR/scSeq_naive_hESCs_Stowers_2021/all_HNES1_DUX4-OE/HNES1_DUX4_overexpression_L001_filtered_feature_bc_matrix.h5")
DUX4_fragpath <- "/Users/stowersj/OneDrive - BABRAHAM/DataR/scSeq_naive_hESCs_Stowers_2021/all_HNES1_DUX4-OE/HNES1_DUX4_overexpression_L001_atac_fragments.tsv.gz"

```

```{r message=FALSE}
# get gene annotations for hg38
annotation <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
seqlevelsStyle(annotation) <- "UCSC"
genome(annotation) <- "hg38"
```

```{r}
# create a Seurat object containing the RNA adata
DUX4 <- CreateSeuratObject(counts = DUX4_counts$`Gene Expression`,assay = "RNA")
```

```{r}
# create ATAC assay and add it to the object [2143 cells]
DUX4[["ATAC"]] <- CreateChromatinAssay(counts = DUX4_counts$Peaks,sep = c(":", "-"),fragments = DUX4_fragpath,annotation = annotation)
```

```{r}
# Quality control
DefaultAssay(DUX4) <- "ATAC"

DUX4 <- NucleosomeSignal(DUX4) #ratio mono-nucleosome fragments (147-294bp) to nucleosome free (<147bp)
DUX4 <- TSSEnrichment(DUX4)

VlnPlot(object = DUX4,features = c("nCount_RNA", "nCount_ATAC", "TSS.enrichment", "nucleosome_signal"),ncol = 4,pt.size = 0.1)

# filter out low quality cells
DUX4 <- subset(x = DUX4, subset = nCount_ATAC < 500000 & nCount_RNA < 150000 & nCount_ATAC > 2500 & nCount_RNA > 1000 
               & nucleosome_signal < 2 & TSS.enrichment > 1)

VlnPlot(object = DUX4,features = c("nCount_RNA", "nCount_ATAC", "TSS.enrichment", "nucleosome_signal"),ncol = 4,pt.size = 0.1)
```

```{r}
# call peaks using MACS2 [15min]
peaks <- CallPeaks(DUX4, macs2.path = "/Users/stowersj/opt/anaconda3/envs/macs2/bin/macs2")

# remove peaks on nonstandard chromosomes and in genomic blacklist regions
peaks <- keepStandardChromosomes(peaks, pruning.mode = "coarse")
peaks <- subsetByOverlaps(x = peaks, ranges = blacklist_hg38_unified, invert = TRUE)

# quantify counts in each peak [20min]
macs2_counts <- FeatureMatrix(fragments = Fragments(DUX4),features = peaks,cells = colnames(DUX4))

# create a new assay using the MACS2 peak set and add it to the Seurat object
DUX4[["peaks"]] <- CreateChromatinAssay(counts = macs2_counts,fragments = DUX4_fragpath,annotation = annotation)
```

```{r}
# Gene expression data processing
# DefaultAssay(DUX4) <- "SCT"
DUX4 <- SCTransform(DUX4, variable.features.n = 100) #replaces NormalizeData(), ScaleData(), FindVariableFeatures()
DUX4 <- RunPCA(DUX4)

head(VariableFeatures(DUX4), 100)
# LabelPoints(plot = VariableFeaturePlot(DUX4), points = top50, repel = TRUE, xnudge = 0, ynudge = 0)
```

```{r}
# DNA accessibility data processing
DefaultAssay(DUX4) <- "peaks"
DUX4 <- FindTopFeatures(DUX4, min.cutoff = 5)
DUX4 <- RunTFIDF(DUX4)
DUX4 <- RunSVD(DUX4)
```

```{r}
# build a joint neighbor graph using both assays
DefaultAssay(DUX4) <- "peaks"
# DefaultAssay(DUX4) <- "SCT"

DUX4 <- FindMultiModalNeighbors(object = DUX4,
                                reduction.list = list("pca", "lsi"), 
                                dims.list = list(1:30, 2:30),
                                modality.weight.name = "RNA.weight",verbose = TRUE)
```

```{r}
DUX4 <- RunUMAP(object = DUX4,nn.name = "weighted.nn",assay = "RNA",verbose = TRUE)
p1 <- DimPlot(DUX4, label = TRUE, repel = TRUE, reduction = "umap") + NoLegend()
DefaultAssay(DUX4) <- "SCT"
p2 <- FeaturePlot(DUX4, features = c("TPRX1"), reduction = "umap")
p3 <- FeaturePlot(DUX4, features = c("ZSCAN4"), reduction = "umap")

p1+p2+p3
```


```{r}
DUX4 <- FindNeighbors(object = DUX4)
DefaultAssay(DUX4) <- "SCT"
DUX4 <- FindClusters(object = DUX4, resolution = 0.7)

DUX4 <- RunUMAP(DUX4, dims = 1:30, n.neighbors = 30L)
DimPlot(DUX4, reduction = "umap", group.by ="seurat_clusters")

VlnPlot(DUX4, features = c("DPPA3", "KLF17", "SUSD2"))
VlnPlot(DUX4, features = c("ZSCAN4", "LEUTX", "DUXA"))
DotPlot(DUX4, features = ZGA_genes,group.by ="seurat_clusters" ) + theme(axis.text.x = element_text(angle = 30, hjust = 1), axis.title.x = element_blank()) + guides(fill=FALSE)

VlnPlot(object = DUX4,features = c("nCount_peaks", "nCount_ATAC", "TSS.enrichment", "nucleosome_signal"),ncol = 4,pt.size = 0.1)

FeaturePlot(DUX4, features = c("TPRX1"), reduction = "umap")
FeaturePlot(DUX4, features = c("ZSCAN4"), reduction = "umap")
```

```{r}
# compute GC content for each peak
DefaultAssay(DUX4) <- "peaks"
DUX4 <- RegionStats(DUX4, genome = BSgenome.Hsapiens.UCSC.hg38)

# link peaks to genes
DUX4 <- LinkPeaks(object = DUX4,peak.assay = "peaks", expression.assay = "SCT",genes.use = c("SOX2", "POU5F1"))

# idents.plot <- c("B naive", "B intermediate", "B memory", "CD14 Mono", "CD16 Mono", "CD8 TEM", "CD8 Naive")

p1 <- CoveragePlot(object = DUX4,region = "SOX2",features = "SOX2",expression.assay = "SCT",
                   # idents = idents.plot,
                   extend.upstream = 500, extend.downstream = 10000)

p2 <- CoveragePlot(object = DUX4,region = "POU5F1",features = "POU5F1",expression.assay = "SCT", 
                   # idents = idents.plot,
                   extend.upstream = 500, extend.downstream = 10000)

patchwork::wrap_plots(p1, p2, ncol = 1)
```

