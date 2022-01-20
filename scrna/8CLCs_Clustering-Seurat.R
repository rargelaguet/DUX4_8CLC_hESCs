library(Seurat)
library(ggplot2)
library(ggrepel)

theme_set(theme_bw(base_size = 10))
theme(plot.title = element_text (hjust=0.5))

ZGA_genes <- c("CCNA1","DUXA","KDM4E","KHDC1L","LEUTX","MBD3L2","MBD3L3","PRAMEF1","PRAMEF12",
               "PRAMEF11","RFPL2","RFPL4A","RFPL4B","SLC34A2","TRIM43","TRIM43B","TRIM49","TRIM49B","ZNF296","ZSCAN4")

##### ---- LOAD DATA ---- #####

##### ---- h0 + FILTER ---- #####

h0.data <- Seurat::Read10X(data.dir = "/Users/stowersj/OneDrive - BABRAHAM/DataR/scSeq_naive_hESCs_MariaR_2020/feature_bc_matrix/h0/")
h0 <- CreateSeuratObject(counts = h0.data, project = "h0", min.cells = 3, min.features = 200)
rm(h0.data)

h0[["percent.mt"]] <- PercentageFeatureSet(h0, pattern = "^MT-")
h0$percent.largest.gene <- apply(h0@assays$RNA@counts, 2, function(x)(100*max(x))/sum(x))

h0 <- subset(h0, subset = nFeature_RNA > 2500 & nFeature_RNA < 10500 & percent.mt < 10 & percent.largest.gene < 20)

##### ----  CLR-NORMALIZE ---- #### /.../

h0 <- NormalizeData(h0, normalization.method = "CLR")

##### ----  HIGHLY VARIABLE FEATURES ---- ####

h0 <- FindVariableFeatures(h0, selection.method = "vst", nfeatures = 500)

# top50 <- head(VariableFeatures(h0), 50)
# LabelPoints(plot = VariableFeaturePlot(h0), points = top50, repel = TRUE, xnudge = 0, ynudge = 0)

##### ----  SCALING ---- ####

all.genes <- rownames(h0)
h0 <- ScaleData(h0, features = all.genes)
rm(all.genes)

### ---- PCA ---- ##

h0 <- RunPCA(h0, features = VariableFeatures(object = h0))

### ---- TSNE ---- ##

# saved.seed <- 13
# set.seed(saved.seed)
# h0 <- RunTSNE(h0, dims = 1:30, seed.use=saved.seed, perplexity= 30)

##### ---- CLUSTER CELLS ---- ####

h0 <- FindNeighbors(h0, dims = 1:30)
h0 <- FindClusters(h0, resolution = 0.5)

##### ---- NON-LIN DEM REDUCTION UMAP  ---- ####

h0 <- RunUMAP(h0, dims = 1:30)

DimPlot(h0, reduction = "pca", group.by ="seurat_clusters")
DimPlot(h0, reduction = "umap", group.by ="seurat_clusters")

####### PLOTS #######

FeaturePlot(h0, features = c("ZSCAN4"), reduction ="umap")
FeaturePlot(h0, features = c("LEUTX"), reduction ="umap")

VlnPlot(h0, features = "ZSCAN4",group.by ="seurat_clusters", assay = "RNA" )
VlnPlot(h0, features = "LEUTX",group.by ="seurat_clusters", assay = "RNA" )
VlnPlot(h0, features = ZGA_genes)
VlnPlot(h0, features = "KLF17",group.by ="seurat_clusters", assay = "RNA" )

DotPlot(h0, features = ZGA_genes,group.by ="seurat_clusters" ) + theme(axis.text.x = element_text(angle = 30, hjust = 1), axis.title.x = element_blank()) + guides(fill=FALSE)

FeatureScatter(h0, feature1 = "DPPA3", feature2 = "CCNA1", group.by = 'seurat_clusters', slot='counts')
FeatureScatter(h0, feature1 = "ZSCAN4", feature2 = "LEUTX", group.by = 'seurat_clusters', slot='counts')
FeatureScatter(h0, feature1 = "MBD3L3", feature2 = "DUXA", group.by = 'seurat_clusters', slot='counts')
FeatureScatter(h0, feature1 = "H3.Y", feature2 = "TPRX1", group.by = 'seurat_clusters', slot='counts')

##### ---- FINDING MARKERS h0_ cluster 5_ 8CLCs ---- ####

h0_8CLC.markers.roc <- FindMarkers(h0, ident.1 = "5", test.use="roc", only.pos = TRUE, logfc.threshold = 0.25)
h0_8CLC.markers.wilcox <- FindMarkers(h0, ident.1 = "5",test.use="wilcox", only.pos = TRUE, logfc.threshold = 0.25)

ggplot(h0_8CLC.markers.roc, aes(x=avg_log2FC, y=myAUC, size=avg_diff, 
                            # col=avg_diff,
                            label=rownames(h0_8CLC.markers.roc))) +
  geom_point(alpha=0.6, color="black") +
  geom_text_repel(size=3, hjust = 0.1, nudge_x = 0.15) +
  labs(x="log2FC", y="AUC") +
  theme(plot.title = element_text(hjust=0.5)) +
  scale_size(range=c(0.1,8))

h0_8CLC.markers.fused <- merge(h0_8CLC.markers.roc, h0_8CLC.markers.wilcox, by=0, all=TRUE) 
ggplot(h0_8CLC.markers.fused, aes(x=avg_log2FC.x, y=(myAUC), size=-log2(p_val), 
                                  # col=avg_diff,
                                  label=Row.names)) +
  geom_point(alpha=0.6, color="black") +
  geom_text_repel(size=3, hjust = 0.2, nudge_x = 0.2) +
  labs(x="avg_log2FC", y="AUC") +
  theme(plot.title = element_text(hjust=0.5)) +
  scale_size(range=c(0.1,8))

