suppressPackageStartupMessages(library(Seurat))
suppressPackageStartupMessages(library(argparse))

here::i_am("rna/clustering/seurat_clustering.R")

######################
## Define arguments ##
######################

p <- ArgumentParser(description='')
p$add_argument('--seurat',         type="character",                            help='Seurat file')
p$add_argument('--metadata',    type="character",                            help='metadata file')
p$add_argument('--nfeatures',    type="integer",    default=1500,             help='Number of highly variable features')
p$add_argument('--npcs',        type="integer",    default=15,               help='Number of PCs')
p$add_argument('--knn',         type="integer",    default=15,               help='Number of k nearest neighbours for the graph')
p$add_argument('--clustering_resolution',    type="double",    default=0.5,               help='Clustering resolution')
p$add_argument('--n_neighbors', type="integer",    default=15,               help='(UMAP) Number of neighbours')
p$add_argument('--min_dist',    type="double",     default=0.3,              help='(UMAP) Minimum distance')
p$add_argument('--outdir',      type="character",                                  help='Output directory')

args <- p$parse_args(commandArgs(TRUE))

#####################
## Define settings ##
#####################

source(here::here("settings.R"))
source(here::here("utils.R"))

## START TEST ##
# args$samples <- opts$samples
# args$seurat <- io$rna.seurat
# args$metadata <- file.path(io$basedir,"results/rna/qc/sample_metadata_after_qc.txt.gz")
# args$nfeatures <- 1500
# args$clustering_resolution <- 0.5
# args$npcs <- 15
# args$knn <- 15
# args$n_neighbors <- 15
# args$min_dist <- 0.3
# args$outdir <- file.path(io$basedir,"/results/rna/clustering")
## END TEST ##

###################
## Load metadata ##
###################

sample_metadata <- fread(args$metadata) %>%
  .[pass_rnaQC==TRUE]
table(sample_metadata$sample)

#################
## Load Seurat ##
#################

seurat <- readRDS(io$rna.seurat)[,sample_metadata$cell]

# Add metadata to the Seurat object
foo <- sample_metadata %>% as.data.frame %>% tibble::column_to_rownames("cell") %>%
  .[colnames(seurat),]
stopifnot(colnames(seurat) == rownames(foo))
seurat@meta.data <- foo

###################
## Normalisation ##
###################

seurat <- SCTransform(seurat, 
  variable.features.n = args$nfeatures, 
  vars.to.regress = c("nFeature_RNA"),
  verbose = FALSE
)

##############################
## Dimensionality reduction ##
##############################

seurat <- RunPCA(seurat, verbose = F)
# seurat <- RunTSNE(seurat, perplexity=50)
seurat <- RunUMAP(seurat, dims=1:args$npcs, n.neighbors = args$n_neighbors, min.dist=args$min_dist, verbose=F)

################
## Clustering ##
################

seurat <- FindNeighbors(seurat, k.param = args$knn)
seurat <- FindClusters(object = seurat, resolution = args$clustering_resolution)

###################
## Visualisation ##
###################

Idents(seurat) <- sprintf("SCT_snn_res.%s",args$clustering_resolution)
# Idents(seurat) <- "eight_cell_like_jasmin"
# Idents(seurat) <- "eight_cell_like_ricard"

# DimPlot(seurat, label = FALSE, reduction = 'tsne', pt.size = 1) + NoLegend() + NoAxes()
p <- DimPlot(seurat, label = FALSE, reduction = 'umap', pt.size = 1) + NoLegend() + NoAxes()


# outfile <- sprintf("%s/umap_by_cluster_nfeatures%s_npcs%s_resolution%s_nneig%s_mindist%s.pdf",args$outdir, args$nfeatures,args$npcs,args$clustering_resolution, args$n_neighbors,args$min_dist)
outfile <- sprintf("%s/umap_by_cluster.pdf",args$outdir)
pdf(outfile, width=7, height=5)
print(p)
dev.off()


####################
## Plot ZGA genes ##
####################

genes.to.plot <- c("KHDC1L","ZSCAN4","TRIM49","TRIM49B","LEUTX","SLC34A2")
p <- FeaturePlot(seurat, features = genes.to.plot, reduction='umap') & NoAxes() & NoLegend()

# outfile <- sprintf("%s/umap_by_ZGAgenes_nfeatures%s_npcs%s_resolution%s_nneig%s_mindist%s.pdf",args$outdir, args$nfeatures,args$npcs,args$clustering_resolution, args$n_neighbors,args$min_dist)
outfile <- sprintf("%s/umap_by_ZGAgenes.pdf",args$outdir)
pdf(outfile, width=10, height=8)
print(p)
dev.off()

##########
## Save ##
##########

tmp <- seurat@meta.data[,c(sprintf("SCT_snn_res.%s",args$clustering_resolution)),drop=F] %>% as.data.table(keep.rownames = T) %>%
  setnames(c("cell","cluster")) %>%
  .[,cluster:=paste0("cluster_",cluster)]

sample_metadata.to_save <- fread(args$metadata) %>% merge(tmp,by="cell", all.x=T)

# mean(is.na(sample_metadata.to_save$cluster))
# table(sample_metadata$sample)

fwrite(sample_metadata.to_save, sprintf("%s/sample_metadata_after_clustering.txt.gz",args$outdir), sep="\t", na="NA", quote=F)
