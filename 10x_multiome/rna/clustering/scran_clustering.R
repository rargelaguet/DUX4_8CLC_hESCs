suppressPackageStartupMessages(library(SingleCellExperiment))
suppressPackageStartupMessages(library(scater))
suppressPackageStartupMessages(library(scran))
suppressPackageStartupMessages(library(argparse))
suppressPackageStartupMessages(library(igraph))
suppressPackageStartupMessages(library(edgeR))

######################
## Define arguments ##
######################

p <- ArgumentParser(description='')
p$add_argument('--sce',         type="character",                            help='SingleCellExperiment file')
p$add_argument('--metadata',    type="character",                            help='metadata file')
p$add_argument('--samples',     type="character",                nargs='+',  help='Query batch(es)')
p$add_argument('--features',    type="integer",    default=1000,             help='Number of cores')
p$add_argument('--npcs',        type="integer",    default=30,               help='Number of cores')
p$add_argument('--n_neighbors', type="integer",    default=30,               help='(UMAP) Number of neighbours')
p$add_argument('--knn',         type="integer",    default=30,               help='Number of k nearest neighbours for the graph')
p$add_argument('--nclusters',    type="integer",    default=10,               help='Number of clusters')
p$add_argument('--min_dist',    type="double",     default=0.3,              help='(UMAP) Minimum distance')
p$add_argument('--outdir',      type="character",                                  help='Output directory')

args <- p$parse_args(commandArgs(TRUE))

#####################
## Define settings ##
#####################

source(here::here("settings.R"))
source(here::here("utils.R"))

## START TEST ##
args$samples <- opts$samples
args$sce <- io$rna.sce
args$metadata <- file.path(io$basedir,"results/rna/qc/sample_metadata_after_qc.txt.gz")
args$features <- 2000
args$npcs <- 25
args$nclusters <- 10
args$colour_by <- c("sample")
args$outdir <- file.path(io$basedir,"/results/rna/clustering")
## END TEST ##

##########################
## Load sample metadata ##
##########################

sample_metadata <- fread(args$metadata) %>%
  .[pass_rnaQC==TRUE & sample%in%args$samples]
table(sample_metadata$sample)

stopifnot(args$colour_by %in% colnames(sample_metadata))

###############
## Load data ##
###############

# Load RNA expression data as SingleCellExperiment object
sce <- load_SingleCellExperiment(args$sce, cells=sample_metadata$cell, normalise = TRUE)

# Add cell metadata to colData
colData(sce) <- sample_metadata %>% tibble::column_to_rownames("cell") %>% DataFrame

#######################
## Feature selection ##
#######################

decomp <- modelGeneVar(sce)
decomp <- decomp[decomp$mean > 0.01,]

hvgs <- decomp[order(decomp$FDR),] %>% head(n=args$features) %>% rownames
# hvgs <- rownames(decomp)[decomp$p.value <= 0.05]

# Subset SingleCellExperiment
sce_filt <- sce[hvgs,]
dim(sce_filt)

##############################
## Dimensionality reduction ##
##############################

# PCA
sce_filt <- runPCA(sce_filt, ncomponents = args$npcs, ntop=args$features)

# UMAP
set.seed(42)
sce_filt <- runUMAP(sce_filt, dimred="PCA", n_neighbors = args$n_neighbors, min_dist = args$min_dist)
reducedDim(sce, "UMAP") <- reducedDim(sce_filt, "UMAP")

##########
## Plot ##
##########

# to.plot <- reducedDim(sce_filt,"UMAP") %>% as.data.table %>% 
#   .[,cell:=colnames(sce_filt)] %>%
#   merge(sample_metadata, by="cell")
# 
# p <- ggplot(to.plot, aes_string(x="V1", y="V2", fill="sample")) +
#   geom_point(size=1.5, shape=21, stroke=0.05) +
#   # scale_fill_manual(values=opts$celltype.colors) +
#   theme_classic() +
#   theme(
#     axis.title = element_blank(),
#     axis.text = element_blank(),
#     axis.ticks = element_blank(),
#     legend.position="none",
#     legend.title=element_blank()
#   )
# 
# pdf(sprintf("%s/umap_by_sample.pdf",args$outdir), width=7, height=5)
# print(p)
# dev.off()


################
## Clustering ##
################

# Build KNN graph
graph <- buildSNNGraph(sce_filt, k = args$knn, use.dimred = "PCA", type = "rank")

# Detect communities in the graph
walktrap <- cluster_walktrap(graph)
# louvain <- cluster_louvain(graph)

# Clustering
cluster <- cut_at(walktrap, no = args$nclusters)
sce_filt$cluster <- factor(cluster)

# Plot and colour by cluster
p <- plotUMAP(sce_filt, colour_by="cluster")
outfile <- sprintf("%s/%s_umap_%d_%d_%d_%s_%s.pdf",args$outdir, paste(args$samples,collapse="-"), args$features,args$npcs,args$n_neighbors,args$min_dist,i)
pdf(outfile, width=7, height=5)
print(p)
dev.off()


##################
## Find markers ##
##################

markers <- findMarkers(sce_filt, 
  groups = sce_filt$cluster, 
  test.type = "wilcox", 
  direction = "up", 
  pval.type = c("some"),
  lfc = 1
)

# create data.frame with marker genes
dt <- names(markers) %>% 
  map(function(i) { markers[[i]] %>% 
        as.data.table(keep.rownames = T) %>% 
        setnames("rn","gene") %>% 
        .[,c("gene","FDR")] %>% 
        .[FDR<1e-2] %>%
        .[,FDR:=round(FDR,digits=5)] %>%
        .[,cluster:=i]
    }) %>% rbindlist
table(dt$cluster)


# save
fwrite(dt, sprintf("%s/marker_genes.txt.gz",args$outdir), sep="\t")

# Plot markers
opts$ngenes.to.plot <- 10

for (i in unique(sce_filt$cluster)) {
  markers.to.plot <- markers[[i]] %>% head(n=opts$ngenes.to.plot) %>% rownames
  for (j in markers.to.plot) {
    p <- plotUMAP(sce, colour_by=j)
    pdf(sprintf("%s/umap_cluster_%s_expr_%s.pdf",args$outdir,i,j), width=7, height=5)
    print(p)
    dev.off()
  }
}

markers <- findMarkers(
  sce_filt, 
  groups = sce_filt$cluster, 
  test.type = "wilcox", 
  lfc = 1
)

#############################
## Differential expression ##
#############################

# opts$groupA <- c("1")
# opts$groupB <- c("5")
# 
# # Convert SCE to DGEList
# sce_edger <- convertTo(sce_filt[,sce_filt$cluster%in%c(opts$groupA,opts$groupB)], type="edgeR")
# dim(sce_edger)
# 
# # Define design matrix (with intercept)
# design <- model.matrix(~droplevels(sce_edger$samples$cluster))
# 
# # Estimate dispersions
# sce_edger  <- estimateDisp(sce_edger,design)
# 
# # Fit GLM
# fit <- glmQLFit(sce_edger,design)
# 
# # Likelihood ratio test
# lrt <- glmQLFTest(fit)
# 
# # Construct output data.frame
# out <- topTags(lrt, n=nrow(lrt))$table %>% as.data.table(keep.rownames=T) %>%
#   setnames(c("gene","logFC","logCPM","LR","p.value","padj_fdr")) %>%
#   .[,c("logCPM","LR"):=NULL] %>%
#   .[, sig := (padj_fdr<=0.01 & abs(logFC)>=1)] %>%
#   setorder(-sig, padj_fdr, na.last=T)
# 
# # Plot differentially expressed genes
# opts$ngenes.to.plot <- 10
# 
# markers.to.plot <- out[sig==T & logFC>1,gene] %>% head(opts$ngenes.to.plot)
# for (i in markers.to.plot) {
#   p <- plotUMAP(sce, colour_by=i)
#   pdf(sprintf("%s/differential/umap_up_%s.pdf",args$outdir,i), width=7, height=5)
#   print(p)
#   dev.off()
# }
# 
# markers.to.plot <- out[sig==T & logFC<(-1),gene] %>% head(opts$ngenes.to.plot)
# for (i in markers.to.plot) {
#   p <- plotUMAP(sce, colour_by=i)
#   pdf(sprintf("%s/differential/umap_down_%s.pdf",args$outdir,i), width=7, height=5)
#   print(p)
#   dev.off()
# }

##########
## TEST ##
##########

# genes.to.plot <- rownames(sce_filt)[grep("Hb",rownames(sce_filt))]
# # genes.to.plot <- c("","","","")
# genes.to.plot <- genes.to.plot[genes.to.plot%in%rownames(sce)]
# for (i in genes.to.plot) {
#   p <- plotUMAP(sce, colour_by=i) +
#     labs(color=i) +
#     scale_color_gradient(low = "gray80", high = "red")
#   print(p)
# }
