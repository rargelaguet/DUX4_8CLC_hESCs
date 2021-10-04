suppressPackageStartupMessages(library(Seurat))
suppressPackageStartupMessages(library(argparse))

######################
## Define arguments ##
######################

p <- ArgumentParser(description='')
p$add_argument('--inputdir',        type="character",                    help='Input directory')
p$add_argument('--outputdir',       type="character",                    help='Output directory')
p$add_argument('--test',            action="store_true",                 help='Testing mode')
p$add_argument('--samples',         type="character",       nargs="+",   help='Samples')
args <- p$parse_args(commandArgs(TRUE))

#####################
## Define settings ##
#####################

here::i_am("rna/processing/1_create_seurat_rna.R")
source(here::here("settings.R"))

# These settings are important for consistency with ArchR, which provides little flexibility to edit cell names
opts$trim.barcode <- FALSE
opts$sample_cell_separator <- "#"

## START TEST ##
args <- list()
args$inputdir <- file.path(io$basedir,"original/rna")
args$outputdir <- file.path(io$basedir,"processed/rna")
args$samples <- opts$samples
args$test <- FALSE
## END TEST ##


##############################
## Load and merge data sets ##
##############################

stopifnot(args$samples%in%opts$samples)
if (args$test) args$samples <- head(args$samples,n=2)

count_mtx <- list()
cell.info <- list()
gene.info <- list()

i <- args$samples[1]
for (i in args$samples) {
  print(i)
    
  count_mtx[[i]] <- Read10X(sprintf("%s/%s",args$inputdir,i))[["Gene Expression"]]

  # Load gene metadata
  # gene.loc <- sprintf("%s/%s/genes.tsv.gz",args$inputdir,i)
  # gene.info[[i]] <- fread(gene.loc, header=F) %>%
  #   setnames(c("ens_id","symbol")) %>%
  #   .[,idx:=1:.N]
  # dim(gene.info[[i]])
  
  # Load cell metadata
  # barcode.loc <- sprintf("%s/%s/barcodes.tsv.gz",args$inputdir,i)
  # cell.info[[i]] <- fread(barcode.loc, header=F) %>%
  #   setnames("barcode") %>%
  #   .[,barcode:=ifelse(rep(opts$trim.barcode,.N),gsub("-1","",barcode),barcode)] %>%
  #   .[,c("sample","cell"):=list(i,sprintf("%s%s%s",i,opts$sample_cell_separator,barcode))]
  # dim(cell.info[[i]])
  
  # Load matrix  
  # matrix.loc <- sprintf("%s/%s/matrix.mtx.gz",args$inputdir,i)
  # count_mtx[[i]] <- Matrix::readMM(matrix.loc)[gene.info[[i]]$idx,]
  # stopifnot(nrow(cell.info[[i]])==ncol(count_mtx[[i]]))
  # rownames(count_mtx[[i]]) <- gene.info[[i]]$symbol
  # colnames(count_mtx[[i]]) <- cell.info[[i]]$cell
  
  # Basic filtering
  # count_mtx[[i]] <- count_mtx[[i]][,colSums(count_mtx[[i]])>=500]
}

print(lapply(count_mtx,dim))

#######################
## Keep common genes ##
#######################

genes <- Reduce("intersect",lapply(count_mtx,rownames))
for (i in 1:length(count_mtx)) {
  count_mtx[[i]] <- count_mtx[[i]][genes,]
}

stopifnot(length(unique(lapply(count_mtx,nrow)))==1)
stopifnot(length(unique(lapply(count_mtx,rownames)))==1)

#################
## Concatenate ##
#################

# Extract unique gene metadata
# gene.info <- gene.info[[1]]

# Concatenate cell metadata
cell.info <- rbindlist(cell.info)
rownames(cell.info) <- cell.info$cell

# Concatenate matrices
count_mtx <- do.call("cbind",count_mtx)
# colnames(count_mtx) <- cell.info$cell

##################
## Filter genes ##
##################

# Keep protein-coding genes
# if (!is.null(opts$subset.proteincoding)){
#     genes <- fread(opts$subset.proteincoding)[,ens_id]
#     genes <- genes[genes %in% mouse.genes]
#     mouse.genes <- mouse.genes[mouse.genes %in% genes]
#     count_mtx <- count_mtx[mouse.genes,]
# }

# Remove duplicated genes
count_mtx <- count_mtx[!duplicated(rownames(count_mtx)),]

# Sanity checks
stopifnot(sum(duplicated(rownames(count_mtx)))==0)
stopifnot(sum(duplicated(colnames(count_mtx)))==0)

##########################
## Create Seurat object ##
##########################

cell.info.to.seurat <- cell.info[cell%in%colnames(count_mtx)] %>% setkey(cell) %>% .[colnames(count_mtx)] %>% as.data.frame
rownames(cell.info.to.seurat) <- cell.info.to.seurat$cell
stopifnot(rownames(cell.info.to.seurat)==colnames(count_mtx))
stopifnot(sum(is.na(rownames(cell.info.to.seurat$cell)))==0)

seurat <- CreateSeuratObject(count_mtx, meta.data = cell.info.to.seurat)

head(seurat@meta.data)

# Add mitochondrial percenatge
seurat[["mitochondrial_percent_RNA"]] <- PercentageFeatureSet(seurat, pattern = "^MT-")

# Add ribosomal RNA content
ribo.genes <- grep(pattern = "^RP[L|S]", x = rownames(seurat), value = TRUE)
seurat[["ribosomal_percent_RNA"]] <- PercentageFeatureSet(seurat, features = ribo.genes)

#####################
## Create metadata ##
#####################

metadata <- seurat@meta.data %>% as.data.table %>% .[,orig.ident:=NULL] %>%
  .[,c("cell","barcode","sample","nFeature_RNA","nCount_RNA","mitochondrial_percent_RNA","ribosomal_percent_RNA")]

# metadata[,stage:=substr(sample,1,4)]
# metadata[,stage:="E8.5"]

##########
## Save ##
##########

fwrite(metadata, paste0(args$outputdir,"/metadata.txt.gz"), quote=F, na="NA", sep="\t")
# fwrite(cell.info, paste0(args$outputdir,"/cell_info.txt.gz"), quote=F, na="NA", sep="\t")
# fwrite(gene.info, paste0(args$outputdir,"/gene_info.txt.gz"), quote=F, na="NA", sep="\t")
saveRDS(seurat, paste0(args$outputdir,"/seurat.rds"))

