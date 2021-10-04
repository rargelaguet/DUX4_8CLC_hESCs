suppressPackageStartupMessages(library(Seurat))
suppressPackageStartupMessages(library(argparse))

######################
## Define arguments ##
######################

p <- ArgumentParser(description='')
p$add_argument('--inputdir',        type="character",                    help='Input directory')
p$add_argument('--outputdir',       type="character",                    help='Output directory')
p$add_argument('--use_soupX',       action="store_true",                 help='use SoupX-corrected matrix?')
p$add_argument('--test',            action="store_true",                 help='Testing mode')
p$add_argument('--samples',         type="character",       nargs="+",   help='Samples')
args <- p$parse_args(commandArgs(TRUE))

#####################
## Define settings ##
#####################

if (grepl("ricard",Sys.info()['nodename'])) {
  source("/Users/ricard/gastrulation_multiome_10x/settings.R")
} else if (grepl("ebi",Sys.info()['nodename'])) {
  source("/homes/ricard/gastrulation_multiome_10x/settings.R")
}

# These settings are important for consistency with ArchR, which provides little flexibility to edit cell names
opts$trim.barcode <- FALSE
opts$sample_cell_separator <- "#"

## START TEST ##
args <- list()
args$inputdir <- paste0(io$basedir,"/original")
args$outputdir <- paste0(io$basedir,"/processed/rna_soupX")
args$samples <- opts$samples
args$use_soupX <- TRUE
args$test <- FALSE
## END TEST ##


##############################
## Load and merge data sets ##
##############################

stopifnot(args$samples%in%opts$samples)
if (args$test) args$samples <- head(args$samples,n=2)

mtx <- list()
cell.info <- list()
gene.info <- list()

for (i in args$samples) {
  print(i)
    
  # Load gene metadata
  gene.loc <- sprintf("%s/%s/filtered_feature_bc_matrix/features.tsv.gz",args$inputdir,i)
  gene.info[[i]] <- fread(gene.loc, header=F, select=c(1,2,3)) %>%
    setnames(c("ens_id","symbol","modality")) %>%
    .[,idx:=1:.N] %>%
    .[modality=="Gene Expression"]
  dim(gene.info[[i]])
  
  # Load cell metadata
  barcode.loc <- sprintf("%s/%s/filtered_feature_bc_matrix/barcodes.tsv.gz",args$inputdir,i)
  cell.info[[i]] <- fread(barcode.loc, header=F) %>%
    setnames("barcode") %>%
    .[,barcode:=ifelse(rep(opts$trim.barcode,.N),gsub("-1","",barcode),barcode)] %>%
    .[,c("sample","cell"):=list(i,sprintf("%s%s%s",i,opts$sample_cell_separator,barcode))]
  dim(cell.info[[i]])
  
  # Load matrix  
  if (args$use_soupX) {
    matrix.loc <- sprintf("%s/%s/filtered_feature_bc_matrix/soupX/soupX_adjusted_matrix.mtx.gz",args$inputdir,i)
  } else {
    matrix.loc <- sprintf("%s/%s/filtered_feature_bc_matrix/matrix.mtx.gz",args$inputdir,i)
  }
  mtx[[i]] <- Matrix::readMM(matrix.loc)[gene.info[[i]]$idx,]
  stopifnot(nrow(cell.info[[i]])==colnames(mtx[[i]]))
  rownames(mtx[[i]]) <- gene.info[[i]]$symbol
  colnames(mtx[[i]]) <- cell.info[[i]]$cell
}

# Sanity checks
# stopifnot(length(unique(lapply(gene.info,nrow)))==1)
# stopifnot(length(unique(lapply(mtx,nrow)))==1)
# stopifnot(length(unique(lapply(mtx,rownames)))==1)

#################
## Concatenate ##
#################

# Extract unique gene metadata
gene.info <- gene.info[[1]]

# Concatenate cell metadata
cell.info <- rbindlist(cell.info)
rownames(cell.info) <- cell.info$cell

# Concatenate matrices
mtx <- do.call("cbind",mtx)
colnames(mtx) <- cell.info$cell

##################
## Filter genes ##
##################

# Keep protein-coding genes
# if (!is.null(opts$subset.proteincoding)){
#     genes <- fread(opts$subset.proteincoding)[,ens_id]
#     genes <- genes[genes %in% mouse.genes]
#     mouse.genes <- mouse.genes[mouse.genes %in% genes]
#     mtx <- mtx[mouse.genes,]
# }

# Remove duplicated genes
gene.info <- gene.info[!duplicated(gene.info$symbol),]
mtx <- mtx[gene.info$symbol,]

# Sanity checks
stopifnot(sum(duplicated(rownames(mtx)))==0)
stopifnot(sum(duplicated(colnames(mtx)))==0)
stopifnot(sum(duplicated(cell.info$cell))==0)
stopifnot(all(colnames(mtx) == cell.info$cell))
stopifnot(all(rownames(mtx) == gene.info$symbol))

##########################
## Create Seurat object ##
##########################

seurat <- CreateSeuratObject(mtx, meta.data = cell.info)

# Add mitochondrial percenatge
seurat[["mitochondrial_percent_RNA"]] <- PercentageFeatureSet(seurat, pattern = "^mt-")

# Add ribosomal RNA content
ribo.genes <- c(grep(pattern = "^Rpl", x = rownames(seurat), value = TRUE), grep(pattern = "^Rps", x = rownames(seurat), value = TRUE))
seurat[["ribosomal_percent_RNA"]] <- PercentageFeatureSet(seurat, features = ribo.genes)

#####################
## Create metadata ##
#####################

metadata <- seurat@meta.data %>% as.data.table %>% .[,orig.ident:=NULL] %>%
  .[,c("cell","barcode","sample","nFeature_RNA","nCount_RNA","mitochondrial_percent_RNA","ribosomal_percent_RNA")]

metadata[,stage:=substr(sample,1,4)]

##########
## Save ##
##########

saveRDS(seurat, paste0(args$outputdir,"/seurat.rds"))
fwrite(cell.info, paste0(args$outputdir,"/cell_info.txt.gz"), quote=F, na="NA", sep="\t")
fwrite(gene.info, paste0(args$outputdir,"/gene_info.txt.gz"), quote=F, na="NA", sep="\t")
fwrite(metadata, paste0(args$outputdir,"/metadata.txt.gz"), quote=F, na="NA", sep="\t")

