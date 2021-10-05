# https://www.ArchRProject.com/bookdown/how-does-archr-make-pseudo-bulk-replicates.html

here::i_am("atac/archR/pseudobulk/2_archR_pseudobulk_celltypes.R")

#####################
## Define settings ##
#####################

source(here::here("settings.R"))
source(here::here("utils.R"))

# I/O
io$outdir <- file.path(io$archR.directory,"pseudobulk"); dir.create(io$outdir, showWarnings=F)

# Options
opts$matrices.to.pseudobulk <- c("PeakMatrix","GeneScoreMatrix_distal","GeneScoreMatrix_TSS")
opts$group.by <- "eight_cell_like_ricard"

########################
## Load cell metadata ##
########################

sample_metadata <- fread(io$metadata) %>%
  .[pass_atacQC==TRUE & !is.na(eight_cell_like_ricard)] %>%
  .[sample%in%opts$samples]

########################
## Load ArchR project ##
########################

source(here::here("atac/archR/load_archR_project.R"))

# Subset ArchR
ArchRProject.filt <- ArchRProject[sample_metadata$cell]
table(getCellColData(ArchRProject.filt,"Sample")[[1]])

###########################
## Update ArchR metadata ##
###########################

foo <- sample_metadata %>% 
  .[cell%in%rownames(ArchRProject.filt)] %>% setkey(cell) %>% .[rownames(ArchRProject.filt)] %>%
  as.data.frame() %>% tibble::column_to_rownames("cell")
stopifnot(all(rownames(foo) == rownames(getCellColData(ArchRProject.filt))))
ArchRProject.filt <- addCellColData(
  ArchRProject.filt,
  data = foo[[opts$group.by]], 
  name = opts$group.by,
  cells = rownames(foo),
  force = TRUE
)

###################################################
## Pseudobulk into a SummarizedExperiment object ##
###################################################

if (is.null(opts$matrices.to.pseudobulk)) {
	opts$matrices.to.pseudobulk <- getAvailableMatrices(ArchRProject.filt)
}


stopifnot(opts$matrices.to.pseudobulk%in%getAvailableMatrices(ArchRProject.filt))

se_list <- list()
for (i in opts$matrices.to.pseudobulk) {
  
  # summarise
  se_list[[i]] <- getGroupSE(ArchRProject.filt, groupBy = opts$group.by, useMatrix = i, divideN = TRUE)
  
  # save
  outfile <- sprintf("%s/pseudobulk_%s_normalised_summarized_experiment.rds",io$outdir,i)
  saveRDS(se_list[[i]], outfile)
}

##########
## Test ##
##########

# tmp <- rowMeans(assay(se_list$GeneScoreMatrix_nodistal_v2))
# names(tmp) <- rowData(se_list$GeneScoreMatrix_nodistal_v2)$name
# sort(tmp) %>% tail
