# https://www.ArchRProject.com/bookdown/how-does-archr-make-pseudo-bulk-replicates.html

here::i_am("atac/archR/pseudobulk/2_archR_pseudobulk_celltypes.R")

#####################
## Define settings ##
#####################

source(here::here("settings.R"))
source(here::here("utils.R"))

# I/O
io$outdir <- file.path(io$archR.directory,"pseudobulk")

# Options
opts$matrices.to.pseudobulk <- "DeviationMatrix_Motif_cisbp_lenient" # c("PeakMatrix","GeneScoreMatrix","GeneScoreMatrix_nodistal")
opts$group.by <- "celltype.predicted"

########################
## Load cell metadata ##
########################

sample_metadata <- fread(io$metadata) %>%
  .[pass_atacQC==TRUE & doublet_call==FALSE & !is.na(celltype.predicted)] %>%
  .[sample%in%opts$samples]

########################
## Load ArchR project ##
########################

source(here::here("atac/archR/load_archR_project.R"))

# Subset ArchR
ArchRProject.filt <- ArchRProject[sample_metadata$cell]
table(getCellColData(ArchRProject.filt,"Sample")[[1]])

# getAvailableMatrices(ArchRProject)

###################################################
## Pseudobulk into a SummarizedExperiment object ##
###################################################

if (is.null(opts$matrices.to.pseudobulk)) {
	opts$matrices.to.pseudobulk <- getAvailableMatrices(ArchRProject.filt)
}

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
