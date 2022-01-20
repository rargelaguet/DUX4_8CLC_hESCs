#########
## I/O ##
#########

setwd("/bi/group/reik/ricard/scripts/gastrulation_multiome_10x")

source(here::here("settings.R"))

# io$Rscript <- "/Library/Frameworks/R.framework/Versions/Current/Resources/bin/Rscript"
io$Rscript <- "Rscript"
io$script <- here::here("rna/differential/differential.R")
io$tmpdir <- file.path(io$basedir,"results/rna/differential/tmp"); dir.create(io$tmpdir, showWarnings=F)
io$outdir <- file.path(io$basedir,"results/rna/differential/test")

##########################
## Load sample metadata ##
##########################

# io$metadata <- paste0(io$basedir,"/results/rna/celltype_denoising/sample_metadata_after_celltype_denoising.txt.gz")
sample_metadata <- fread(io$metadata) %>% 
  .[pass_rnaQC==TRUE & doublet_call==FALSE & !is.na(celltype.mapped)]

#############
## Options ##
#############

# Define samples
# opts$samples <- c(
#   "E7.5_rep1",
#   "E7.5_rep2",
#   "E8.5_rep1",
#   "E8.5_rep2"
# )

# Define stages
opts$stages <- c(
  "E7.5",
  "E8.0",
  "E8.5"
)

# Testing mode
opts$test_mode <- FALSE

# Define cell types
opts$groups <- c("Epiblast","Primitive_Streak")
# opts$groups <- names(which(table(sample_metadata[stage%in%opts$stages,celltype.mapped])>=25))

###################################
## Run all pair-wise comparisons ##
###################################

for (i in 1:length(opts$groups)) {
  groupA <- opts$groups[[i]]
  for (j in i:length(opts$groups)) {
    if (i!=j) {
      groupB <- opts$groups[[j]]
      outfile <- sprintf("%s/%s_vs_%s.txt.gz", io$outdir,groupA,groupB)
      
      if (!file.exists(outfile)) {
        
        # Define LSF command
        # if (grepl("ricard",Sys.info()['nodename'])) {
        #   lsf <- ""
        # } else if (grepl("ebi",Sys.info()['nodename'])) {
        #   lsf <- sprintf("bsub -M 9000 -n 1 -o %s/%s_vs_%s.txt", io$tmpdir,groupA,groupB)
        # }  else if (grepl("babraham",Sys.info()['nodename'])) {
        #   lsf <- sprintf("ssub -m 5G -c 1 -o %s/%s_vs_%s.txt", io$tmpdir,groupA,groupB)
        # }
        lsf <- ""
        cmd <- sprintf("%s %s %s --stages %s --groupA %s --groupB %s --outfile %s", lsf, io$Rscript, io$script, paste(opts$stages, collapse=" "), groupA, groupB, outfile)
        if (isTRUE(opts$test_mode)) cmd <- paste0(cmd, " --test_mode")
        
        # Run
        print(cmd)
        system(cmd)
      }
    }
  }
}


##############################
## Run selected comparisons ##
##############################

# opts$comparisons <- list(
#   c("groupA"="Mixed_mesoderm",          "groupB"="ExE_ectoderm")
#   # c("groupA"="Blood_progenitors_2", "groupB"="Caudal_Mesoderm"),
#   # c("groupA"="Allantois",           "groupB"="Haematoendothelial_progenitors"),
#   # c("groupA"="Blood_progenitors_2", "groupB"="NMP"),
#   # c("groupA"="Erythroid1",          "groupB"="Anterior_Primitive_Streak"),
#   # c("groupA"="Erythroid1",          "groupB"="Haematoendothelial_progenitors")
# )
# 
# for (comparison in opts$comparisons) {
#   groupA <- comparison[["groupA"]]; groupB <- comparison[["groupB"]]
#   for (test in opts$statistical.test) {
#     outfile <- sprintf("%s/%s_vs_%s.txt.gz", io$outdir,groupA,groupB)
#     
#     # Define LSF command
#     if (grepl("ricard",Sys.info()['nodename'])) {
#       lsf <- ""
#     } else if (grepl("ebi",Sys.info()['nodename'])) {
#       lsf <- sprintf("bsub -M 15000 -n 1 -q research-rh74 -o %s/%s_vs_%s.txt", io$tmpdir,groupA,groupB)
#     }
#     cmd <- sprintf("%s Rscript %s --stages %s --groupA %s --groupB %s --test %s --outfile %s", lsf, io$script, paste(opts$stages, collapse=" "), groupA, groupB, test, outfile)
#     if (isTRUE(opts$test_mode)) cmd <- paste0(cmd, " --test_mode")
#     
#     # Run
#     print(cmd)
#     system(cmd)
#   }
# }
# 
