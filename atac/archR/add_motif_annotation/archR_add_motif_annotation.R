#####################
## Define settings ##
#####################

source(here::here("settings.R"))

########################
## Load ArchR Project ##
########################

source(here::here("atac/archR/load_archR_project.R"))

#####################
## Define settings ##
#####################

# I/O
# io$metadata <- paste0(io$basedir,"/processed/atac/archR/sample_metadata_after_archR.txt.gz")
# io$outdir <- paste0(io$basedir,"/results/atac/archR/motif_annotations")

##########################
## Add motif annotation ##
##########################

# cisbp (stringent threshold)
ArchRProject <- addMotifAnnotations(
  ArchRProject,
  motifSet = "cisbp",
  cutOff = 5e-05,
  name = "Motif_cisbp",
  cutOff = 5e-05,
  width = 7,
  force = TRUE
)

# cisbp (lenient threshold)
ArchRProject <- addMotifAnnotations(
  ArchRProject,
  motifSet = "cisbp",
  name = "Motif_cisbp_lenient",
  cutOff = 1e-04,
  width = 7,
  force = TRUE
)

# homer
# ArchRProject <- addMotifAnnotations(
#   ArchRProject,
#   motifSet = "homer",
#   cutOff = opts$motif.pvalue.cutoff,
#   name = "Motif_homer",
#   force = TRUE
# )


# JASPAR2018 human
# ArchRProject <- addMotifAnnotations(
#   ArchRProject, 
#   motifSet = "JASPAR2018",      
#   collection = "CORE",  
#   species = "Homo sapiens",
#   cutOff = opts$motif.pvalue.cutoff,   
#   name = "Motif_JASPAR2018_human",
#   force = TRUE
# )

# JASPAR2018 mouse
# ArchRProject <- addMotifAnnotations(
#   ArchRProject, 
#   motifSet = "JASPAR2018",      
#   collection = "CORE",  
#   species = "Mus musculus",
#   cutOff = opts$motif.pvalue.cutoff,   
#   name = "Motif_JASPAR2018_human",
#   force = TRUE
# )

# JASPAR2020 human
ArchRProject <- addMotifAnnotations(
  ArchRProject, 
  motifSet = "JASPAR2020",      
  collection = "CORE",  
  species = "Homo sapiens",
  cutOff = 5e-05,   
  name = "Motif_JASPAR2020_human",
  force = TRUE
)

# JASPAR2020 mouse
# ArchRProject <- addMotifAnnotations(
#   ArchRProject, 
#   motifSet = "JASPAR2020",      
#   collection = "CORE",  
#   species = "Mus musculus",
#   cutOff = opts$motif.pvalue.cutoff,   
#   name = "Motif_JASPAR2020_mouse",
#   force = TRUE
# )

################################
## Save peakAnnotation object ##
################################

saveRDS(ArchRProject@peakAnnotation, sprintf("%s/Annotations/peakAnnotation.rds",io$archR.directory))

#############
## Explore ##
#############

# Motif_cisbp <- ArchRProject@peakAnnotation[["Motif_cisbp"]][["motifSummary"]]
# Motif_JASPAR2020_human <- ArchRProject@peakAnnotation[["Motif_JASPAR2020_human"]][["motifSummary"]]
# Motif_JASPAR2020_mouse <- ArchRProject@peakAnnotation[["Motif_JASPAR2020_mouse"]][["motifSummary"]]
# Motif_homer <- ArchRProject@peakAnnotation[["Motif_homer"]][["motifSummary"]]
