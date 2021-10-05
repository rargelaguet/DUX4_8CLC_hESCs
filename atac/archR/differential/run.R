#########
## I/O ##
#########

if (grepl("ricard",Sys.info()['nodename'])) {
  source("/Users/ricard/gastrulation_multiome_10x/settings.R")
  io$script <- "/Users/ricard/gastrulation_multiome_10x/atac/archR/differential/archr_differential_accessibility_peaks.R"
} else if(grepl("ebi",Sys.info()['nodename'])){
  source("/homes/ricard/gastrulation_multiome_10x/settings.R")
  io$script <- "/homes/ricard/gastrulation_multiome_10x/atac/archR/differential/archr_differential_accessibility_peaks.R"
  io$tmpdir <- paste0(io$basedir,"/results/atac/archR/differential/tmp")
} else {
  stop("Computer not recognised")
}
io$outdir <- paste0(io$basedir,"/results/atac/archR/differential")

#############
## Options ##
#############

# Statistical test
opts$statistical.test <- "wilcoxon"

# Define matrix
opts$matrix <- "PeakMatrix"

opts$ignore.small.celltypes <- TRUE

# Cell types
opts$celltypes <- c(
  "Epiblast",
  "Primitive_Streak",
  "Caudal_epiblast",
  "PGC",
  "Anterior_Primitive_Streak",
  "Notochord",
  "Def._endoderm",
  "Gut",
  "Nascent_mesoderm",
  "Mixed_mesoderm",
  "Intermediate_mesoderm",
  "Caudal_Mesoderm",
  "Paraxial_mesoderm",
  "Somitic_mesoderm",
  "Pharyngeal_mesoderm",
  "Cardiomyocytes",
  "Allantois",
  "ExE_mesoderm",
  "Mesenchyme",
  "Haematoendothelial_progenitors",
  "Endothelium",
  "Blood_progenitors_1",
  "Blood_progenitors_2",
  "Erythroid1",
  "Erythroid2",
  "Erythroid3",
  "NMP",
  "Rostral_neurectoderm",
  "Caudal_neurectoderm",
  "Neural_crest",
  "Forebrain_Midbrain_Hindbrain",
  "Spinal_cord",
  "Surface_ectoderm",
  "Visceral_endoderm",
  "ExE_endoderm",
  "ExE_ectoderm",
  "Parietal_endoderm"
)

########################
## Load cell metadata ##
########################

sample_metadata <- fread(io$metadata) %>%
  .[pass_atacQC==TRUE] %>%
  .[sample%in%opts$samples & celltype.predicted%in%opts$celltypes]
table(sample_metadata$celltype.predicted) %>% sort


# subset celltypes with sufficient number of cells
if (opts$ignore.small.celltypes) {
  opts$min.cells <- 50
  sample_metadata <- sample_metadata %>%
    .[,N:=.N,by=c("celltype.predicted")] %>% .[N>opts$min.cells] %>% .[,N:=NULL]
}
opts$celltypes <- unique(sample_metadata$celltype.predicted)# %>% head(n=3)

#########
## Run ##
#########

for (i in 1:length(opts$celltypes)) {
  for (j in i:length(opts$celltypes)) {
    if (i!=j) {
      groupA <- opts$celltypes[[i]]
      groupB <- opts$celltypes[[j]]
      # print(sprintf("i=%s (%s), j=%s (%s)",i,groupA,j,groupB))
      
      outfile <- sprintf("%s/%s_%s_vs_%s.txt.gz", io$outdir,opts$matrix,groupA,groupB)
      if (!file.exists(outfile)) {
        # Define LSF command
        if (grepl("ricard",Sys.info()['nodename'])) {
          lsf <- ""
        } else if (grepl("ebi",Sys.info()['nodename'])) {
          lsf <- sprintf("bsub -M 7000 -n 1 -o %s/%s_%s_vs_%s.txt", io$tmpdir,opts$matrix,groupA,groupB)
        }
        cmd <- sprintf("%s Rscript %s --groupA %s --groupB %s --matrix %s --test %s --outfile %s", lsf, io$script, groupA, groupB, opts$matrix, opts$statistical.test, outfile)
        
        # Run
        print(cmd)
        system(cmd)
      }
    }
  }
}

