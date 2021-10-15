
################
## Define I/O ##
################

if (grepl("ricard",Sys.info()['nodename'])) {
  source("/Users/ricard/gastrulation_multiome_10x/settings.R")
  io$script <- "/Users/ricard/gastrulation_multiome_10x/atac/archR/dimensionality_reduction/automated/archR_dimensionality_reduction.R"
} else {
  source("/homes/ricard/gastrulation_multiome_10x/settings.R")
  io$script <- "/homes/ricard/gastrulation_multiome_10x/atac/archR/dimensionality_reduction/automated/archR_dimensionality_reduction.R"
}
io$outdir <- paste0(io$basedir,"/results/atac/archR/dimensionality_reduction/PeakMatrix"); dir.create(io$outdir)
io$tmpdir <- paste0(io$outdir,"/tmp"); dir.create(io$tmpdir)

####################
## Define options ##
####################

# Samples
opts$samples <- c(
  "E7.5_rep1",
  "E7.5_rep2",
  "E8.0_rep1",
  "E8.0_rep2",
  "E8.5_rep1",
  "E8.5_rep2"
)

# input Matrix
# opts$matrix <- "PeakMatrix"
opts$matrix <- "PeakMatrix"

# Number of highly variable peaks
# opts$nfeatures <- c(1000,2000,3000)
opts$nfeatures <- c(25000,50000)

# Number of LSI dimensions
# opts$ndims <- c(25,50)
opts$ndims <- c(50)

# Variable to do MNN batch correction on
opts$batch.variable <- "stage"

# UMAP hyperparameters
# opts$n_neighbors <- c(20,30,40)
opts$n_neighbors <- c(30,45)

# opts$min_dist <- c(0.20,0.30,0.40)
opts$min_dist <- c(0.30,0.45)


##############################
## Run one sample at a time ##
##############################

# # LSF params
# opts$memory <- 7000

# # opts$colour_by <- c("celltype.mapped","log_nFrags_atac","doublet_call")
# opts$colour_by <- c("celltype.predicted","log_nFrags_atac")

# for (i in opts$samples) {
#   outdir <- sprintf("%s/%s",io$outdir,i); dir.create(outdir, showWarnings = F)
#   for (j in opts$nfeatures) {
#     for (k in opts$ndims) {
      
#       # Define LSF command
#       if (grepl("ricard",Sys.info()['nodename'])) {
#         lsf <- ""
#       } else if (grepl("ebi",Sys.info()['nodename'])) {
#         lsf <- sprintf("bsub -M %s -n 1 -o %s/%s_%d_%d.txt", opts$memory, io$tmpdir,i,j,k)
#       }
#       cmd <- sprintf("%s Rscript %s --samples %s --nfeatures %d --ndims %d --matrix %s --n_neighbors %s --min_dist %s --colour_by %s --outdir %s",
#                      lsf, io$script, i, j, k, opts$matrix, paste(opts$n_neighbors,collapse=" "), paste(opts$min_dist,collapse=" "), paste(opts$colour_by,collapse=" "), outdir)
      
#       # Run
#       print(cmd)
#       system(cmd)
#     }
#   }
# }

###########################################################
## Run one stage at a time (no batch correction applied) ##
###########################################################

# # LSF params
# opts$memory <- 10000

# # opts$colour_by <- c("celltype.mapped","sample","log_nFrags_atac","doublet_call")
# opts$colour_by <- c("celltype.predicted","sample")

# for (i in opts$stages) {
#   samples <- opts$samples[grep(i,opts$samples)]
#   outdir <- sprintf("%s/%s",io$outdir,i); dir.create(outdir, showWarnings = F)
#   for (j in opts$nfeatures) {
#     for (k in opts$ndims) {

#       # Define LSF command
#       if (grepl("ricard",Sys.info()['nodename'])) {
#         lsf <- ""
#       } else if (grepl("ebi",Sys.info()['nodename'])) {
#         lsf <- sprintf("bsub -M %s -n 1 -o %s/%s_%d_%d.txt", opts$memory, io$tmpdir,i,j,k)
#       }
#       cmd <- sprintf("%s Rscript %s --samples %s --nfeatures %d --ndims %d --matrix %s --n_neighbors %s --min_dist %s --colour_by %s --outdir %s",
#                      lsf, io$script, paste(samples,collapse=" "), j, k, opts$matrix, paste(opts$n_neighbors,collapse=" "), paste(opts$min_dist,collapse=" "), paste(opts$colour_by,collapse=" "), outdir)

#       # Run
#       print(cmd)
#       system(cmd)
#     }
#   }
# }

##################################################
## Run all samples at once, no batch correction ##
##################################################

# # LSF params
# opts$memory <- 15000

# outdir <- sprintf("%s/all_cells",io$outdir); dir.create(outdir, showWarnings = F)
# # opts$colour_by <- c("celltype.mapped","sample","stage","log_nFrags_atac","doublet_call")
# opts$colour_by <- c("celltype.predicted","stage")

# for (j in opts$nfeatures) {
#   for (k in opts$ndims) {
    
#     # Define LSF command
#     if (grepl("ricard",Sys.info()['nodename'])) {
#       lsf <- ""
#     } else if (grepl("ebi",Sys.info()['nodename'])) {
#       lsf <- sprintf("bsub -M %s -n 1 -o %s/allsamples_%d_%d.txt", opts$memory, io$tmpdir,j,k)
#     }
#     cmd <- sprintf("%s Rscript %s --samples %s --nfeatures %d --ndims %d --matrix %s --n_neighbors %s --min_dist %s --colour_by %s --outdir %s",
#                    lsf, io$script, paste(opts$samples,collapse=" "), j, k, opts$matrix, paste(opts$n_neighbors,collapse=" "), paste(opts$min_dist,collapse=" "), paste(opts$colour_by,collapse=" "), outdir)
    
#     # Run
#     print(cmd)
#     system(cmd)
#   }
# }

######################################################################
## Run all samples at once, no batch correction, removing ExE cells ##
######################################################################

# LSF params
opts$memory <- 15000

outdir <- sprintf("%s/all_cells_remove_ExE_celltypes",io$outdir); dir.create(outdir, showWarnings = F)
# opts$colour_by <- c("celltype.mapped","sample","stage","log_nFrags_atac","doublet_call")
opts$colour_by <- c("celltype.predicted","stage")

for (j in opts$nfeatures) {
  for (k in opts$ndims) {
    
    # Define LSF command
    if (grepl("ricard",Sys.info()['nodename'])) {
      lsf <- ""
    } else if (grepl("ebi",Sys.info()['nodename'])) {
      lsf <- sprintf("bsub -M %s -n 1 -o %s/allsamples_removeExE_%d_%d.txt", opts$memory, io$tmpdir,j,k)
    }
    cmd <- sprintf("%s Rscript %s --samples %s --nfeatures %d --ndims %d --matrix %s --n_neighbors %s --min_dist %s --colour_by %s --remove_ExE_celltypes --outdir %s",
                   lsf, io$script, paste(opts$samples,collapse=" "), j, k, opts$matrix, paste(opts$n_neighbors,collapse=" "), paste(opts$min_dist,collapse=" "), paste(opts$colour_by,collapse=" "), outdir)
    
    # Run
    print(cmd)
    system(cmd)
  }
}

##################################################################
## Run all samples at once, Harmony batch correction per sample ##
##################################################################

# # LSF params
# opts$memory <- 15000

# outdir <- sprintf("%s/all_cells_harmony",io$outdir); dir.create(outdir, showWarnings = F)
# # opts$colour_by <- c("celltype.mapped","sample","stage","log_nFrags_atac","doublet_call")
# opts$colour_by <- c("celltype.predicted","stage")

# for (j in opts$nfeatures) {
#   for (k in opts$ndims) {
    
#     # Define LSF command
#     if (grepl("ricard",Sys.info()['nodename'])) {
#       lsf <- ""
#     } else if (grepl("ebi",Sys.info()['nodename'])) {
#       lsf <- sprintf("bsub -M %s -n 1 -o %s/allsamples_harmony_%d_%d.txt", opts$memory, io$tmpdir,j,k)
#     }
#     cmd <- sprintf("%s Rscript %s --samples %s --nfeatures %d --ndims %d  --matrix %s --batch.variable %s --batch.method Harmony --n_neighbors %s --min_dist %s --colour_by %s --outdir %s",
#                    lsf, io$script, paste(opts$samples,collapse=" "), j, k, opts$matrix, opts$batch.variable, paste(opts$n_neighbors,collapse=" "), paste(opts$min_dist,collapse=" "), paste(opts$colour_by,collapse=" "), outdir)
    
#     # Run
#     # print(cmd)
#     # system(cmd)
#   }
# }

##############################################################
## Run all samples at once, MNN batch correction per sample ##
##############################################################

# outdir <- sprintf("%s/all_cells_mnn",io$outdir); dir.create(outdir, showWarnings = F)
# # opts$colour_by <- c("celltype.mapped","sample","stage","log_nFrags_atac","doublet_call")
# opts$colour_by <- c("celltype.mapped","celltype.predicted","sample","stage","log_nFrags_atac")

# for (j in opts$nfeatures) {
#   for (k in opts$ndims) {

#     # Define LSF command
#     if (grepl("ricard",Sys.info()['nodename'])) {
#       lsf <- ""
#     } else if (grepl("ebi",Sys.info()['nodename'])) {
#       lsf <- sprintf("bsub -M %s -n 1 -o %s/allsamples_mnn_%d_%d.txt", opts$memory, io$tmpdir,j,k)
#     }
#     cmd <- sprintf("%s Rscript %s --samples %s --nfeatures %d --ndims %d --matrix %s --batch.variable %s --batch.method MNN --n_neighbors %s --min_dist %s --colour_by %s --outdir %s",
#                    lsf, io$script, paste(opts$samples,collapse=" "), j, k, opts$matrix, opts$batch.variable, paste(opts$n_neighbors,collapse=" "), paste(opts$min_dist,collapse=" "), paste(opts$colour_by,collapse=" "), outdir)

#     # Run
#     print(cmd)
#     system(cmd)
#   }
# }