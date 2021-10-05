#####################
## Define settings ##
#####################

# Load default settings
source(here::here("settings.R"))
source(here::here("utils.R"))

# Options

# I/O
io$outdir <- file.path(io$basedir,"results/rna/differential/pseudobulk")

##############################
## Load pseudobulk RNA data ##
##############################

rna.sce <- readRDS(io$rna.pseudobulk.sce)[,opts$celltypes]

rna.dt <- logcounts(rna.sce) %>%
  as.data.table(keep.rownames = T) %>%
  setnames("rn","gene") %>%
  melt(id.vars="gene", variable.name="celltype", value.name="expr")

#############################
## Differential expression ##
#############################

for (i in 1:length(opts$celltypes)) {
  for (j in i:length(opts$celltypes)) {
    if (i!=j) {
      
      rna_filt.dt <- rna.dt[celltype%in%c(opts$celltypes[[i]],opts$celltypes[[j]])] %>% 
        .[,expr:=round(expr,2)] %>%
        .[,logFC:=round(groupA-groupB,2)] %>%
        dcast(gene~celltype,value.var="expr") %>%
        setnames(c("gene","groupA","groupB"))
      
      # save      
      outfile <- sprintf("%s/%s_vs_%s.txt.gz", io$outdir,opts$celltypes[[i]],opts$celltypes[[j]])
      fwrite(rna_filt.dt, outfile, sep="\t")
    }
  }
}
