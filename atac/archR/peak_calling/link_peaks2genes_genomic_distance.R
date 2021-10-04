
#####################
## Define settings ##
#####################

if (grepl("ricard",Sys.info()['nodename'])) {
  source("/Users/ricard/gastrulation_multiome_10x/settings.R")
  source("/Users/ricard/gastrulation_multiome_10x/utils.R")
} else if (grepl("ebi",Sys.info()['nodename'])) {
  source("/homes/ricard/gastrulation_multiome_10x/settings.R")
  source("/homes/ricard/gastrulation_multiome_10x/utils.R")
} else {
  stop("Computer not recognised")
}

# I/O
io$outdir <- paste0(io$basedir,"/results/atac/archR/peak_calling/peaks2genes")

# Options
opts <- list()
opts$gene_window <- 1e5  # maximum window length for the overlap

###############
## Load data ##
###############

# Load gene metadata
gene_metadata <- fread(io$gene_metadata) %>% 
  .[,chr:=as.factor(sub("chr","",chr))] %>%
  setnames("symbol","gene") %>%
  .[, c("chr","start","end","gene","ens_id","strand")]

# Load peak metadata
peakSet.dt <- fread(io$archR.peak.metadata) %>%
  .[,chr:=as.factor(sub("chr","",chr))] %>%
  .[,c("chr","start","end")] %>%
  .[,peak:=sprintf("%s_%s_%s",chr,start,end)] %>%
  setkey(chr,start,end)


#############
## Overlap ##
#############

gene_metadata.ov <- gene_metadata

gene_metadata.ov[strand=="+",c("gene.start","gene.end"):=list(start,end)]
gene_metadata.ov[strand=="-",c("gene.start","gene.end"):=list(end,start)]

gene_metadata.ov[strand=="+",c("start","end"):=list (gene.start-opts$gene_window, gene.end+opts$gene_window)]
gene_metadata.ov[strand=="-",c("end","start"):=list (gene.start+opts$gene_window, gene.end-opts$gene_window)]

gene_metadata.ov %>% .[,strand:=NULL] %>% setkey(chr,start,end)

ov <- foverlaps(
  peakSet.dt,
  gene_metadata.ov,
  nomatch = NA
) %>%  .[,c("start","end"):=NULL] %>%
  setnames(c("i.start","i.end"),c("peak.start","peak.end")) %>%
  # setnames(c("start","end"),c("gene.start","gene.end")) %>%
  .[,peak.mean:=(peak.start+peak.end)/2] %>%
  # .[,c("start_dist","end_dist"):=list( abs(gene.end-peak.mean), abs(gene.start-peak.mean))] %>%
  # .[,dist:=abs(gene.end-peak.mean)] %>%
  .[,dist:=max(abs(gene.end-peak.mean), abs(gene.start-peak.mean)),by=c("gene","ens_id","peak")] %>%
  .[peak.mean>gene.start & peak.mean<gene.end,dist:=0]
  
  # .[,c("start_dist","end_dist"):=list( gene.end-peak.start, gene.start-peak.end)] %>%
  # .[,c("start_dist","end_dist"):=list( ifelse(end_dist<0 & start_dist>0,0,start_dist), ifelse(end_dist<0 & start_dist>0,0,end_dist) )] %>%
  # .[,dist:=ifelse(abs(start_dist)<abs(end_dist),abs(start_dist),abs(end_dist))] %>% .[,c("start_dist","end_dist"):=NULL]

# ov[peak=="1_170646867_170647467"]

# Select nearest gene
ov_nearest <- ov %>%
  .[.[,.I[dist==min(dist)], by=c("peak")]$V1] %>%
  .[complete.cases(.)] %>%
  .[!duplicated(peak)]

# Sanity check  
# ov_nearest$gene[(duplicated(ov_nearest$peak))]

##########
## Save ##
##########

fwrite(ov, paste0(io$outdir,"/peaks2genes_all.txt.gz"), sep="\t", na="NA")
fwrite(ov_nearest, paste0(io$outdir,"/peaks2genes_nearest.txt.gz"), sep="\t", na="NA")
