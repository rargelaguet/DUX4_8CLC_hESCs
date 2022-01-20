# Save matrix.mtx.gz
mtx <- getMatrixFromProject(ArchRProject.filt, useMatrix = "PeakMatrix", binarize = T)@assays@data[[1]]
outfile <- paste0(io$archR.directory,"/PeakMatrix/matrix.mtx.gz")
Matrix::writeMM(mtx, file=outfile)

# Save barcodes.tsv.gz
barcodes <- colnames(mtx)
outfile <- paste0(io$archR.directory,"/PeakMatrix/barcodes.tsv.gz")
write.table(barcodes, outfile, col.names = F, row.names = F, quote=F)

# save features.tsv.gz
features <- getPeakSet(ArchRProject.filt) %>% as.data.table() %>% setnames(c("seqnames"),c("chr")) %>%
  .[,c("chr","start","end")] %>% 
  .[,peak_idx:=sprintf("%s:%s-%s",chr,start,end)] %>% 
  .[,peak_idx]
outfile <- paste0(io$archR.directory,"/PeakMatrix/features.tsv.gz")
write.table(features, outfile, col.names = F, row.names = F, quote=F)
