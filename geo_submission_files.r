library(lumi)
library(limma)
library(beadarray)

### Generate raw and normalized data to submit to geo

setwd('Z:/SERVICES/Microarrays/GX-Nugen/Angela Riedel/Sept2014/AR_Sep14_Results/')
raw_data<-lumiR("raw_data.txt")
norm_data<-lumiN(raw_data)

detection_pvals_raw_data=detection(raw_data)
detection_pvals_norm_data=detection(norm_data)

raw_table=as.data.frame(matrix(0,nrow=nrow(raw_data),ncol=(ncol(raw_data)*2)))
nc=1
for (i in 1:length(colnames(raw_data))){
  sample=colnames(raw_data)[i]
  raw_data_col=which(colnames(exprs(raw_data))==sample)
  raw_table[,nc]=exprs(raw_data)[,raw_data_col]
  colnames(raw_table)[nc]=sample
  detect_pval_index=which(colnames(detection_pvals_raw_data)==sample)
  raw_table[,nc+1]=detection_pvals_raw_data[,raw_data_col]
  colnames(raw_table)[nc+1]='Detection Pval'
  nc=2*i+1
  print(i)
  print(nc)
}
raw_table=cbind(ID_REF=rownames(exprs(raw_data)), raw_table)

norm_table=as.data.frame(matrix(0,nrow=nrow(norm_data),ncol=(ncol(norm_data)*2)))
nc=1
for (i in 1:length(colnames(norm_data))){
  sample=colnames(norm_data)[i]
  norm_data_col=which(colnames(exprs(norm_data))==sample)
  norm_table[,nc]=exprs(norm_data)[,norm_data_col]
  colnames(norm_table)[nc]=sample
  detect_pval_index=which(colnames(detection_pvals_norm_data)==sample)
  norm_table[,nc+1]=detection_pvals_norm_data[,norm_data_col]
  colnames(norm_table)[nc+1]='Detection Pval'
  nc=2*i+1
  print(i)
  print(nc)
}

norm_table=cbind(ID_REF=rownames(exprs(norm_data)), norm_table)

dir.create("geo_submission_files", showWarnings=F)
setwd('geo_submission_files/')
write.csv(raw_table, 'matrix_non_normalized.csv', row.names=FALSE)
write.csv(norm_table, 'matrix_normalized.csv', row.names=FALSE)
setwd("..")

## Geo files
dir.create("geo_files", showWarnings=F)
setwd('geo_files/')
produceGEOPlatformFile(x.lumi = raw_data)
produceGEOSampleInfoTemplate(lumiNormalized = norm_data)
produceGEOSubmissionFile(lumiNormalized = norm_data, lumiRaw = raw_data, sampleInfo = 'GEOsampleInfo.txt')
