library(edgeR)
library(biomaRt)
# DATASETS_PATH = paste("Datasets-",format(Sys.time(), "%Y%m%d/"), sep="")
RESULT_PATH = paste("Results-",format(Sys.time(), "%Y%m%d/"), sep="")
dir.create(RESULT_PATH)
for( region in c("VTA","CPU","PFC","HIP","NAC","BLA") ){
  # region = "VTA"
  # counts <- read.csv(paste0("../data/0921/0825_", region, "_sample_averaged_profile.csv") , header = TRUE, row.names = 1)
  counts <- read.csv(paste0("0826_", region, "_sample_maintained_profile_annot.csv") , header = TRUE, row.names = 1)
  # annoted_averaged_gex_lv-%s.csv
  # counts <- read.csv( paste("rescued_gex_lv-", region, ".csv", sep="") , header = TRUE, row.names = 1)
  # geneLengths <- read.csv("DATASETS/geneLengths.csv",header = TRUE)
  
  total_mapped = colSums(counts)
  # RPKM_counts = rpkm(counts, geneLengths$length, lib.size=total_mapped) # Get RPKM
  
  snames <- colnames(counts) # Sample names
  snames
  sim  <- substr(snames, 2, 3)
  time  <- substr(snames, 5, 8)
  group <- interaction(sim, time)
  group <- interaction(time)
  conditions = group
  # conditions = c("con", "trt")
  group = conditions
  # conditions <- factor(VTA_metatab$group)#factor out the group: it enables to find automatically which group belongs to
  test.expression = data.frame(matrix(NA, nrow=nrow(counts), ncol=0))
  
  # counts = counts[filtered.rownames,]
  counts[counts < 0] <- 0
  counts.mat <- as.matrix(counts)
  
  # design <- model.matrix(~0+group,VTA_metatab)
  design <- model.matrix(~0+group)
  #design
  
  #cont.matrix <- makeContrasts(deg=group99-group00, levels=design)
  contrasts = c()
  
  for (i in 1:99) {
    contrasts = append(contrasts, paste0(i, "-0"))
  }
  
  cont.matrix <- makeContrasts(contrasts=contrasts,
                               levels=design)
  
  voom<-voom(counts.mat,design,plot=TRUE)
  write.csv(voom, paste(RESULT_PATH, "rescued_voom_result_", region, ".csv",sep=""))
  fit<-lmFit(voom,design)
  vfit <- contrasts.fit(fit, cont.matrix)
  efit <- eBayes(vfit)
  plotSA(efit)
  
  coefficients=colnames(cont.matrix)
  summarytab <- data.frame()
  for (i in 1:length(coefficients)) {
    tab <- topTable(efit, number=length(row.names(counts.mat)), coef=i, adjust="BH")
    # tab <- merge(tab, gene_id_DB, by.x=0, by.y=0)
    # names(tab)[1]="gene_id"
    filename = paste(RESULT_PATH, region,"_rescued_",coefficients[i],".csv",sep="")
    write.table(tab,filename,quote=F,sep = "\t", row.names = T, col.names = T)
  }
  
  
}