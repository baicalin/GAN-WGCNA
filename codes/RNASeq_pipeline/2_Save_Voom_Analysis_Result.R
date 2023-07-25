if(0)"
if (!require('BiocManager', quietly = TRUE))
    install.packages('BiocManager')

BiocManager::install('edgeR')
"

install.packages('BiocManager')
BiocManager::install('edgeR')

library(edgeR)
library(biomaRt)
DATASETS_PATH = "Datasets/" # paste("Datasets-",format(Sys.time(), "%Y%m%d/"), sep="")
RESULT_PATH = "Results/" # paste("Results-",format(Sys.time(), "%Y%m%d/"), sep="")

dir.create(RESULT_PATH)
for( region in c("VTA","CPU","PFC","HIP","NAC","BLA") ){
    # region = "VTA"
    counts <- read.csv( paste(DATASETS_PATH, region, "_Count_Dataset.csv", sep="") , header = TRUE, row.names = 1)
    geneLengths <- read.csv("Datasets/geneLengths.csv",header = TRUE)

    total_mapped = colSums(counts)
    RPKM_counts = rpkm(counts, geneLengths$length, lib.size=total_mapped) # Get RPKM

    snames <- colnames(counts) # Sample names
    snames
    sa  <- substr(snames, 10, 10)
    ch  <- substr(snames, 11, 11)
    group <- interaction(sa, ch)
    conditions = group

    # conditions <- factor(VTA_metatab$group)#factor out the group: it enables to find automatically which group belongs to
    test.expression = data.frame(matrix(NA, nrow=nrow(RPKM_counts), ncol=0))

    for (i in levels(conditions)) {
        RPKM_counts_subset <- RPKM_counts[,conditions==i]
        count.by.min <- RPKM_counts_subset >= 1
        binary.mat <- as.data.frame(ifelse(count.by.min == "TRUE",1,0))

        cre.threshold <- floor((80/100)*ncol(RPKM_counts_subset))

        filtered = rowSums(binary.mat) >= cre.threshold
        test.expression[[i]] <- filtered
    }
    final.filter <- apply(test.expression,1,any)
    filtered.rpkm <- RPKM_counts[final.filter,]

    filtered.rownames <- row.names(filtered.rpkm)
    counts = counts[filtered.rownames,]
    counts.mat <- as.matrix(counts)

    # design <- model.matrix(~0+group,VTA_metatab)
    design <- model.matrix(~0+group)
    #design

    cont.matrix <- makeContrasts(CNvSN=groupC.N-groupS.N,
                                 SCvSS=groupS.C-groupS.S,
                                 CSvSS=groupC.S-groupS.S,
                                 CCvSS=groupC.C-groupS.S,
                                 levels=design)

    voom<-voom(counts.mat,design,plot=TRUE)
    write.csv(voom, paste(RESULT_PATH, "Voom_Result_", region, ".csv",sep=""))
    fit<-lmFit(voom,design)
    vfit <- contrasts.fit(fit, cont.matrix)
    efit <- eBayes(vfit)
    plotSA(efit)
    # For voom-limma
    # coefficients=colnames(cont.matrix)
    # summarytab <- data.frame()
    # for (i in 1:length(coefficients)) {
    #     tab <- topTable(efit, number=length(row.names(counts.mat)), coef=i, adjust="BH")
    #     # tab <- merge(tab, gene_id_DB, by.x=0, by.y=0)
    #     # names(tab)[1]="gene_id"
    #     filename = paste(RESULT_PATH, region,"_",coefficients[i],".csv",sep="")
    #     write.table(tab,filename,quote=F,sep = "\t", row.names = T, col.names = T)
}

