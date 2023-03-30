library(WGCNA)
library(ggplot2)
library(RFLPtools)
library(ape)
library(org.Mm.eg.db)
library(AnnotationDbi)
library(GSEABase)
library(GOstats)
library(plyr)
library(reshape2)
require(phytools)
library(hash)
library(png)
library(jpeg)

DATA_PATH = "Dataset/"
RESULT_PATH = "Results/temporal/"
enableWGCNAThreads(nThreads = 16)

### Estimate power value for each region
for ( region in c("PFC", "NAC", "VTA", "BLA","CPU", "HIP") ) {

  WGCNA_matrix = read.csv(paste0(DATA_PATH, region, "_sample_averaged_profile.csv"), row.names = "Gene")
  WGCNA_matrix = t(WGCNA_matrix[])
  s = abs(bicor(WGCNA_matrix, nThreads = 16))
  powers = c(c(1:10), seq(from = 12, to=300, by=10))
  sft = pickSoftThreshold(WGCNA_matrix, powerVector = powers, verbose = 5)


  plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
       xlab='Soft Threshold (power)',ylab='Scale Free Topology Model Fit,signed R^2',
       type='n', main = paste0('Scale independence of ', region));
  text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
       labels=powers,cex=1,col='red'); abline(h=0.90,col='red')
  # pdf(file= paste0(RESULT_PATH, region, "-Soft Threshold plot.pdf"))
  # dev.off()
  save.image(file = paste0(RESULT_PATH, region, "-work_space_full-part1.RData"))

}

### Selected power values for each region
sft_hash = hash("PFC"=155, "NAC"=101,"VTA"=43,"BLA"=89, "CPU"=59, "HIP"=129)

### Plotting power value selection - focused on Selected power values
for ( region in c("PFC", "NAC", "VTA", "BLA","CPU", "HIP") ) {
  
  WGCNA_matrix = read.csv(paste0(DATA_PATH, "/", region, "_sample_averaged_profile.csv"), row.names = "Gene")
  WGCNA_matrix = t(WGCNA_matrix[])
  powers = c(1:(sft_hash[[region]] + 30))
  # powers = c(c(1:10), seq(from = sft_hash[[region]] - 30, to=sft_hash[[region]] + 30, by=1))
  sft = pickSoftThreshold(WGCNA_matrix, powerVector = powers, verbose = 5)
  
  plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
       xlab='Soft Threshold (power)',ylab='Scale Free Topology Model Fit,signed R^2',
       type='n', main = paste0('Scale independence of ', region));
  text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
       labels=powers,cex=1,col='red'); abline(h=0.90,col='red')
  write.csv(sft$fitIndices, file= paste0("20210610_", region, "_sft.csv"))
  
}

### Module detection based on temporal profile
for ( region in c("PFC", "NAC", "VTA", "BLA","CPU", "HIP") ) {

  WGCNA_matrix = read.csv(paste0(DATA_PATH, DATA_DATE, "_", region, "_sample_averaged_profile.csv"), row.names = "Gene")
  WGCNA_matrix = t(WGCNA_matrix[])

  softPower = sft_hash[[region]]
  adjacency = adjacency(WGCNA_matrix, power = softPower, type = "signed")
  TOM = TOMsimilarity(adjacency, TOMType="signed") # specify network type
  dissTOM = 1-TOM

  geneTree = hclust(as.dist(dissTOM), method="average")
  plot(geneTree, xlab="", sub="", main= paste0("Gene Clustering on TOM-based dissimilarity-", region), labels= FALSE, hang=0.04)

  #This sets the minimum number of genes to cluster into a module
  minModuleSize = 30
  modules = cutreeDynamic(dendro= geneTree, distM= dissTOM, deepSplit=2, pamRespectsDendro= FALSE, minClusterSize = minModuleSize)
  module.colours = labels2colors(modules)

  plotDendroAndColors(geneTree, module.colours, c("Dynamic Tree Cut", "Merged dynamic"),
                      dendroLabels= FALSE, hang=0.03, addGuide= TRUE, guideHang=0.05,
                      main =paste0('Dendrogram with modules of nodes in colors-', region, "(", softPower, ")")
                      )
  save.image(file = paste0(RESULT_PATH, region, "-work_space_dendro.RData"))

}

### GO for each module and calculate corr
for ( region in c("PFC", "NAC", "VTA", "BLA","CPU", "HIP") ) {


  load(paste0(RESULT_PATH, region, "-work_space_dendro.RData"))

  ref_genes = colnames(WGCNA_matrix)

  GO = toTable(org.Mm.egGO)
  SYMBOL = toTable(org.Mm.egENSEMBL)
  GO_data_frame = data.frame(GO$go_id, GO$Evidence,SYMBOL$ensembl_id[match(GO$gene_id,SYMBOL$gene_id)])
  GO_data_frame = na.omit(GO_data_frame)
  GO_ALLFrame = GOAllFrame(GOFrame(GO_data_frame, organism = 'Mouse mucus'))

  gsc <- GeneSetCollection(GO_ALLFrame, setType = GOCollection())

  GSEAGO = vector('list',length(unique(modules)))
  print(length(unique(module.colours))-1)
  # for(i in 0:(length(unique(modules))-1)){
    # GSEAGO[[i+1]] = summary(hyperGTest(GSEAGOHyperGParams(name = 'Mouse mucus GO',
    #                                                       geneSetCollection = gsc, geneIds = colnames(WGCNA_matrix)[modules==i+1],
    #                                                       universeGeneIds = ref_genes, ontology = 'BP', pvalueCutoff = 0.05,
    #                                                       conditional = FALSE, testDirection = 'over')))
  #   print(i)
  # }
  for(i in 0:(length(unique(modules))-1)){
    GSEAGO[[i+1]] = summary(hyperGTest(GSEAGOHyperGParams(name = 'Homo sapiens GO',
                                                          geneSetCollection = gsc, geneIds = colnames(WGCNA_matrix)[modules==i],
                                                          universeGeneIds = ref_genes, ontology = 'BP', pvalueCutoff = 0.05,
                                                          conditional = FALSE, testDirection = 'over')))
    print(i)
  }

  # smaller ontology size; more specific ontology term
  cutoff_size = 50
  GO_module_name = rep(NA,length(unique(modules)))

  for (i in 1:length(unique(modules))){
      print(i)
      GO_module_name[i] =
      GSEAGO[[i]][GSEAGO[[i]]$Size<cutoff_size,
                  ][which(GSEAGO[[i]][GSEAGO[[i]]$Size<cutoff_size,]$Count==max(GSEAGO[[i]][GSEAGO[[i]]$
                                                                                              Size<cutoff_size,]$Count)),7]
  }

  # Correlate traits --------------------------------------------------------

  datTraits = read.csv("Sample_trait.csv")
  head(datTraits)
  #form a data frame analogous to expression data that will hold the clinical traits.
  rownames(datTraits) = datTraits$Sample
  datTraits$Sample = NULL
  table(rownames(datTraits)==rownames(WGCNA_matrix)) #should return TRUE if datasets align correctly, otherwise your names are out of order
  head(datTraits)

  #Define number of genes and samples
  nGenes = ncol(WGCNA_matrix)
  nSamples = nrow(WGCNA_matrix)

  #Recalculate MEs with color labels
  MEs0 = moduleEigengenes(WGCNA_matrix,  module.colours, verbose=5)$eigengenes
  MEs = orderMEs(MEs0)
  moduleTraitCor = cor(MEs, datTraits, use= "p")
  moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples)

  textMatrix= paste(signif(moduleTraitCor, 2), "\n(",
                    signif(moduleTraitPvalue, 1), ")", sep= "")
  dim(textMatrix)= dim(moduleTraitCor)

  highest_corr_list = c()
  for (index in rownames(moduleTraitCor)) {
    highest_corr_list = append(highest_corr_list, max(moduleTraitCor[index,]))

  }
  max(highest_corr_list)

  lst <- sort(highest_corr_list, index.return=TRUE, decreasing=TRUE)
  top_5_corr_GO_indexs = lapply(lst, `[`, lst$x %in% head(unique(lst$x),5))


  highest_corr_list_abs = c()
  for (index in rownames(moduleTraitCor)) {
    highest_corr_list_abs = append(highest_corr_list_abs, sum(abs(moduleTraitCor[index,])))
  }
  max(highest_corr_list_abs)

  lst <- sort(highest_corr_list_abs, index.return=TRUE, decreasing=TRUE)
  top_5_corr_abs_GO_indexs = lapply(lst, `[`, lst$x %in% head(unique(lst$x),5))

  summerized_GO = append(top_5_corr_GO_indexs$ix, top_5_corr_abs_GO_indexs$ix)

  empty_x_arr = c()

  for (x in names(datTraits)) {
    empty_x_arr = append(empty_x_arr, " ")
  }

  Heatmap_GO_label = match(colnames(MEs), paste0('ME',unique(module.colours)))
  Heatmap_GO_label = GO_module_name[Heatmap_GO_label]

  write.csv(MEs0, file= paste0(RESULT_PATH, region, "-MEs0-tmp.csv"))
  write.csv(module.colours, file= paste0(RESULT_PATH, region, "-module_colours-tmp.csv"))
  write.csv(GO_module_name, file= paste0(RESULT_PATH, region, "-GO_term-tmp.csv"))
  write.csv(moduleTraitCor, file= paste0(RESULT_PATH, region, "-moduleTraitCor-tmp.csv"))
  write.csv(moduleTraitPvalue, file= paste0(RESULT_PATH, region, "-moduleTraitPvalue-tmp.csv"))
  save.image(file = paste0(RESULT_PATH, "-work_space_corr.RData"))

}

### Saving missing values

for ( region in c("PFC", "VTA", "NAC","CPU", "HIP", "BLA") ) {

load(paste0(RESULT_PATH, region, "-work_space_corr.RData"))
MM = as.data.frame(WGCNA_matrix, MEs, use = "p")
top_H_Gene = chooseTopHubInEachModule(WGCNA_matrix, module.colours)
datKME = signedKME(WGCNA_matrix, MEs0, outputColumnName = "MM.")

write.csv(MM, file= paste0(RESULT_PATH, region, "-MM.csv"))
write.csv(top_H_Gene, file= paste0(RESULT_PATH, region, "-topHG.csv"))
write.csv(datKME, file= paste0(RESULT_PATH, region, "-KME-tmp.csv"))

}

### Saving missing values = MEDiss

for ( region in c("PFC", "VTA", "NAC","CPU", "HIP", "BLA") ) {


  load(paste0(RESULT_PATH, region, "-work_space_corr.RData"))
  Heatmap_GO_label = match(colnames(MEs0), paste0('ME',unique(module.colours)))
  Heatmap_GO_label = GO_module_name[Heatmap_GO_label]
  colnames(MEs0) = Heatmap_GO_label
  # Calculate dissimilarity of module eigengenes
  MEDiss = 1-cor(MEs0);
  # Cluster module eigengenes
  METree = hclust(as.dist(MEDiss), method = "average");
  # Plot the result
  # sizeGrWindow(7, 6)
  plot(METree, main = paste0("Dendrogram of module eigengenes in ", region), xlab = "", sub = "")

  write.csv(MEDiss, file= paste0(RESULT_PATH, region, "-MEDiss-tmp.csv"))


}

### Network Heatmap plot

for ( region in c("PFC", "VTA", "NAC","CPU", "HIP", "BLA") ) {
  load(paste0(RESULT_PATH, region, "-work_space_corr.RData"))
  heatmap(adjacency[1:100, 1:100])
  TOMplot(
    dissTOM,
    geneTree,
    module.colours,
    main = paste0("Network Heatmap plot,", region)
  )
}
