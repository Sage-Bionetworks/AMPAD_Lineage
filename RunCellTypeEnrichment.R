library(dplyr)
GeneNameFile <- 'GN_pv1.rda'
load(GeneNameFile)

load('E:/SageDocs/AMP-AD_project_stuff/CellTypeData/mouseMarkerGenes.rda')
TCX_genes <- mouseMarkerGenes$Cortex
for(i in 1:length(names(TCX_genes))){
  TCX_genes[[i]] <- toupper(TCX_genes[[i]])
}


library(R.matlab)
library(utilityFunctions)
ClusterList <- 'TCX_F_pv1_k15.mat'
RM <- readMat(ClusterList)

GN_hgnc <- convertEnsemblToHgnc(GN)


l <- list()



for (i in 1:max(RM$IDX)){
  
  In <- which(GN_hgnc$ensembl_gene_id %in% GN[RM$IDX ==i])   
  l[[i]] <- GN_hgnc$external_gene_name[In]

}

names(l) <- as.character(c(1:length(l)))

refGeneSet <- l %>%
  unlist %>%
  unique


res <- list()
res$fisher <- utilityFunctions::outerSapplyParallel(utilityFunctions::fisherWrapper,
                                                    l,
                                                    TCX_genes,
                                                    refGeneSet)
#pvalues are odd rows, odds ratios are even rows
res$pval <- res$fisher[which(1:nrow(res$fisher)%%2==1),]
rownames(res$pval) <- names(TCX_genes)
res$OR <- res$fisher[which(1:nrow(res$fisher)%%2==0),]
rownames(res$OR) <- names(TCX_genes)

#for each cluster identify most enriched cell type 

df <- list()
df$cluster <- names(l)
df$CellType <- names(l)
df$pValue <- names(l)
TCX_names <- names(TCX_genes)

for (i in 1:length(l)){
  tmpIn <- which(res$pval[,i] == min(res$pval[,i]))
  df$CellType[i] <- TCX_names[tmpIn]
  df$pValue[i] <- min(res$pval[,i])
}

df <- data.frame(df)

saveRDS(df,'TCX_k15_pv1.rds')


