#Loading pre-computed lineage
load("~/Documents/DriverPrediction/LineageMisc/Data/TCX_Monocle.RData")
source('MonoclePlottingCustom.R')
source('MiscPreprocessing.R')
source('LineageFunctions.R')



#Create gene clusters based on state expression patterns 
MonRun2 = MonRun
ScaledDat = ScaleMatrix(MonRun2@assayData$exprs)
MonRun$State2 <- MonRun$State
MonRun$State2[MonRun$State == 6] <- 4
MonRun$State2[MonRun$State == 7] <- 6

#getting gene symbol
gene_short_name <- Make.Gene.Symb(GeneNamesAD)

UnqStates <- unique(MonRun$State2)

cNames <- paste(UnqStates,'_mean',sep = '')

library(Matrix)

M <- matrix(rep(0,length(GeneNamesAD)*length(UnqStates)),nrow = length(GeneNamesAD),ncol = length(UnqStates))
colnames(M) <- cNames
rownames(M) <- gene_short_name


#get branch means 
for (i in 1:length(UnqStates)){
  
  for (j in 1:length(GeneNamesAD)){
    
    M[j,i] <- mean(ScaledDat[j,which(MonRun$State2==UnqStates[i])])
    
  }
  
}

M2 <- ScaleMatrix(M) #scaled matrix
M3 <- 1*(M2>0.5) #binarized scaled matrix 

library(pheatmap)

#performing biclustering
pheatmap(M2)
write.csv(M2df,'./Data/TCX_F_pv1_bclust_norm2.csv')




