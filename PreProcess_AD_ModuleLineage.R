setwd('E:/SageDocs/PredictingDriverGenes/LineageMisc/')

Dat <- read.delim('MAYO_CBE_TCX_logCPM.tsv',stringsAsFactors = F)
Dat2 <- read.delim('MAYO_CBE_TCX_Covariates.tsv',stringsAsFactors = F)

AMP_mods <-  read.csv('TCX_AMPAD_Modules.csv')

GeneNames <- Dat$ensembl_gene_id
GeneNamesAD <- AMP_mods$GeneID

Names <- colnames(Dat)

for (i in 1:length(Names)){
  
  Names[i] <- substring(Names[i],2)
  
}

colnames(Dat) <- Names
cNames <- Dat2$SampleID
l <- length(Names)

#deleting columns not in the covariate list 
for (i in 1:l){
  if (!(Names[i] %in% cNames)){
    Dat[,i] <- NULL
  }
}




#Normalize all columns 
source('MiscPreprocessing.R')

DatNorm <- ColNorm(Dat)
In <- which(GeneNames %in% GeneNamesAD)
DatNorm2 <- DatNorm[Get.High.Var.Genes(DatNorm),]
library(Rtsne)
Temp <- Rtsne(t(DatNorm2))

plot(Temp$Y[,1],Temp$Y[,2], col = as.factor(Dat2$Sex))

#Keeping only TCX data 
In_BR <- grep('TCX',Dat2$Tissue.Diagnosis)
DatNorm3 <- DatNorm2[,In_BR]
Dat3 <- Dat2[In_BR,]
Temp <- Rtsne(t(DatNorm3))
plot(Temp$Y[,1],Temp$Y[,2], col = as.factor(Dat3$Sex))

#Keeping only female data 
In_S <- grep('FEMALE',Dat3$Sex)
DatNorm4 <- DatNorm3[,In_S]
Dat4 <- Dat3[In_S,]
Temp <- Rtsne(t(DatNorm4))
plot(Temp$Y[,1],Temp$Y[,2], col = as.factor(Dat4$Tissue.Diagnosis))

#Performing lineage inference with Monocle2
source('LineageFunctions.R')
temp <- exp(DatNorm4)-1
rownames(temp) <- NULL
colnames(temp) <- NULL
MonRun <- RunMonocleDefault(temp, Dat4$AgeAtDeath)
plot_cell_trajectory(MonRun, color_by = "Labels")