source('MiscPreprocessing.R')



Dat <- read.delim('Data/MAYO_CBE_TCX_logCPM.tsv',stringsAsFactors = F)
Dat2 <- read.delim('Data/MAYO_CBE_TCX_Covariates.tsv',stringsAsFactors = F)

Cov <- read.csv('Data/mayo_igap_snps.csv',stringsAsFactors = F)
Cov[,2:22] <- round(Cov[,2:22])

#AMP_mods <-  read.csv('Data/TCX_AMPAD_Modules.csv')
AMP_mods <-  read.csv('Data/TCX_DE.csv')
#AMP_mods <-  read.csv('DLPFC_DE.csv')
In <- which(AMP_mods$logPV >= 1)
AMP_mods <- AMP_mods[In,]




#Normalize all columns 

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
temp <- rep(T,l)
for (i in 1:l){
  if (!(Names[i] %in% cNames)){
    temp[i] <- F
  }
}

In <- which(temp)
#print(temp)
Dat <- Dat[,In]

#deleting extra rows in covariate list
Names <- Names[In]
l <- length(cNames)
temp <- rep(T,l)
for (i in 1:l){
  if (!(cNames[i] %in% Names)){
    temp[i] <- F
  }
}
In <- which(temp)
Dat2 <- Dat2[In,]


DatNorm <- ColNorm(Dat)
In_genes <- which(GeneNames %in% GeneNamesAD)
#In_genes <- Get.High.Var.Genes(Dat)
DatNorm2 <- DatNorm[In_genes,]
GeneNamesAD <- GeneNames[In_genes]


In_BR <- grep('TCX',Dat2$Tissue.Diagnosis)
#In_BR <- grep('DLPFC',Dat2$Tissue.Diagnosis)
DatNorm3 <- DatNorm2[,In_BR]
Dat3 <- Dat2[In_BR,]

Sex <- 'FEMALE'
In_S <- which(Dat3$Sex == Sex)
DatNorm4 <- DatNorm3[,In_S]
Dat4 <- Dat3[In_S,]

In_cov <- which(Cov$ID %in% Dat4$Donor_ID)
Cov <- Cov[In_cov,]
In_cov <- c()
for(i in 1:length(Dat4$Donor_ID)){
  temp <- which(Cov$ID == Dat4$Donor_ID[i])
  In_cov <- c(In_cov,temp[1])
}
Cov <- Cov[In_cov,]

for (i in 23:26){
  Cov[,i] <- (Cov[,i] - min(Cov[,i]))/(max(Cov[,i])-min(Cov[,i]))
}

source('LineageFunctions.R')
temp <- DatNorm4
temp2 <- cbind(Dat4,Cov)

gene_short_name <- Make.Gene.Symb(GeneNamesAD)

rownames(temp) <- NULL
colnames(temp) <- NULL
MonRun <- RunMonocleTobit(temp, temp2, C_by = 'Pseudotime',gene_short_name = gene_short_name)


plot_cell_trajectory(MonRun, color_by = "Tissue.Diagnosis")


#Create gene clusters based on state expression patterns 
MonRun2 = MonRun
ScaledDat = ScaleMatrix(MonRun2@assayData$exprs)
#MonRun2@assayData$exprs <- ScaledDat
MonRun$State2 <- MonRun$State
MonRun$State2[MonRun$State == 6] <- 4

#Loading Viral dataset 
load('Data/MAYO_TCX_RNA_workspace.RData')

# Identify samples that are in the not in the viral samples and the samples that are in it. 
DatV <- MAYO_RNA_workspace$virusLevelCounts
In_f <- which(colnames(MAYO_RNA_workspace$virusLevelCounts) %in% Dat4$SampleID)
In_nf <- which(!(Dat4$SampleID %in% colnames(MAYO_RNA_workspace$virusLevelCounts)[In_f]))

DatV <- DatV[,In_f]

#fill missing samples with dummy (0) values
for (i in 1:length(In_nf)){
  DatV[[Dat4$SampleID[In_nf[i]]]] <- rep(0,dim(DatV)[1])
}

#re-arrange columns to match MonRun 
In_ord <- c()

for(i in 1:length(Dat4$SampleID)){
  t <- which(colnames(DatV) == Dat4$SampleID[i])
  In_ord <- c(In_ord,t)
} 

DatV <- DatV[,In_ord]

DatV2 <- as.data.frame(t(DatV))
In_V <- which(colSums(DatV2)>=5)
DatV2 <- DatV2[,In_V]

cNamesV <- colnames(DatV2)
cNamesV2 <- paste(c(1:length(cNamesV)),'V',sep='')

for(i in 1:length(cNamesV)){
  MonRun[[cNamesV2[i]]] <- DatV2[[cNamesV[i]]]
}


V_corr <- rep(0,length(cNamesV))
#check the correlation of different viruses with pseudotime
for(i in 1:length(cNamesV)){
  V_corr[i] <- cor(MonRun$Pseudotime,DatV2[[cNamesV[i]]],method = 'pearson')[1]
}


