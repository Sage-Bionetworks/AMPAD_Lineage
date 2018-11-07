source('MiscPreprocessing.R')


Dat <- read.delim('Data/ROSMAP_DLPFC_logCPM.tsv',stringsAsFactors = F)
#Dat2 <- read.delim('Data/ROSMAP_DLPFC_Covariates.tsv',stringsAsFactors = F)
Dat2 <- read.csv('Data/ROSMAP_DLPFC_Covariates3.csv',stringsAsFactors = F)

#AMP_mods <-  read.csv('TCX_AMPAD_Modules.csv')
#AMP_mods <-  read.csv('TCX_DE.csv')
AMP_mods <-  read.csv('Data/DLPFC_DE.csv')
In <- which(AMP_mods$LogPV.DLPFC >= 1)
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

DatNorm2 <- DatNorm2[,Dat2$Batch<7]
Dat2 <- Dat2[Dat2$Batch<7,] 



DatNorm3 <- DatNorm2
Dat3 <- Dat2

#Keeping only female data 
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
#temp2 <- cbind(Dat4,Cov)
temp2 <- Dat4
temp2$APOE4 <- as.character(temp2$APOE4)
temp2$braaksc <- as.character(temp2$braaksc)
temp2$ceradsc <- as.character(temp2$ceradsc)
temp2$cogdx.1 <- as.character(temp2$cogdx.1)

gene_short_name <- Make.Gene.Symb(GeneNamesAD)

rownames(temp) <- NULL
colnames(temp) <- NULL
MonRun <- RunMonocleTobit(temp, temp2, C_by = 'Pseudotime',gene_short_name = gene_short_name)

plot_cell_trajectory(MonRun, color_by = "Diagnosis")


#Create differential expression based on state 
MonRun2 = MonRun
ScaledDat = ScaleMatrix(MonRun2@assayData$exprs)
#MonRun2@assayData$exprs <- ScaledDat
MonRun$State2 <- MonRun$State
MonRun$State2[MonRun$State == 4] <- 3
MonRun$State2[MonRun$State == 7] <- 3
MonRun$State2[MonRun$State == 6] <- 5


#find centers for each state 
UnqStates <-as.vector(unique(MonRun$State))
C_s1 <- rep(0,length(UnqStates))
C_s2 <- rep(0,length(UnqStates))

for (i in 1:length(UnqStates)){
  In_s <- which(MonRun$State==UnqStates[i])
  C_s1[i] <- mean(MonRun@reducedDimS[1,In_s])
  C_s2[i] <- mean(MonRun@reducedDimS[2,In_s])
  
}

D_mean <- rep(0,length(MonRun$State))

#For each sample, find the distance to the state mean 
for(i in 1:length(MonRun$State)){
  In_S <- which(UnqStates==MonRun$State[i])
  D_mean[i] <- sqrt((MonRun@reducedDimS[1,i]-C_s1[In_S])**2 + (MonRun@reducedDimS[2,i]-C_s2[In_S])**2)
  
}

#Normalize each sample's distance based on max distance for that state
D_meanNorm <- rep(0,length(MonRun$State))

for(i in 1:length(MonRun$State)){
  In_S <- which(MonRun$State==MonRun$State[i])
  D_meanNorm[i] <- (D_mean[i]-min(D_mean[In_S]))/max(D_mean[In_S])
  
}

#Make new dataframe with SampleID, D_mean, D_meanNorm, State, Pseudotime
l <- list()

MonRun$D_mean <- D_mean
MonRun$D_meanNorm <- D_meanNorm

l$SampleID <- MonRun$SampleID
l$D_mean <- MonRun$D_mean
l$D_meanNorm <- MonRun$D_meanNorm
l$State <- MonRun$State
l$Pseudotime <- MonRun$Pseudotime



l <- as.data.frame(l)

write.csv(l,'Data/DLPFC_F_pv1_StateDist.csv')