source('MiscPreprocessing.R')
source('LineageFunctions.R')


load('Data/DLPFC_Monocle.RData')
miDat <- read.table('../DataFiles/ROSMAP_arraymiRNA.gct',skip = 2, stringsAsFactors = F, header =T)

miNames <- miDat$Genes

Dat <- miDat[,3:dim(miDat)[2]]

In <- which(colnames(Dat) %in% MonRun$mirna_id)

Dat <- Dat[,In]

Diag <- c()
BRAAK <- c()
CERAD <- c()
CogDX <- c()
PST <- c()
APOE4 <- c()

for(i in 1:dim(Dat)[2]){
  
  In <- which(MonRun$mirna_id == colnames(Dat)[i])[1]
  Diag <- c(Diag,MonRun$Diagnosis[In])
  BRAAK <- c(BRAAK,MonRun$braaksc[In])
  CERAD <- c(CERAD,MonRun$ceradsc[In])
  CogDX <- c(CogDX,MonRun$cogdx.1[In])
  PST <- c(PST,MonRun$Pseudotime[In])
  APOE4 <- c(APOE4,MonRun$APOE4[In])
  
  
}

l <- list()

l$mirna_id <- colnames(Dat)
l$Diag <- Diag
l$BRAAK <- BRAAK
l$CERAD <- CERAD
l$CogDX <- CogDX
l$PST <- PST 
l$APOE4 <- APOE4

df <- as.data.frame(l)

DatNorm <- ColNorm(Dat)
In_genes <- Get.High.Var.Genes(DatNorm, frac = 0.3)
DatNorm <- DatNorm[In_genes,]
miNames <- miNames[In_genes]

temp <- DatNorm
temp2 <- df
rownames(temp) <- NULL
colnames(temp) <- NULL

MonRun2 <- RunMonocleTobit(temp, temp2, C_by = 'Pseudotime')
plot_cell_trajectory(MonRun2, color_by = "Diag")

#Calculate correlation with RNA-Seq pseudotime  
PST_corr <- rep(0,dim(Dat)[1])

for(i in 1:dim(Dat)[1]){
  
  PST_corr[i] <- (cor(x = MonRun2$PST, y = as.numeric(DatNorm[i,])
                      , method = 'spearman'))[1]
  
}

In <- which(abs(PST_corr)>0.4)
miDat$Genes[In]










