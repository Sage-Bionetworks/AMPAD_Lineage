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
diff_test_res <- differentialGeneTest(MonRun,fullModelFormulaStr = "~State2")
d <- dim(diff_test_res)

diff_test_res$maxType <- rep(0,d[1])
diff_test_res$minType <- rep(0,d[1])

l <- as.array(unique(MonRun$State2))

for(i in 1:d[1]){
  
  t <- rep(0,length(l))
  
  for(j in 1:length(l)){
    
    t[j] <- mean(MonRun2@assayData$exprs[i,MonRun$State2==l[j]])
    
  }
  
  diff_test_res$maxType[i] <- which.max(t)
  diff_test_res$minType[i] <- which.min(t)
  
}

In_sort <- sort(diff_test_res$qval, method = 'radix', index.return = T, 
                na.last = T)

diff_test_res2 <- diff_test_res[In_sort$ix,]
saveRDS(diff_test_res2,'Data/DLPFC_F_pv1_Mon_State_DE.rds')


#Calculate Pseudotime correlation 
CorrMat <- rep(0,length(gene_short_name))

for( i in 1:length(gene_short_name)){
  
  CorrMat[i] <- abs(cor(x = MonRun$Pseudotime, y = as.vector(ScaledDat[i,])
                        , method = 'pearson'))[1]
  
}

CorrDF <- as.data.frame(CorrMat)
CorrDF$gene_short_name <- gene_short_name
saveRDS(CorrDF,'Data/DLPFC_F_pv1_PS_corr.rds')


#Pseudotime enrichment 
IGAP_list <- read.csv('Data/IGAP_gene_summary.csv')
PS <- seq(0,.8,.05)

MinEnr <- rep(0,length(PS))
MeanEnr <- rep(0,length(PS))

for(i in 1:length(PS)){
  
  CutOff <- PS[i]
  PS_genes <- CorrDF$gene_short_name[which(CorrMat>CutOff & CorrMat<CutOff+.1)]
  
  '%ni%' <- Negate('%in%')
  
  In_igap1 <- which(IGAP_list$Names %in% PS_genes)
  In_igap2 <- which(IGAP_list$Names %ni% PS_genes)
  
  MinEnr[i] <- t.test(IGAP_list$Mean[In_igap1],IGAP_list$Mean[In_igap2])$p.value
  MeanEnr[i] <- t.test(log10(IGAP_list$Min[In_igap1]),log10(IGAP_list$Min[In_igap2]))$p.value
  
  
}

plot(PS+.05,log10(MinEnr),type='l',col='Blue')
lines(PS+.05,log10(MeanEnr),type='l',col='Red')

p <- ggplot(MonRun@phenoData@data, aes(x=Diagnosis, y=Pseudotime)) + geom_violin()
p + stat_summary(fun.y=mean, geom="point", shape=23, size=2)



#DE of resistant branch 
BEAM_res <- BEAM(MonRun, branch_point = 2, cores = 4)
In_sort <- sort(BEAM_res$qval, method = 'radix', index.return = T, na.last = T)
BEAM_res <- BEAM_res[In_sort$ix,]
saveRDS(diff_test_res2,'Data/TCX_F_pv1_Mon_Branch2_DE.rds')