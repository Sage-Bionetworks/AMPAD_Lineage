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


p <- ggplot(MonRun@phenoData@data, aes(x=Diagnosis, y=Pseudotime,fill=Diagnosis)) + geom_boxplot()
p + stat_summary(fun.y=mean, geom="point", shape=23, size=2)  + geom_jitter(width = 0.2)

p <- ggplot(MonRun@phenoData@data, aes(x=APOE4, y=Pseudotime,fill=APOE4)) + geom_boxplot()
p + stat_summary(fun.y=mean, geom="point", shape=23, size=2)  + geom_jitter(width = 0.2)

p <- ggplot(MonRun@phenoData@data, aes(x=cogdx.1, y=Pseudotime,fill=cogdx.1)) + geom_boxplot()
p + stat_summary(fun.y=mean, geom="point", shape=23, size=2)  + geom_jitter(width = 0.2)

p <- ggplot(MonRun@phenoData@data, aes(x=braaksc, y=Pseudotime,fill=braaksc)) + geom_boxplot()
p + stat_summary(fun.y=mean, geom="point", shape=23, size=2)  + geom_jitter(width = 0.2)

p <- ggplot(MonRun@phenoData@data, aes(x=ceradsc, y=Pseudotime,fill=ceradsc)) + geom_boxplot()
p + stat_summary(fun.y=mean, geom="point", shape=23, size=2)  + geom_jitter(width = 0.2)

#Create differential expression based on state 
MonRun2 = MonRun
ScaledDat = ScaleMatrix(MonRun2@assayData$exprs)
#MonRun2@assayData$exprs <- ScaledDat
MonRun$State2 <- MonRun$State
#MonRun$State2[MonRun$State == 4] <- 3
#MonRun$State2[MonRun$State == 7] <- 3
#MonRun$State2[MonRun$State == 6] <- 5
MonRun$State2[MonRun$State == 4] <- 3
MonRun$State2[MonRun$State == 7] <- 3
MonRun$State2[MonRun$State == 6] <- 4
MonRun$State2[MonRun$State == 5] <- 4
MonRun$State2[MonRun$State == 8] <- 5
MonRun$State2[MonRun$State == 9] <- 6



l2 <- list()
l2$gene_names <- gene_short_name
l2$p_2 <- rep(0,length(gene_short_name))
l2$p_3 <- rep(0,length(gene_short_name))
l2$p_4 <- rep(0,length(gene_short_name))
l2$p_5 <- rep(0,length(gene_short_name))
l2$p_6 <- rep(0,length(gene_short_name))

l2$d_2 <- rep(0,length(gene_short_name))
l2$d_3 <- rep(0,length(gene_short_name))
l2$d_4 <- rep(0,length(gene_short_name))
l2$d_5 <- rep(0,length(gene_short_name))
l2$d_6 <- rep(0,length(gene_short_name))


for (i in 1:length(gene_short_name)){
  l <- list()
  l$x <- as.vector(t(temp[i,]))
  l$s <- as.character(MonRun$State2)
  
  df <- as.data.frame(l)
  res.aov <- aov(x ~ s, data = df)
  tk <- TukeyHSD(res.aov)
  
  l2$p_2[i] <- tk$s[1,4]
  l2$p_3[i] <- tk$s[2,4]
  l2$p_4[i] <- tk$s[3,4]
  l2$p_5[i] <- tk$s[4,4]
  l2$p_6[i] <- tk$s[5,4]
  
  l2$d_2[i] <- tk$s[1,1]
  l2$d_3[i] <- tk$s[2,1]
  l2$d_4[i] <- tk$s[3,1]
  l2$d_5[i] <- tk$s[4,1]
  l2$d_6[i] <- tk$s[5,1]
  
}


df2 <- as.data.frame(l2)
saveRDS(df2,'Data/DLPFC_F_pv1_Anova_DE.rds')



