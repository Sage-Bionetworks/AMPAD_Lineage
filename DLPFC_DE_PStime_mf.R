source('MiscPreprocessing.R')


library(synapser)


tcxCPMObj <- synapser::synGet('syn8456638')
Dat <- read.delim(tcxCPMObj$path,stringsAsFactors = F)

#synapse id of dat2 file: syn8466814
tcxCovObj <- synapser::synGet('syn11024258')
Dat2 <- read.delim(tcxCovObj$path,stringsAsFactors = F)



foobar <- synapser::synGet('syn11695124')
foobar2 <- data.table::fread(foobar$path,data.table=F)

de_file <- synapser::synGet('syn8456721')
de1 <- data.table::fread(de_file$path,data.table=F)
de3 <- dplyr::filter(de1,Model=='Diagnosis',Comparison=='AD-CONTROL')
AMP_mods <- data.frame(GeneID=de3$ensembl_gene_id,logPV= - log(de3$adj.P.Val)/log(10), stringsAsFactors=F)
AMP_mods <- dplyr::filter(AMP_mods,GeneID%in% foobar2$GeneID)
#subsetting genes based on differential expression 
#AMP_mods <-  read.csv('Data/DLPFC_DE.csv')
In <- which(AMP_mods$logPV >= 1)
AMP_mods <- AMP_mods[In,]
#In <- which(AMP_mods$LogPV.DLPFC >= 1)
#AMP_mods <- AMP_mods[In,]




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


DatNorm4 <- DatNorm3
Dat4 <- Dat3

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
temp2$APOE4 <- as.character(temp2$apoe_genotype)
temp2$cogdx <- as.character(temp2$cogdx)

gene_short_name <- Make.Gene.Symb(GeneNamesAD)

rownames(temp) <- NULL
colnames(temp) <- NULL
MonRun <- RunMonocleTobit(temp, temp2, C_by = 'Pseudotime',gene_short_name = gene_short_name)

plot_cell_trajectory(MonRun, color_by = "Diagnosis")


p <- ggplot(MonRun@phenoData@data, aes(x=Diagnosis, y=Pseudotime,fill=Diagnosis)) + geom_boxplot()
p + stat_summary(fun.y=mean, geom="point", shape=23, size=2)  + geom_jitter(width = 0.2)

p <- ggplot(MonRun@phenoData@data, aes(x=APOE4, y=Pseudotime,fill=APOE4)) + geom_boxplot()
p + stat_summary(fun.y=mean, geom="point", shape=23, size=2)  + geom_jitter(width = 0.2)

p <- ggplot(MonRun@phenoData@data, aes(x=cogdx, y=Pseudotime,fill=cogdx)) + geom_boxplot()
p + stat_summary(fun.y=mean, geom="point", shape=23, size=2)  + geom_jitter(width = 0.2)

#p <- ggplot(MonRun@phenoData@data, aes(x=braaksc, y=Pseudotime,fill=braaksc)) + geom_boxplot()
#p + stat_summary(fun.y=mean, geom="point", shape=23, size=2)  + geom_jitter(width = 0.2)

#p <- ggplot(MonRun@phenoData@data, aes(x=ceradsc, y=Pseudotime,fill=ceradsc)) + geom_boxplot()
#p + stat_summary(fun.y=mean, geom="point", shape=23, size=2)  + geom_jitter(width = 0.2)

## plot t-SNE 
library(Rtsne)
Temp <- Rtsne(t(DatNorm4))
df <- as.data.frame(list(Y_1=Temp$Y[,1], Y_2 = Temp$Y[,2],
                         Diagnosis = Dat4$Diagnosis))
ggplot(df,aes(x=Y_1,y=Y_2,color=Diagnosis))+geom_point()

## plot PCA 
Temp <- prcomp(DatNorm4, center = TRUE,scale. = TRUE)
df <- as.data.frame(list(Y_1=Temp$rotation[,1], 
                         Y_2 = Temp$rotation[,2],
                         Diagnosis = Dat4$Diagnosis))
ggplot(df,aes(x=Y_1,y=Y_2,color=Diagnosis))+geom_point()