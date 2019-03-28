source('MiscPreprocessing.R')



#Read data
Dat <- read.delim('Data/ROSMAP_DLPFC_logCPM.tsv',stringsAsFactors = F)
Dat2 <- read.csv('Data/ROSMAP_DLPFC_Covariates3.csv',stringsAsFactors = F)

#Sub-setting genes
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

DatNorm2 <- DatNorm[In_genes,]
GeneNamesAD <- GeneNames[In_genes]

#Removing bad batches
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


#Running monocle
source('LineageFunctions.R')
temp <- DatNorm4
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

##Scale and smooth data 
#scaling data
ScDat <- DatNorm4

ScDat2 <- ScDat + 0.0 

#fitting df=3 smoothing spline through data
for(i in 1:dim(ScDat)[1]){
  
  ft <- smooth.spline(MonRun@phenoData@data$Pseudotime, ScDat[i,], df = 5)
  ScDat2[i,] <- predict(ft,MonRun@phenoData@data$Pseudotime)$y
  
  if(i%%1000==0){
    print(i)
  }
  
}

ScDat2 <- ScaleMatrix(ScDat2)


#calculate average marker expression for Neurons, Astrocytes, Microglia and Oligodendrocytes 
load('Data/mouseMarkerGenes.rda')
BR_genes <- mouseMarkerGenes$Cortex
for(i in 1:length(names(BR_genes))){
  BR_genes[[i]] <- toupper(BR_genes[[i]])
}
In_a <- which(gene_short_name%in% BR_genes$Astrocyte)
In_n <- which(gene_short_name%in% BR_genes$Pyramidal)
In_m <- which(gene_short_name%in% BR_genes$Microglia)
In_o <- which(gene_short_name%in% BR_genes$Oligo)

M_a <- colMeans(ScDat2[In_a,])
M_n <- colMeans(ScDat2[In_n,])
M_m <- colMeans(ScDat2[In_m,])
M_o <- colMeans(ScDat2[In_o,])

Srtd <- sort(MonRun@phenoData@data$Pseudotime,index.return = T)
Srtd <- Srtd$ix

#plotting 
jpeg('./PaperFigs/DLPFC_CellTypes_new.jpg',width = 600, height = 350)
plot(MonRun@phenoData@data$Pseudotime[Srtd], M_a[Srtd], type = 'l', lwd = 3, col = 'Red', 
     xlab = 'Pseudotime', ylab = 'Mean Norm. Expression', ylim = c(0,1),cex.lab=1.3)
lines(MonRun@phenoData@data$Pseudotime[Srtd], M_n[Srtd], lwd = 3, col = 'Blue')
lines(MonRun@phenoData@data$Pseudotime[Srtd], M_m[Srtd], lwd = 3, col = 'Green')
lines(MonRun@phenoData@data$Pseudotime[Srtd], M_o[Srtd], lwd = 3, col = 'Orange')
legend(0, 0.85, legend=c("Astrocytes", "Neurons","Microglia","Oligodendrocytes"),
       col=c("Red", "Blue","Green","Orange"), lty=1, cex=1.2, pt.cex = 1)
dev.off()