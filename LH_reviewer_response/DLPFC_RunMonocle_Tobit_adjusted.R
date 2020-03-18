setwd("~/AMPAD_Lineage/LH_reviewer_response/DLPFC_LH_reviewer_response")
source('MiscPreprocessing.R')

#Loading data
synapser::synLogin()

#load rosmap filtered counts logCPM:
tcxCPMObj <- synapser::synGet('syn8456638')
Dat <- read.delim(tcxCPMObj$path,stringsAsFactors = F)

#synapse id of dat2 file (rosmap covariates): syn8466814
tcxCovObj <- synapser::synGet('syn11024258')
Dat2 <- read.delim(tcxCovObj$path,stringsAsFactors = F)

foobar <- synapser::synGet('syn11695124')
foobar2 <- data.table::fread(foobar$path,data.table=F)

#load ROSMAP_DLPFC_DiffExpression.tsv
de_file <- synapser::synGet('syn8456721')
de1 <- data.table::fread(de_file$path,data.table=F)
de3 <- dplyr::filter(de1,Model=='Diagnosis',Comparison=='AD-CONTROL')
AMP_mods <- data.frame(GeneID=de3$ensembl_gene_id,logPV= - log(de3$adj.P.Val)/log(10), stringsAsFactors=F)
AMP_mods <- dplyr::filter(AMP_mods,GeneID%in% foobar2$GeneID)
#subsetting genes based on differential expression 
#AMP_mods <-  read.csv('Data/DLPFC_DE.csv')
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
DatNorm2 <- DatNorm[In_genes,]
GeneNamesAD <- GeneNames[In_genes]

#removing bad batches
DatNorm2 <- DatNorm2[,Dat2$Batch<7]
Dat2 <- Dat2[Dat2$Batch<7,] 



DatNorm3 <- DatNorm2
Dat3 <- Dat2

#Keeping only female data 
#Sex <- 'FEMALE'
In_S <- which(Dat3$msex == 0)
DatNorm4 <- DatNorm3[,In_S]
Dat4 <- Dat3[In_S,]

# In_cov <- which(Cov$ID %in% Dat4$Donor_ID)
# Cov <- Cov[In_cov,]
# In_cov <- c()
# for(i in 1:length(Dat4$Donor_ID)){
#   temp <- which(Cov$ID == Dat4$Donor_ID[i])
#   In_cov <- c(In_cov,temp[1])
# }
# Cov <- Cov[In_cov,]
# 
# for (i in 23:26){
#   Cov[,i] <- (Cov[,i] - min(Cov[,i]))/(max(Cov[,i])-min(Cov[,i]))
# }

source('LineageFunctions.R')
temp <- DatNorm4
temp2 <- Dat4
#temp2$APOE4 <- as.character(temp2$APOE4)
#temp2$braaksc <- as.character(temp2$braaksc)
#temp2$ceradsc <- as.character(temp2$ceradsc)
#temp2$cogdx.1 <- as.character(temp2$cogdx.1)

#converting ENSG to gene symbols
#convert to gene symbol
convertEnsemblToHgnc <- function(ensemblIds){
  
  ensembl=biomaRt::useMart('ENSEMBL_MART_ENSEMBL',
                           dataset = 'hsapiens_gene_ensembl',
                           host='www.ensembl.org')
  
  genes<-getBM(attributes = c('ensembl_gene_id','external_gene_name'),
               filters='ensembl_gene_id',
               values=ensemblIds,
               mart=ensembl)
  return(genes)
}
Make.Gene.Symb <- function(GeneENSG){
  
  #source('convertEnsemblToHgnc.R')
  GeneConv <- convertEnsemblToHgnc(GeneENSG)
  Symb <- as.character(c(1:length(GeneENSG)))
  
  for (i in 1:length(GeneENSG)){
    In <- which(GeneConv$ensembl_gene_id == GeneENSG[i])
    if (length(In)>0){
      Symb[i] <- GeneConv$external_gene_name[In]
    }
  }
  
  return(Symb)
  
}
gene_short_name <- Make.Gene.Symb(GeneNamesAD)


# tcxLineageTimes <- synapser::synTableQuery("select * from syn17023721")$asDataFrame()
dlpfcLineageTimes <- synapser::synTableQuery("select * from syn17023795")$asDataFrame()
# tcxCovObj <- synapser::synGet("syn8466814")
# tcxCov <- data.table::fread(tcxCovObj$path,data.table=F)
# dlpfcCovObj <- synapser::synGet("syn11024258")
# dlpfcCov <- data.table::fread(dlpfcCovObj$path,data.table=F)
# tcxLineageTimes <- tcxLineageTimes[,-c(1:3)]
# dlpfcLineageTimes <- dlpfcLineageTimes[,-c(1:3)]
# dlpfc <- dplyr::left_join(dlpfcLineageTimes,dlpfcCov,by='SampleID')
# tcx <- dplyr::left_join(tcxLineageTimes,tcxCov,by='SampleID')
# 
# tcxAD <- rep(NA,nrow(tcx))
# tcxAD[tcx$Tissue.Diagnosis == 'TCX.AD'] <- 1
# tcxAD[tcx$Tissue.Diagnosis == 'TCX.CONTROL'] <- 0
# tcxDf <- data.frame(diagnosis=tcxAD,
#                     pseudotime=tcx$Pseudotime,
#                     stringsAsFactors = FALSE)
# tcxDf$pseudotime <- scale(tcxDf$pseudotime,center=F)
# summary(glm(diagnosis ~ pseudotime,tcxDf,family='binomial'))
# 
# dlpfcAD <- rep(NA,nrow(dlpfc))
# dlpfcAD[dlpfc$Diagnosis == 'AD'] <- 1
# dlpfcAD[dlpfc$Diagnosis == 'CONTROL'] <- 0
# dlpfcDf <- data.frame(diagnosis = dlpfcAD,
#                       pseudotime=dlpfc$Pseudotime,
#                       stringsAsFactors = FALSE)
# dlpfcDf$pseudotime <- scale(dlpfcDf$pseudotime,center=F)
# 
# tcxDf$BrainRegion <- 'TCX'
# dlpfcDf$BrainRegion <- 'DLPFC'
# combinedDf <- rbind(tcxDf,dlpfcDf)
# combinedDf <- combinedDf[!is.na(combinedDf$diagnosis),]
# combinedDf$diagnosis[combinedDf$diagnosis==1] <- 'AD'
# combinedDf$diagnosis[combinedDf$diagnosis==0] <- 'Control'
# combinedDf$diagnosis <- factor(combinedDf$diagnosis,levels=c('AD','Control'))
# combinedDf$BrainRegion <- as.factor(combinedDf$BrainRegion)
# 
# 
# rosmapObj <- synapser::synGet('syn3191087')
# rosmap <- data.table::fread(rosmapObj$path,data.table=F)
# 
# rosmapIdObj <- synapser::synGet('syn3382527')
# rosmapId <- data.table::fread(rosmapIdObj$path,data.table=F)
# rosmapId <- dplyr::select(rosmapId,projid,rnaseq_id)
# 
# rosmapRNAid<-dplyr::left_join(rosmapId,rosmap)
# dlpfcComplete<-dplyr::left_join(dlpfcLineageTimes,rosmapRNAid,by=c('SampleID'='rnaseq_id'))
# dlpfcComplete$Pseudotime <- scale(dlpfcComplete$Pseudotime,center=F)
# 
# dlpfcComplete <- dplyr::select(dlpfcComplete,Pseudotime,braaksc,ceradsc,cogdx,SampleID)
# 
# 
# dlpfcComplete$braaksc <- factor(dlpfcComplete$braaksc,levels = c(0:6))
# dlpfcComplete$ceradsc <- factor(dlpfcComplete$ceradsc,levels = c(1:4))
# dlpfcComplete$cogdxNew <- rep(NA,nrow(dlpfcComplete))
# dlpfcComplete$cogdxNew[dlpfcComplete$cogdx==1] <- 'NCI'
# dlpfcComplete$cogdxNew[dlpfcComplete$cogdx==2] <- 'MCI'
# dlpfcComplete$cogdxNew[dlpfcComplete$cogdx==4] <- 'LOAD'
# dlpfcComplete$cogdxNew <- factor(dlpfcComplete$cogdxNew,levels = c('NCI','MCI','LOAD'))
# 
# dlpfcComplete <- dlpfcComplete[!duplicated(dlpfcComplete),]
# temp2 <- dplyr::left_join(temp2,dlpfcComplete,by='SampleID')
# 
# temp2$cogdx.x <- as.factor(temp2$cogdx.x)


#running monocle
rownames(temp) <- NULL
colnames(temp) <- NULL

detach("package:synapser", unload = TRUE)
unloadNamespace("PythonEmbedInR") 
library(monocle)
MonRun <- RunMonocleTobit(temp, temp2, C_by = 'Pseudotime',gene_short_name = gene_short_name)
plot_cell_trajectory(MonRun,color_by = "Diagnosis",show_branch_points=F,use_color_gradient = F,cell_size = 0.5)

#compare my pseudotime estimates to what is store in the project
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

write.csv(l,'DLPFC_F_pv1_StatePhenotype_replicated.csv')

#tiff(filename='figure2A.tiff',height=85,width=85,units='mm',res=300)
#tiff(file='~/Desktop/MANUSCRIPT/figure2a.tiff',height=85,width=100,units='mm',res=300)
g<- plot_cell_trajectory(MonRun,color_by = "Diagnosis",show_branch_points=F,use_color_gradient = F,cell_size = 0.5)
g <- g + ggplot2::scale_color_viridis_d()
g
#dev.off()

source('MonoclePlottingCustom.R')

#jpeg('./PaperFigs/DLPFC_diag_new.jpg',width = 600, height = 350)
#tiff(file='~/Desktop/MANUSCRIPT/figure3a1.tiff',height=85,width=100,units='mm',res=300)
g<- plot_cell_trajectory(MonRun,color_by = "braaksc",show_branch_points=F,use_color_gradient = F,cell_size = 0.5)
g <- g + ggplot2::scale_color_viridis_d()
g <- g + ggplot2::labs(color="Braak Score")
g
#dev.off()

MonRun$apoe_genotype <- factor(MonRun$apoe_genotype,levels=c(0,1,2))
tiff(file='~/Desktop/MANUSCRIPT/figureS6a.tiff',height=85,width=100,units='mm',res=300)
g<- plot_cell_trajectory(MonRun,color_by = "apoe_genotype",show_branch_points=F,use_color_gradient = F,cell_size = 0.5)
g <- g + ggplot2::scale_color_viridis_d()
g <- g + ggplot2::labs(color="APOE e4 Dosage")
g
dev.off()


tiff(file='~/Desktop/MANUSCRIPT/figure3a2.tiff',height=85,width=100,units='mm',res=300)
g<- plot_cell_trajectory(MonRun,color_by = "ceradsc",show_branch_points=F,use_color_gradient = F,cell_size = 0.5)
g <- g + ggplot2::scale_color_viridis_d()
g <- g + ggplot2::labs(color="CERAD Score")
g
dev.off()

tiff(file='~/Desktop/MANUSCRIPT/figure3a3.tiff',height=85,width=100,units='mm',res=300)
g<- plot_cell_trajectory(MonRun,color_by = "cogdx.x",show_branch_points=F,use_color_gradient = F,cell_size = 0.5)
g <- g + ggplot2::scale_color_viridis_d()
g <- g + ggplot2::labs(color="Cognitive Diagnosis")
g
dev.off()
#dev.off()

tiff(file='~/Desktop/MANUSCRIPT/figure3b1.tiff',height=85,width=100,units='mm',res=300)
g <- ggplot2::ggplot(MonRun@phenoData@data, aes(x=braaksc, y=scale(Pseudotime,center=F),fill=braaksc)) 
g <- g + ggplot2::geom_boxplot()
g <- g + ggplot2::stat_summary(fun.y=mean, geom="point", shape=23, size=2)
g <- g + ggplot2::theme(axis.text=element_text(size=15), axis.title=element_text(size=15,face="bold"),
                        legend.text=element_text(size=15)) 
g <- g + ggplot2::scale_fill_viridis_d()
g <- g + ggplot2::labs(fill="Braak\nScore",y="Pseudotime",x="Braak Score")
g
dev.off()


tiff(file='~/Desktop/MANUSCRIPT/figure3b2.tiff',height=85,width=100,units='mm',res=300)
g <- ggplot2::ggplot(MonRun@phenoData@data, aes(x=ceradsc, y=scale(Pseudotime,center=F),fill=ceradsc)) 
g <- g + ggplot2::geom_boxplot()
g <- g + ggplot2::stat_summary(fun.y=mean, geom="point", shape=23, size=2)
g <- g + ggplot2::theme(axis.text=element_text(size=15), axis.title=element_text(size=15,face="bold"),
                        legend.text=element_text(size=15)) 
g <- g + ggplot2::scale_fill_viridis_d()
g <- g + ggplot2::labs(fill="CERAD\nScore",y="Pseudotime",x="CERAD Score")
g
dev.off()


tiff(file='~/Desktop/MANUSCRIPT/figure3b3.tiff',height=85,width=100,units='mm',res=300)
g <- ggplot2::ggplot(MonRun@phenoData@data, aes(x=cogdx.x, y=scale(Pseudotime,center=F),fill=cogdx.x)) 
g <- g + ggplot2::geom_boxplot()
g <- g + ggplot2::stat_summary(fun.y=mean, geom="point", shape=23, size=2)
g <- g + ggplot2::theme(axis.text=element_text(size=15), axis.title=element_text(size=15,face="bold"),
                        legend.text=element_text(size=15)) 
g <- g + ggplot2::scale_fill_viridis_d()
g <- g + ggplot2::labs(fill="Cognitive\nDiagnosis",y="Pseudotime",x="Cognitive Diagnosis")
g
dev.off()



MonRun$State2 <- MonRun$State
MonRun$State2[MonRun$State == 4] <- 3
MonRun$State2[MonRun$State == 7] <- 3
MonRun$State2[MonRun$State == 6] <- 4
MonRun$State2[MonRun$State == 5] <- 4
MonRun$State2[MonRun$State == 8] <- 5
MonRun$State2[MonRun$State == 9] <- 6

MonRun$State2 <- as.numeric(MonRun$State2)
MonRun$State2 <- as.factor(MonRun$State2)

tiff(file='~/Desktop/MANUSCRIPT/figureS6.tiff',height=85,width=100,units='mm',res=300)
g<- plot_cell_trajectory(MonRun,color_by = "State2",show_branch_points=F,use_color_gradient = F,cell_size = 0.5)
g <- g + ggplot2::scale_color_viridis_d()
g <- g + ggplot2::labs(color="Lineage\nState")
g
dev.off()


apoeDf <- data.frame(apoe=MonRun$apoe_genotype==1 | MonRun$apoe_genotype==2,
                     State=MonRun$State2,
                     stringsAsFactors=F)

summary(glm(apoe ~ State,apoeDf,family='binomial'))



######Change the function from LineageFunctions.R to adjust for RIN####
RunMonocleTobit_adj <- function(Dat, Labels, max_components=2, meth = 'DDRTree',  C_by = NULL, 
                                gene_short_name = NULL){ 
  
  library(monocle)
  
  HSMM_expr_matrix <- Dat
  names(HSMM_expr_matrix)<-seq(1,dim(Dat)[2])
  
  if(is.null(gene_short_name)){
    gene_short_name <- c(1:dim(Dat)[1])
  }
  
  
  gene_short_name <- data.frame(gene_short_name)
  Labels <- data.frame(Labels)
  rownames(Labels) <- seq(1,dim(Dat)[2])
  
  pd <- new("AnnotatedDataFrame", data = Labels)
  fd <- new("AnnotatedDataFrame", data = gene_short_name)
  
  
  
  HSMM <- newCellDataSet(as.matrix(HSMM_expr_matrix),
                         phenoData = pd,
                         featureData = fd,
                         expressionFamily=tobit())
  
  HSMM <- reduceDimension(HSMM, max_components=max_components, reduction_method = meth, residualModelFormulaStr = "~RINcontinuous")
  HSMM <- orderCells(HSMM)
  if(is.null(C_by)){
    plot_cell_trajectory(HSMM, color_by="Labels")
  }
  else{
    plot_cell_trajectory(HSMM, color_by=C_by)
  }
  
  
  return(HSMM)
  
}
MonRun_adj <- RunMonocleTobit_adj(temp, temp2, C_by = 'Pseudotime',gene_short_name = gene_short_name)
plot_cell_trajectory(MonRun_adj, color_by = "Diagnosis")
head(MonRun$Pseudotime)
head(MonRun_adj$Pseudotime)

#Create differential expression based on state 
MonRun_adj2 = MonRun_adj
ScaledDat = ScaleMatrix(MonRun_adj2@assayData$exprs)
#MonRun2@assayData$exprs <- ScaledDat
MonRun_adj$State2 <- MonRun_adj$State
MonRun_adj$State2[MonRun_adj$State == 4] <- 3
MonRun_adj$State2[MonRun_adj$State == 7] <- 3
MonRun_adj$State2[MonRun_adj$State == 6] <- 5


#find centers for each state 
UnqStates <-as.vector(unique(MonRun_adj$State))
C_s1 <- rep(0,length(UnqStates))
C_s2 <- rep(0,length(UnqStates))

for (i in 1:length(UnqStates)){
  In_s <- which(MonRun_adj$State==UnqStates[i])
  C_s1[i] <- mean(MonRun_adj@reducedDimS[1,In_s])
  C_s2[i] <- mean(MonRun_adj@reducedDimS[2,In_s])
  
}

D_mean <- rep(0,length(MonRun_adj$State))

#For each sample, find the distance to the state mean 
for(i in 1:length(MonRun_adj$State)){
  In_S <- which(UnqStates==MonRun_adj$State[i])
  D_mean[i] <- sqrt((MonRun_adj@reducedDimS[1,i]-C_s1[In_S])**2 + (MonRun_adj@reducedDimS[2,i]-C_s2[In_S])**2)
  
}

#Normalize each sample's distance based on max distance for that state
D_meanNorm <- rep(0,length(MonRun_adj$State))

for(i in 1:length(MonRun_adj$State)){
  In_S <- which(MonRun_adj$State==MonRun_adj$State[i])
  D_meanNorm[i] <- (D_mean[i]-min(D_mean[In_S]))/max(D_mean[In_S])
  
}

#Make new dataframe with SampleID, D_mean, D_meanNorm, State, Pseudotime
l_adj <- list()

MonRun_adj$D_mean <- D_mean
MonRun_adj$D_meanNorm <- D_meanNorm

l_adj$SampleID <- MonRun_adj$SampleID
l_adj$D_mean <- MonRun_adj$D_mean
l_adj$D_meanNorm <- MonRun_adj$D_meanNorm
l_adj$State <- MonRun_adj$State
l_adj$Pseudotime <- MonRun_adj$Pseudotime

l_adj <- as.data.frame(l_adj)
write.csv(l_adj,'DLPFC_F_pv1_StatePhenotype_RINadj.csv')
