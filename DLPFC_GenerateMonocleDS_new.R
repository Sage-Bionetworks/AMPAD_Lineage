source('MiscPreprocessing.R')

#Loading data
synapser::synLogin()

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
gene_short_name <- Make.Gene.Symb(GeneNamesAD)


tcxLineageTimes <- synapser::synTableQuery("select * from syn17023721")$asDataFrame()
dlpfcLineageTimes <- synapser::synTableQuery("select * from syn17023795")$asDataFrame()
tcxCovObj <- synapser::synGet("syn8466814")
tcxCov <- data.table::fread(tcxCovObj$path,data.table=F)
dlpfcCovObj <- synapser::synGet("syn11024258")
dlpfcCov <- data.table::fread(dlpfcCovObj$path,data.table=F)
tcxLineageTimes <- tcxLineageTimes[,-c(1:3)]
dlpfcLineageTimes <- dlpfcLineageTimes[,-c(1:3)]
dlpfc <- dplyr::left_join(dlpfcLineageTimes,dlpfcCov,by='SampleID')
tcx <- dplyr::left_join(tcxLineageTimes,tcxCov,by='SampleID')

tcxAD <- rep(NA,nrow(tcx))
tcxAD[tcx$Tissue.Diagnosis == 'TCX.AD'] <- 1
tcxAD[tcx$Tissue.Diagnosis == 'TCX.CONTROL'] <- 0
tcxDf <- data.frame(diagnosis=tcxAD,
                    pseudotime=tcx$Pseudotime,
                    stringsAsFactors = FALSE)
tcxDf$pseudotime <- scale(tcxDf$pseudotime,center=F)
summary(glm(diagnosis ~ pseudotime,tcxDf,family='binomial'))

dlpfcAD <- rep(NA,nrow(dlpfc))
dlpfcAD[dlpfc$Diagnosis == 'AD'] <- 1
dlpfcAD[dlpfc$Diagnosis == 'CONTROL'] <- 0
dlpfcDf <- data.frame(diagnosis = dlpfcAD,
                      pseudotime=dlpfc$Pseudotime,
                      stringsAsFactors = FALSE)
dlpfcDf$pseudotime <- scale(dlpfcDf$pseudotime,center=F)

tcxDf$BrainRegion <- 'TCX'
dlpfcDf$BrainRegion <- 'DLPFC'
combinedDf <- rbind(tcxDf,dlpfcDf)
combinedDf <- combinedDf[!is.na(combinedDf$diagnosis),]
combinedDf$diagnosis[combinedDf$diagnosis==1] <- 'AD'
combinedDf$diagnosis[combinedDf$diagnosis==0] <- 'Control'
combinedDf$diagnosis <- factor(combinedDf$diagnosis,levels=c('AD','Control'))
combinedDf$BrainRegion <- as.factor(combinedDf$BrainRegion)


rosmapObj <- synapser::synGet('syn3191087')
rosmap <- data.table::fread(rosmapObj$path,data.table=F)

rosmapIdObj <- synapser::synGet('syn3382527')
rosmapId <- data.table::fread(rosmapIdObj$path,data.table=F)
rosmapId <- dplyr::select(rosmapId,projid,rnaseq_id)

rosmapRNAid<-dplyr::left_join(rosmapId,rosmap)
dlpfcComplete<-dplyr::left_join(dlpfcLineageTimes,rosmapRNAid,by=c('SampleID'='rnaseq_id'))
dlpfcComplete$Pseudotime <- scale(dlpfcComplete$Pseudotime,center=F)

dlpfcComplete <- dplyr::select(dlpfcComplete,Pseudotime,braaksc,ceradsc,cogdx,SampleID)


dlpfcComplete$braaksc <- factor(dlpfcComplete$braaksc,levels = c(0:6))
dlpfcComplete$ceradsc <- factor(dlpfcComplete$ceradsc,levels = c(1:4))
dlpfcComplete$cogdxNew <- rep(NA,nrow(dlpfcComplete))
dlpfcComplete$cogdxNew[dlpfcComplete$cogdx==1] <- 'NCI'
dlpfcComplete$cogdxNew[dlpfcComplete$cogdx==2] <- 'MCI'
dlpfcComplete$cogdxNew[dlpfcComplete$cogdx==4] <- 'LOAD'
dlpfcComplete$cogdxNew <- factor(dlpfcComplete$cogdxNew,levels = c('NCI','MCI','LOAD'))

dlpfcComplete <- dlpfcComplete[!duplicated(dlpfcComplete),]
temp2 <- dplyr::left_join(temp2,dlpfcComplete,by='SampleID')

temp2$cogdx.x <- as.factor(temp2$cogdx.x)
#running monocle
rownames(temp) <- NULL
colnames(temp) <- NULL
MonRun <- RunMonocleTobit(temp, temp2, C_by = 'Pseudotime',gene_short_name = gene_short_name)

#tiff(filename='figure2A.tiff',height=85,width=85,units='mm',res=300)
tiff(file='~/Desktop/MANUSCRIPT/figure2a.tiff',height=85,width=100,units='mm',res=300)
g<- plot_cell_trajectory(MonRun,color_by = "Diagnosis",show_branch_points=F,use_color_gradient = F,cell_size = 0.5)
g <- g + ggplot2::scale_color_viridis_d()
g
dev.off()

source('MonoclePlottingCustom.R')

#jpeg('./PaperFigs/DLPFC_diag_new.jpg',width = 600, height = 350)
tiff(file='~/Desktop/MANUSCRIPT/figure3a1.tiff',height=85,width=100,units='mm',res=300)
g<- plot_cell_trajectory(MonRun,color_by = "braaksc",show_branch_points=F,use_color_gradient = F,cell_size = 0.5)
g <- g + ggplot2::scale_color_viridis_d()
g <- g + ggplot2::labs(color="Braak Score")
g
dev.off()

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


