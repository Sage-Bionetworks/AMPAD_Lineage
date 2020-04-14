
#Upload pseudotimes calculated in males only
tcxLineageTimesLH <- read.csv("~/AMPAD_Lineage/LH_reviewer_response/Monly_MF_analyses/TCX_M_PStimes.csv", stringsAsFactors = FALSE)
dlpfcLineageTimesLH <- read.csv("~/AMPAD_Lineage/LH_reviewer_response/Monly_MF_analyses/DLPFC_M_PStimes.csv", stringsAsFactors = FALSE)
tcxCovObj <- synapser::synGet("syn8466814")
tcxCov <- data.table::fread(tcxCovObj$path,data.table=F)
dlpfcCovObj <- synapser::synGet("syn11024258")
dlpfcCov <- data.table::fread(dlpfcCovObj$path,data.table=F)
tcxLineageTimesLH <- tcxLineageTimesLH[,-c(1)]
dlpfcLineageTimesLH <- dlpfcLineageTimesLH[,-c(1)]
dlpfc <- dplyr::left_join(dlpfcLineageTimesLH,dlpfcCov,by='SampleID')
tcx <- dplyr::left_join(tcxLineageTimesLH,tcxCov,by='SampleID')


#Next run a logistic regression
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
summary(glm(diagnosis ~ pseudotime,dlpfcDf,family='binomial'))


#Next make a ggplot boxplot and save.  First combine into a single df
tcxDf$BrainRegion <- 'TCX'
dlpfcDf$BrainRegion <- 'DLPFC'
combinedDf <- rbind(tcxDf,dlpfcDf)
combinedDf <- combinedDf[!is.na(combinedDf$diagnosis),]
combinedDf$diagnosis[combinedDf$diagnosis==1] <- 'AD'
combinedDf$diagnosis[combinedDf$diagnosis==0] <- 'Control'
combinedDf$diagnosis <- factor(combinedDf$diagnosis,levels=c('AD','Control'))
combinedDf$BrainRegion <- as.factor(combinedDf$BrainRegion)
```
#run ggplot2

tiff(file='~/AMPAD_Lineage/LH_reviewer_response/Monly_MF_analyses/M_ADcasecontrol_PS_boxplots.tiff',height=85,width=100,units='mm',res=300)
g <- ggplot2::ggplot(combinedDf,ggplot2::aes(x=BrainRegion,
                                             y=pseudotime,
                                             color=diagnosis))
g <- g + ggplot2::geom_boxplot() 
g <- g + ggplot2::geom_point(position=ggplot2::position_jitterdodge())
#g <- g + ggplot2::scale_color_viridis_d()
g <- g + ggplot2::scale_color_manual(values=viridis::viridis(3)[1:2])
g
dev.off()




#GWAS LOAD distributions
#(needed to tweak the function to change the host site)

convertEnsemblToHgnc <- function(ensemblIds){
  
  ensembl=biomaRt::useMart('ENSEMBL_MART_ENSEMBL',
                           dataset = 'hsapiens_gene_ensembl',
                           host='useast.ensembl.org')
  
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


#add expression data to tcx data frame, calculate correlations
tcxCPMObj <- synapser::synGet('syn8466816')
Dat <- data.table::fread(tcxCPMObj$path,data.table=F)
sampleIds <- colnames(Dat)[-1]
geneIds <- Dat$ensembl_gene_id
Dat <- Dat[,-1]
Dat <- t(Dat)
colnames(Dat) <- geneIds
Dat <- data.frame(Dat,stringsAsFactors=F)
Dat$sampleId <- sampleIds
tcx <- dplyr::left_join(tcx,Dat,by=c('SampleID'='sampleId'))
tcx2 <- dplyr::select(tcx,dplyr::starts_with("ENSG"))
corvec <- cor(tcx2,tcx$Pseudotime,method='spearman')
corDf <- data.frame(geneid=colnames(tcx2),cor=corvec,stringsAsFactors=F)
corDf2 <- corDf
corDf2$external_gene_name <- Make.Gene.Symb(corDf2$geneid)
corDf2$cor <- NULL
corDf2 <- dplyr::left_join(corDf,corDf2,by=c('geneid'))


#and add expression data to dlpfc data frame, calculate correlations
dlpfcCPMObj <- synapser::synGet('syn8456638')
Dat <- data.table::fread(dlpfcCPMObj$path,data.table=F)
sampleIds <- colnames(Dat)[-1]
geneIds <- Dat$ensembl_gene_id
Dat <- Dat[,-1]
Dat <- t(Dat)
colnames(Dat) <- geneIds
Dat <- data.frame(Dat,stringsAsFactors=F)
Dat$sampleId <- sampleIds
dlpfc <- dplyr::left_join(dlpfc,Dat,by=c('SampleID'='sampleId'))
dlpfc2 <- dplyr::select(dlpfc,dplyr::starts_with("ENSG"))
corvec <- cor(dlpfc2,dlpfc$Pseudotime,method='spearman')
corDfdlpfc <- data.frame(geneid=colnames(dlpfc2),cor=corvec,stringsAsFactors=F)
corDfdlpfc2 <- corDfdlpfc
corDfdlpfc2$external_gene_name <- Make.Gene.Symb(corDfdlpfc2$geneid)
corDfdlpfc2$cor <- NULL
corDfdlpfc2 <- dplyr::left_join(corDfdlpfc,corDfdlpfc2,by=c('geneid'))


#input desired gwas genes
ad_gwas <- c("CR1",
             "BIN1",
             "INPP5D",
             "HLA-DRB1",
             "TREM2",
             "MEF2C",
             "NME8",
             "CD2AP",
             "NYAP1",
             "EPHA1",
             "PTK2B",
             "CLU",
             "SPI1",
             "MS4A2",
             "PICALM",
             "SORL1",
             "FERMT2",
             "SLC24A4",
             "ABCA7",
             "APOE",
             "CASS4",
             "ECHDC3",
             "ACE",
             "NDUFAF6",
             "ECHDC3",
             "ADAMTS20",
             "SPPL2A",
             "ADAM10",
             "IQCK",
             "MIR142",
             "ACE",
             "ADAMTS1",
             "SUCLG2P4",
             "FST",
             "OARD1",
             "WWOX",
             "MAF",
             "CD55",
             "YOD1",
             "HLA-DRB1",
             "PSMB8",
             "C4A",
             "GPSM3",
             "HLA-DPA1",
             "HLA-DQA1",
             "HLA-DRA",
             "HLA-DRB5",
             "PSMB9",
             "CD2AP",
             "AGFG2",
             "PILRA",
             "EPHB4",
             "C7orf43",
             "GAL3ST4",
             "ZKSCCAN1",
             "FAM131B",
             "PSMC3",
             "ACP2",
             "C1QTNF4",
             "CELF1",
             "MTCH2",
             "NDUFS3",
             "NUP160",
             "MS4A6A",
             "MS4A7",
             "MS4A4A",
             "EED",
             "PICALM",
             "STYX",
             "RIN3",
             "HMHA1",
             "CNN2",
             "WDR18",
             "CASS4")

corDf <- corDf2
corDfdlpfc <- corDfdlpfc2
mean(abs(corDf$cor))
mean(abs(corDf[corDf$external_gene_name %in% ad_gwas,]$cor))
mean(abs(corDfdlpfc$cor))
mean(abs(corDf[corDfdlpfc$external_gene_name %in% ad_gwas,]$cor))


corDf$adGwas <- corDf$external_gene_name %in% ad_gwas
corDfdlpfc$adGwas <- corDfdlpfc$external_gene_name %in% ad_gwas

corDf$brainRegion <- 'TCX'
corDfdlpfc$brainRegion <- 'DLPFC'
corDfcombined <- rbind(corDf,corDfdlpfc)
colnames(corDfcombined)[4] <- 'LOADGWASGene'

tiff(file='~/AMPAD_Lineage/LH_reviewer_response/Monly_MF_analyses/M_gwasload_boxplot.tiff',height=85,width=100,units='mm',res=300)
g <- ggplot2::ggplot(corDfcombined,ggplot2::aes(x=brainRegion,y=cor,fill=LOADGWASGene))
g <- g + ggplot2::geom_boxplot()
g <- g + ggplot2::scale_fill_viridis_d()
#g <- g + ggplot2::geom_point(position=ggplot2::position_jitterdodge())
g <- g + ggplot2::labs(x = 'Brain Region',y='Correlation with pseudotime',fill='LOAD\nGWAS\nGene')
g
# tiff(file='~/Desktop/MANUSCRIPT/figure2d.tiff',height=85,width=100,units='mm',res=300)

dev.off()


##############REPEAT FOR MALES/FEMALES COMBINED############################

#Upload pseudotimes calculated in males/females combined
tcxLineageTimesLH <- read.csv("~/AMPAD_Lineage/LH_reviewer_response/Monly_MF_analyses/TCX_MF_PStimes.csv", stringsAsFactors = FALSE)
dlpfcLineageTimesLH <- read.csv("~/AMPAD_Lineage/LH_reviewer_response/Monly_MF_analyses/DLPFC_MF_PStimes.csv", stringsAsFactors = FALSE)
tcxCovObj <- synapser::synGet("syn8466814")
tcxCov <- data.table::fread(tcxCovObj$path,data.table=F)
dlpfcCovObj <- synapser::synGet("syn11024258")
dlpfcCov <- data.table::fread(dlpfcCovObj$path,data.table=F)
tcxLineageTimesLH <- tcxLineageTimesLH[,-c(1)]
dlpfcLineageTimesLH <- dlpfcLineageTimesLH[,-c(1)]
dlpfc <- dplyr::left_join(dlpfcLineageTimesLH,dlpfcCov,by='SampleID')
tcx <- dplyr::left_join(tcxLineageTimesLH,tcxCov,by='SampleID')


#Next run a logistic regression
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
summary(glm(diagnosis ~ pseudotime,dlpfcDf,family='binomial'))


#Next make a ggplot boxplot and save.  First combine into a single df
tcxDf$BrainRegion <- 'TCX'
dlpfcDf$BrainRegion <- 'DLPFC'
combinedDf <- rbind(tcxDf,dlpfcDf)
combinedDf <- combinedDf[!is.na(combinedDf$diagnosis),]
combinedDf$diagnosis[combinedDf$diagnosis==1] <- 'AD'
combinedDf$diagnosis[combinedDf$diagnosis==0] <- 'Control'
combinedDf$diagnosis <- factor(combinedDf$diagnosis,levels=c('AD','Control'))
combinedDf$BrainRegion <- as.factor(combinedDf$BrainRegion)

#run ggplot2

tiff(file='~/AMPAD_Lineage/LH_reviewer_response/Monly_MF_analyses/MF_ADcasecontrol_PS_boxplots.tiff',height=85,width=100,units='mm',res=300)
g <- ggplot2::ggplot(combinedDf,ggplot2::aes(x=BrainRegion,
                                             y=pseudotime,
                                             color=diagnosis))
g <- g + ggplot2::geom_boxplot() 
g <- g + ggplot2::geom_point(position=ggplot2::position_jitterdodge())
#g <- g + ggplot2::scale_color_viridis_d()
g <- g + ggplot2::scale_color_manual(values=viridis::viridis(3)[1:2])
g
dev.off()




#GWAS LOAD distributions
#(needed to tweak the function to change the host site)

convertEnsemblToHgnc <- function(ensemblIds){
  
  ensembl=biomaRt::useMart('ENSEMBL_MART_ENSEMBL',
                           dataset = 'hsapiens_gene_ensembl',
                           host='useast.ensembl.org')
  
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


#add expression data to tcx data frame, calculate correlations
tcxCPMObj <- synapser::synGet('syn8466816')
Dat <- data.table::fread(tcxCPMObj$path,data.table=F)
sampleIds <- colnames(Dat)[-1]
geneIds <- Dat$ensembl_gene_id
Dat <- Dat[,-1]
Dat <- t(Dat)
colnames(Dat) <- geneIds
Dat <- data.frame(Dat,stringsAsFactors=F)
Dat$sampleId <- sampleIds
tcx <- dplyr::left_join(tcx,Dat,by=c('SampleID'='sampleId'))
tcx2 <- dplyr::select(tcx,dplyr::starts_with("ENSG"))
corvec <- cor(tcx2,tcx$Pseudotime,method='spearman')
corDf <- data.frame(geneid=colnames(tcx2),cor=corvec,stringsAsFactors=F)
corDf2 <- corDf
corDf2$external_gene_name <- Make.Gene.Symb(corDf2$geneid)
corDf2$cor <- NULL
corDf2 <- dplyr::left_join(corDf,corDf2,by=c('geneid'))


#and add expression data to dlpfc data frame, calculate correlations
dlpfcCPMObj <- synapser::synGet('syn8456638')
Dat <- data.table::fread(dlpfcCPMObj$path,data.table=F)
sampleIds <- colnames(Dat)[-1]
geneIds <- Dat$ensembl_gene_id
Dat <- Dat[,-1]
Dat <- t(Dat)
colnames(Dat) <- geneIds
Dat <- data.frame(Dat,stringsAsFactors=F)
Dat$sampleId <- sampleIds
dlpfc <- dplyr::left_join(dlpfc,Dat,by=c('SampleID'='sampleId'))
dlpfc2 <- dplyr::select(dlpfc,dplyr::starts_with("ENSG"))
corvec <- cor(dlpfc2,dlpfc$Pseudotime,method='spearman')
corDfdlpfc <- data.frame(geneid=colnames(dlpfc2),cor=corvec,stringsAsFactors=F)
corDfdlpfc2 <- corDfdlpfc
corDfdlpfc2$external_gene_name <- Make.Gene.Symb(corDfdlpfc2$geneid)
corDfdlpfc2$cor <- NULL
corDfdlpfc2 <- dplyr::left_join(corDfdlpfc,corDfdlpfc2,by=c('geneid'))


#input desired gwas genes
ad_gwas <- c("CR1",
             "BIN1",
             "INPP5D",
             "HLA-DRB1",
             "TREM2",
             "MEF2C",
             "NME8",
             "CD2AP",
             "NYAP1",
             "EPHA1",
             "PTK2B",
             "CLU",
             "SPI1",
             "MS4A2",
             "PICALM",
             "SORL1",
             "FERMT2",
             "SLC24A4",
             "ABCA7",
             "APOE",
             "CASS4",
             "ECHDC3",
             "ACE",
             "NDUFAF6",
             "ECHDC3",
             "ADAMTS20",
             "SPPL2A",
             "ADAM10",
             "IQCK",
             "MIR142",
             "ACE",
             "ADAMTS1",
             "SUCLG2P4",
             "FST",
             "OARD1",
             "WWOX",
             "MAF",
             "CD55",
             "YOD1",
             "HLA-DRB1",
             "PSMB8",
             "C4A",
             "GPSM3",
             "HLA-DPA1",
             "HLA-DQA1",
             "HLA-DRA",
             "HLA-DRB5",
             "PSMB9",
             "CD2AP",
             "AGFG2",
             "PILRA",
             "EPHB4",
             "C7orf43",
             "GAL3ST4",
             "ZKSCCAN1",
             "FAM131B",
             "PSMC3",
             "ACP2",
             "C1QTNF4",
             "CELF1",
             "MTCH2",
             "NDUFS3",
             "NUP160",
             "MS4A6A",
             "MS4A7",
             "MS4A4A",
             "EED",
             "PICALM",
             "STYX",
             "RIN3",
             "HMHA1",
             "CNN2",
             "WDR18",
             "CASS4")

corDf <- corDf2
corDfdlpfc <- corDfdlpfc2
mean(abs(corDf$cor))
mean(abs(corDf[corDf$external_gene_name %in% ad_gwas,]$cor))
mean(abs(corDfdlpfc$cor))
mean(abs(corDf[corDfdlpfc$external_gene_name %in% ad_gwas,]$cor))


corDf$adGwas <- corDf$external_gene_name %in% ad_gwas
corDfdlpfc$adGwas <- corDfdlpfc$external_gene_name %in% ad_gwas

corDf$brainRegion <- 'TCX'
corDfdlpfc$brainRegion <- 'DLPFC'
corDfcombined <- rbind(corDf,corDfdlpfc)
colnames(corDfcombined)[4] <- 'LOADGWASGene'

tiff(file='~/AMPAD_Lineage/LH_reviewer_response/Monly_MF_analyses/MF_gwasload_boxplot.tiff',height=85,width=100,units='mm',res=300)
g <- ggplot2::ggplot(corDfcombined,ggplot2::aes(x=brainRegion,y=cor,fill=LOADGWASGene))
g <- g + ggplot2::geom_boxplot()
g <- g + ggplot2::scale_fill_viridis_d()
#g <- g + ggplot2::geom_point(position=ggplot2::position_jitterdodge())
g <- g + ggplot2::labs(x = 'Brain Region',y='Correlation with pseudotime',fill='LOAD\nGWAS\nGene')
g
# tiff(file='~/Desktop/MANUSCRIPT/figure2d.tiff',height=85,width=100,units='mm',res=300)

dev.off()
