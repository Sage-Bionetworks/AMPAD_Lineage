tcxLineageTimesLH <- read.csv("~/AMPAD_Lineage/LH_reviewer_response/TCX_F_PStimes_adjusted.csv", stringsAsFactors = FALSE)
dlpfcLineageTimesLH <- read.csv("~/AMPAD_Lineage/LH_reviewer_response/DLPFC_F_PStimes_adjusted.csv", stringsAsFactors = FALSE)
tcxCovObj <- synapser::synGet("syn8466814")
tcxCov <- data.table::fread(tcxCovObj$path,data.table=F)
dlpfcCovObj <- synapser::synGet("syn11024258")
dlpfcCov <- data.table::fread(dlpfcCovObj$path,data.table=F)
tcxLineageTimesLH <- tcxLineageTimesLH[,-c(1)]
dlpfcLineageTimesLH <- dlpfcLineageTimesLH[,-c(1)]
dlpfc <- dplyr::left_join(dlpfcLineageTimesLH,dlpfcCov,by='SampleID')
tcx <- dplyr::left_join(tcxLineageTimesLH,tcxCov,by='SampleID')



#needed to tweak the convert function to a different host (biomaRt stopped working)
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

tcxCPMObj <- synapser::synGet('syn8466816')
#Dat <- read.delim(tcxCPMObj$path,stringsAsFactors = F)
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
map <- utilityFunctions::convertEnsemblToHgnc(corDf$geneid)
corDf <- dplyr::left_join(corDf,map,by=c('geneid'='ensembl_gene_id'))


tcxCPMObj <- synapser::synGet('syn8466816')
#Dat <- read.delim(tcxCPMObj$path,stringsAsFactors = F)
Dat_tcx <- data.table::fread(tcxCPMObj$path,data.table=F)
sampleIds <- colnames(Dat_tcx)[-1]
geneIds <- Dat_tcx$ensembl_gene_id
Dat_tcx <- Dat_tcx[,-1]
Dat_tcx <- t(Dat_tcx)
colnames(Dat_tcx) <- geneIds
Dat_tcx <- data.frame(Dat_tcx,stringsAsFactors=F)
Dat_tcx$sampleId <- sampleIds
tcx <- dplyr::left_join(tcx,Dat_tcx,by=c('SampleID'='sampleId'))
tcx2 <- dplyr::select(tcx,dplyr::starts_with("ENSG"))
#corvec <- cor(tcx2,tcx$Pseudotime_unadj,method='spearman')
#corvec <- cor(tcx2,tcx$Pseudotime_PCs,method='spearman')
#corvec <- cor(tcx2,tcx$Pseudotime_RIN,method='spearman')
corvec <- cor(tcx2,tcx$Pseudotime_PMI,method='spearman')
corDf <- data.frame(geneid=colnames(tcx2),cor=corvec,stringsAsFactors=F)
#gene_names <- corDf$geneid
#map <- utilityFunctions::convertEnsemblToHgnc(corDf$geneid)
corDf2 <- corDf
corDf2$external_gene_name <- Make.Gene.Symb(corDf2$geneid)
corDf2$cor <- NULL
corDf2 <- dplyr::left_join(corDf,corDf2,by=c('geneid'))

#corDf <- dplyr::left_join(corDf,pval,c('external_gene_name'='GeneID'))


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
#to look at PC-adjusted pseudotime, need to remove NAs to run the correlation
#dlpfc_pc <- subset(dlpfc, !is.na(dlpfc$Pseudotime_PCs))
#dlpfc2_pc <- dplyr::select(dlpfc_pc,dplyr::starts_with("ENSG"))
#corvec <- cor(dlpfc2,dlpfc$Pseudotime_unadj,method='spearman')
#corvec <- cor(dlpfc2_pc,dlpfc_pc$Pseudotime_PCs,method='spearman')
#corvec <- cor(dlpfc2_pc,dlpfc_pc$Pseudotime_RIN,method='spearman')
corvec <- cor(dlpfc2,dlpfc$Pseudotime_PMI,method='spearman')
corDfdlpfc <- data.frame(geneid=colnames(dlpfc2),cor=corvec,stringsAsFactors=F)
corDfdlpfc2 <- corDfdlpfc
corDfdlpfc2$external_gene_name <- Make.Gene.Symb(corDfdlpfc2$geneid)
corDfdlpfc2$cor <- NULL
corDfdlpfc2 <- dplyr::left_join(corDfdlpfc,corDfdlpfc2,by=c('geneid'))




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
mean(abs(corDf$cor))
mean(abs(corDf[corDf$external_gene_name %in% ad_gwas,]$cor))

corDfdlpfc <- corDfdlpfc2
mean(abs(corDfdlpfc$cor))
mean(abs(corDfdlpfc[corDfdlpfc$external_gene_name %in% ad_gwas,]$cor))


corDf$adGwas <- corDf$external_gene_name %in% ad_gwas
corDfdlpfc$adGwas <- corDfdlpfc$external_gene_name %in% ad_gwas
corDf$brainRegion <- 'TCX'
corDfdlpfc$brainRegion <- 'DLPFC'
corDfcombined <- rbind(corDf,corDfdlpfc)
colnames(corDfcombined)[4] <- 'LOADGWASGene'
g <- ggplot2::ggplot(corDfcombined,ggplot2::aes(x=brainRegion,y=cor,fill=LOADGWASGene))
g <- g + ggplot2::geom_boxplot()
g <- g + ggplot2::scale_fill_viridis_d()
#g <- g + ggplot2::geom_point(position=ggplot2::position_jitterdodge())
g <- g + ggplot2::labs(x = 'Brain Region',y='Correlation with PMI-adj pseudotime',fill='LOAD\nGWAS\nGene')

#g <- g + ggplot2::geom_point(position=ggplot2::position_jitterdodge())
g
# tiff(file='~/Desktop/MANUSCRIPT/figure2d.tiff',height=85,width=100,units='mm',res=300)
g
dev.off()

