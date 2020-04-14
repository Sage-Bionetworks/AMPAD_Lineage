source('MiscPreprocessing.R')
source('TCX_MonocleFunction_SetRootState.R')

#Loading data 
#synapse id of dat file: syn8466816

synapser::synLogin()

tcxCPMObj <- synapser::synGet('syn8466816')
Dat <- read.delim(tcxCPMObj$path,stringsAsFactors = F)


#synapse id of dat2 file: syn8466814
tcxCovObj <- synapser::synGet('syn8466814')
Dat2 <- read.delim(tcxCovObj$path,stringsAsFactors = F)
#add genetic principal components to the metadata file:
tcx_pcs <- read.csv("~/AMPAD_Lineage/LH_reviewer_response/TCX_principal_components_F.csv")
tcx_pcs$SampleID <- as.character(tcx_pcs$SampleID)
Dat2 <- dplyr::left_join(Dat2, tcx_pcs, by="SampleID")

#subsetting genes based on differential expression
#synapse id of DE file: syn8468023?
de_file <- synapser::synGet('syn8468023')
de_file2 <- synapser::synGet('syn18475579')

de1 <- data.table::fread(de_file$path,data.table=F)
de2 <- data.table::fread(de_file2$path,data.table=F)
de3 <- dplyr::filter(de1,Model=='Diagnosis',Comparison=='AD-CONTROL',Tissue.ref=='TCX')
#AMP_mods <-  read.csv('Data/TCX_DE.csv')
#AMP_mods <- data.frame(GeneID=de3$ensembl_gene_id,logPV= - log(de3$adj.P.Val)/log(10), stringsAsFactors=F)
AMP_mods <- de2[,-1]
In <- which(AMP_mods$logPV >= 1)
#change the DE p-value threshold to 0.01 (-logPV >= 2)
#In <- which(AMP_mods$logPV >= 2)
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

#Subsetting based on brain region
In_BR <- grep('TCX',Dat2$Tissue.Diagnosis)
DatNorm3 <- DatNorm2[,In_BR]
Dat3 <- Dat2[In_BR,]

#subsetting based on gender
Sex <- 'FEMALE'
In_S <- which(Dat3$Sex == Sex)
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
#temp2 <- cbind(Dat4,Cov)
temp2 <- Dat4

temp2$Diagnosis <- temp2$Tissue.SourceDiagnosis

temp2$Diagnosis[temp2$Tissue.SourceDiagnosis=='TCX.AD'] <- 'AD'
temp2$Diagnosis[temp2$Tissue.SourceDiagnosis=='TCX.CONTROL'] <- 'Control'
temp2$Diagnosis[temp2$Tissue.SourceDiagnosis=='TCX.PATH_AGE'] <- 'PA'
temp2$Diagnosis[temp2$Tissue.SourceDiagnosis=='TCX.PSP'] <- 'PSP'
#create a log PMI column (PMI is heavily right-skewed)
temp2$log_pmi <- log(temp2$PMI + 1)

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


#run monocle
rownames(temp) <- NULL
colnames(temp) <- NULL
#saveRDS(temp, file="TCX_countmatrix_Mono.rds")
#saveRDS(temp2, file="TCX_cell_metadata_Mono.rds")
#saveRDS(gene_short_name, file="TCX_gene_metadata_Mono.rds")
#read in data files to simply run monocle:
temp <- readRDS(file="TCX_countmatrix_Mono.rds")
temp2 <- readRDS(file="TCX_cell_metadata_Mono.rds")
gene_short_name <- readRDS(file="TCX_gene_metadata_Mono.rds")

# must unload synapser because causes multiple definitions of S4 vectors
detach("package:synapser", unload=TRUE)
unloadNamespace("PythonEmbedInR") 

library(monocle)
#Original results
MonRun <- RunMonocleTobit(temp, temp2, C_by = 'Pseudotime',gene_short_name = gene_short_name)
plot_cell_trajectory(MonRun, color_by = "Tissue.Diagnosis")

#need to specify the correct root for the unadjusted model and then adjust for each cofactor in different 
#models. Source code for monocle tobit functions called 
#below in 'MonocleFunction_SetRootState.R'

MonRun <- RunMonocleTobit_unadj(temp, temp2, C_by = 'Pseudotime',gene_short_name = gene_short_name)
plot_cell_trajectory(MonRun, color_by = "Tissue.Diagnosis")
#tiff(file='~/AMPAD_Lineage/LH_reviewer_response/TCX_celltraj_unadj.tiff',height=85,width=100,units='mm',res=300)
g<- plot_cell_trajectory(MonRun,color_by = "Diagnosis",show_branch_points=F,use_color_gradient = F,cell_size = 1)
g <- g + ggplot2::scale_color_viridis_d()
g <- g + ggplot2::labs(color="Diagnosis")
g
#dev.off()


MonRun_RIN <- RunMonocleTobit_RIN(temp, temp2, C_by = 'Pseudotime',gene_short_name = gene_short_name)
g<- plot_cell_trajectory(MonRun_RIN,color_by = "Diagnosis",show_branch_points=F,use_color_gradient = F,cell_size = 1)
g <- g + ggplot2::scale_color_viridis_d()
g <- g + ggplot2::labs(color="Diagnosis")
g

MonRun_PMI <- RunMonocleTobit_PMI(temp, temp2, C_by = 'Pseudotime',gene_short_name = gene_short_name)
g<- plot_cell_trajectory(MonRun_PMI,color_by = "Diagnosis",show_branch_points=F,use_color_gradient = F,cell_size = 1)
g <- g + ggplot2::scale_color_viridis_d()
g <- g + ggplot2::labs(color="Diagnosis")
g

MonRun_PCs <- RunMonocleTobit_PCs(temp, temp2, C_by = 'Pseudotime',gene_short_name = gene_short_name)
g<- plot_cell_trajectory(MonRun_PCs,color_by = "Diagnosis",show_branch_points=F,use_color_gradient = F,cell_size = 1)
g <- g + ggplot2::scale_color_viridis_d()
g <- g + ggplot2::labs(color="Diagnosis")
g

MonRun_ALL <- RunMonocleTobit_ALL(temp, temp2, C_by = 'Pseudotime',gene_short_name = gene_short_name)
g<- plot_cell_trajectory(MonRun_ALL,color_by = "Diagnosis",show_branch_points=F,use_color_gradient = F,cell_size = 1)
g <- g + ggplot2::scale_color_viridis_d()
g <- g + ggplot2::labs(color="Diagnosis")
g

#output pseudotimes in a dataframe similar to the existing pseudotimes in synapse (copied from TCX_BranchPhenotype.R)
#Make new dataframe with SampleID, State, Pseudotime
x <- list()


x$SampleID <- MonRun$SampleID
x$State_unadj <- MonRun$State
x$Pseudotime_unadj <- MonRun$Pseudotime
x$State_RIN <- MonRun_RIN$State
x$Pseudotime_RIN <- MonRun_RIN$Pseudotime
x$State_PMI <- MonRun_PMI$State
x$Pseudotime_PMI <- MonRun_PMI$Pseudotime
x$State_PCs <- MonRun_PCs$State
x$Pseudotime_PCs <- MonRun_PCs$Pseudotime
x$State_ALL <- MonRun_ALL$State
x$Pseudotime_ALL <- MonRun_ALL$Pseudotime
x <- as.data.frame(x)


write.csv(x,'~/AMPAD_Lineage/LH_reviewer_response/TCX_F_PStimes_adjusted.csv')


