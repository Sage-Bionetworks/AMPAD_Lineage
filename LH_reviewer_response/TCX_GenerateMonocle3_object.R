
synapser::synLogin()

#MAYO_CBE_TCX_logCPM.tsv (filtered log counts)
tcxCPMObj <- synapser::synGet('syn8466816')
Dat <- read.delim(tcxCPMObj$path,stringsAsFactors = F)
rownames(Dat) <- Dat$ensembl_gene_id

#Covariates (cell_metadata) file for all patients
tcxCovObj <- synapser::synGet('syn8466814')
Dat2 <- read.delim(tcxCovObj$path,stringsAsFactors = F)

#subsetting genes based on differential expression
de_file <- synapser::synGet('syn8468023')
de_file2 <- synapser::synGet('syn18475579')

de1 <- data.table::fread(de_file$path,data.table=F)
de2 <- data.table::fread(de_file2$path,data.table=F)
de3 <- dplyr::filter(de1,Model=='Diagnosis',Comparison=='AD-CONTROL',Tissue.ref=='TCX')
#AMP_mods <-  read.csv('Data/TCX_DE.csv')
#AMP_mods <- data.frame(GeneID=de3$ensembl_gene_id,logPV= - log(de3$adj.P.Val)/log(10), stringsAsFactors=F)
AMP_mods <- de2[,-1]
In <- which(AMP_mods$logPV >= 1)
AMP_mods <- AMP_mods[In,]


#Dat column names need to match the row names of gene_metadata (ensemble IDs)
#Dat 2 row names need to match Dat column names (patient IDs)
#Need to create gene_metadata file with row names=ensemble ids, matched to one column (gene_short_name)
#gene_metadata <- subset(de3, select=c("ensembl_gene_id", "hgnc_symbol"))


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

ColNorm <- function(Dat){
  
  M = max(colSums(Dat))
  l <- length(colnames(Dat))
  
  for( i in 1:l){
    
    Dat[,i] = Dat[,i]*(M/sum(Dat[,i]))
    
  }
  
  return(Dat)
}

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
GeneNamesAD <- as.data.frame(GeneNamesAD)
gene_short_name <- as.data.frame(gene_short_name)
gene_metadata <- cbind(GeneNamesAD, gene_short_name)
DatNorm5<-DatNorm4
rownames(DatNorm5) <- gene_metadata$GeneNamesAD
counts <- Matrix(as.matrix(DatNorm5), sparse=TRUE)
rownames(gene_metadata) <- gene_metadata$GeneNamesAD
rownames(temp2) <- temp2$SampleID

#to save the components of the monocle object for later use:
saveRDS(counts, file="LH_TCX_processed_countmatrix.rds")
saveRDS(gene_metadata, file="LH_TCX_processed_gene_metadata.rds")
saveRDS(temp2, file="LH_TCX_processed_cell_metadata.rds")
saveRDS(cds_tcx, file="LH_CDS_TCX_Monocle3Object.rds")

#for monocle: expression matrix=counts, cell_metadata=temp2, gene_metadata=gene_metadata

# must first unload synapser because causes multiple definitions of S4Vectors
detach("package:synapser", unload=TRUE)
unloadNamespace("PythonEmbedInR") 
#run monocle and get Monocle Object: cds_tcx


cds_tcx <- new_cell_data_set(counts, cell_metadata=temp2, gene_metadata=gene_metadata)

cds_tcx <- preprocess_cds(cds_tcx, method="PCA", norm_method="size_only", num_dim=30)
cds_tcx <- reduce_dimension(cds_tcx)
plot_pc_variance_explained(cds_tcx)
cds_tcx <- cluster_cells(cds_tcx)
cds_tcx@clusters$UMAP$partitions[cds_tcx@clusters$UMAP$partitions == "2"] <- "1"


cds_tcx <- learn_graph(cds_tcx, use_partition=FALSE)
cds_tcx >- order_cells(cds_tcx, root_pr_nodes=get_earliest_principal_node(cds_tcx))

plot_cells(cds_tcx, xcolor_cells_by="pseudotime")
plot_cells(cds_tcx, color_cells_by="partition")