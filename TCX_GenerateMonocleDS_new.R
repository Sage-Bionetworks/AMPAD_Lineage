source('MiscPreprocessing.R')

#Loading data 
#synapse id of dat file: syn8466816

synapser::synLogin()

tcxCPMObj <- synapser::synGet('syn8466816')
Dat <- read.delim(tcxCPMObj$path,stringsAsFactors = F)

#synapse id of dat2 file: syn8466814
tcxCovObj <- synapser::synGet('syn8466814')
Dat2 <- read.delim(tcxCovObj$path,stringsAsFactors = F)

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


#convert to gene symbol
gene_short_name <- Make.Gene.Symb(GeneNamesAD)

#run monocle
rownames(temp) <- NULL
colnames(temp) <- NULL
MonRun <- RunMonocleTobit(temp, temp2, C_by = 'Pseudotime',gene_short_name = gene_short_name)
#save.image(file='temp41720.rda')


set.seed(2)
res <- vector('list',ncol(temp))
corvec <- rep(0,ncol(temp))
for (i in 1:ncol(temp)){
  print(i)
  #ind1<-sample(1:ncol(temp),ncol(temp),replace = T)
  res[[i]] <- RunMonocleTobit(temp[,-i], temp2[-i,], C_by = 'Pseudotime',gene_short_name = gene_short_name)
  corvec[i] <- abs(cor(MonRun$Pseudotime[-i],res[[i]]$Pseudotime,method='spearman'))
  hist(corvec[1:i])
}

for (i in 1:ncol(temp)){
  corvec[i] <- (cor(MonRun$Pseudotime[-i],res[[i]]$Pseudotime,method='spearman'))
  #hist(corvec[1:i])
}
cordf <- data.frame(cor=corvec,brainRegion='TCX',stringsAsFactors = F)
g <- ggplot2::ggplot(cordf,ggplot2::aes(cor,fill='red'))
g <- g + ggplot2::geom_density(adjust=.1)
g <- g + cowplot::theme_cowplot()
g

write.csv(cordf,file='tcx_jacknifed_correlations.csv',quote=F,row.names=F)

save(res,file='tcx_jacknifed_lineages.rda')


set.seed(47)
tcx_pca <- svd(scale(t(temp)))$u[,1:2]
system.time(tcx_tsne <- tsne::tsne(scale(t(temp))))
tcx_umap <- umap::umap(scale(t(temp)))

tcx_dimred_df <- data.frame(Pseudotime=MonRun$Pseudotime,
                              PCA1=tcx_pca[,1],
                              PCA2=tcx_pca[,2],
                              tSNE1=tcx_tsne[,1],
                              tSNE2=tcx_tsne[,2],
                              UMAP1=tcx_umap$layout[,1],
                              UMAP2=tcx_umap$layout[,2],
                              stringsAsFactors=F)


tiff(file='~/Desktop/figure_tcx_pseudotime_pca1.tiff',height=85,width=100,units='mm',res=300)
g <- ggpubr::ggscatter(tcx_dimred_df,
                       x='Pseudotime',
                       y='PCA1',
                       add='reg.line',
                       add.params = list(color='blue',fill='lightgray'),
                       conf.int=TRUE)
g <- g + ggpubr::stat_cor(method='pearson',label.x=0,label.y=-0.2)
g
dev.off()

tiff(file='~/Desktop/figure_tcx_pseudotime_pca2.tiff',height=85,width=100,units='mm',res=300)
g <- ggpubr::ggscatter(tcx_dimred_df,
                       x='Pseudotime',
                       y='PCA2',
                       add='reg.line',
                       add.params = list(color='blue',fill='lightgray'),
                       conf.int=TRUE)
g <- g + ggpubr::stat_cor(method='pearson',label.x=0,label.y=-0.35)
g
dev.off()

tiff(file='~/Desktop/figure_tcx_pseudotime_tsne1.tiff',height=85,width=100,units='mm',res=300)
g <- ggpubr::ggscatter(tcx_dimred_df,
                       x='Pseudotime',
                       y='tSNE1',
                       add='reg.line',
                       add.params = list(color='blue',fill='lightgray'),
                       conf.int=TRUE)
g <- g + ggpubr::stat_cor(method='pearson',label.x=0,label.y=-10)
g
dev.off()

tiff(file='~/Desktop/figure_tcx_pseudotime_tsne2.tiff',height=85,width=100,units='mm',res=300)
g <- ggpubr::ggscatter(tcx_dimred_df,
                       x='Pseudotime',
                       y='tSNE2',
                       add='reg.line',
                       add.params = list(color='blue',fill='lightgray'),
                       conf.int=TRUE)
g <- g + ggpubr::stat_cor(method='pearson',label.x=0,label.y=-10)
g
dev.off()


tiff(file='~/Desktop/figure_tcx_pseudotime_umap1.tiff',height=85,width=100,units='mm',res=300)
g <- ggpubr::ggscatter(tcx_dimred_df,
                       x='Pseudotime',
                       y='UMAP1',
                       add='reg.line',
                       add.params = list(color='blue',fill='lightgray'),
                       conf.int=TRUE)
g <- g + ggpubr::stat_cor(method='pearson',label.x=0,label.y=-2.5)
g
dev.off()

tiff(file='~/Desktop/figure_tcx_pseudotime_umap2.tiff',height=85,width=100,units='mm',res=300)
g <- ggpubr::ggscatter(tcx_dimred_df,
                       x='Pseudotime',
                       y='UMAP2',
                       add='reg.line',
                       add.params = list(color='blue',fill='lightgray'),
                       conf.int=TRUE)
g <- g + ggpubr::stat_cor(method='pearson',label.x=0,label.y=3)
g
dev.off()




#plot_cell_trajectory(MonRun, color_by = "Tissue.Diagnosis")

#MonRun$apoe_genotype <- factor(MonRun$apoe_genotype,levels=c(0,1,2))
fxn3 <- function(x){
  if(x=='TCX.0'){
    return(0)
  }else if (x=='TCX.1'){
    return(1)
  }else if (x=='TCX.2'){
    return(2)
  }
}
MonRun$apoe <- sapply(MonRun$Tissue.APOE4,fxn3)
MonRun$apoe <- factor(MonRun$apoe,levels=c(0,1,2))

tiff(file='~/Desktop/MANUSCRIPT/figureS6b.tiff',height=85,width=100,units='mm',res=300)
g<- plot_cell_trajectory(MonRun,color_by = "apoe",show_branch_points=F,use_color_gradient = F,cell_size = 0.5)
g <- g + ggplot2::scale_color_viridis_d()
g <- g + ggplot2::labs(color="APOE e4 Dosage")
g
dev.off()


tiff(file='~/Desktop/MANUSCRIPT/figure2aTCX.tiff',height=85,width=100,units='mm',res=300)
g<- plot_cell_trajectory(MonRun,color_by = "Diagnosis",show_branch_points=F,use_color_gradient = F,cell_size = 0.5)
g <- g + ggplot2::scale_color_viridis_d()
g
dev.off()

MonRun$State2 <- MonRun$State
MonRun$State2[MonRun$State == 6] <- 4
MonRun$State2[MonRun$State == 7] <- 6

MonRun$State2[MonRun$State2 == 1] <- 7
MonRun$State2[MonRun$State2 == 6] <- 1
MonRun$State2[MonRun$State2 == 7] <- 6

tiff(file='~/Desktop/MANUSCRIPT/figure4a.tiff',height=85,width=100,units='mm',res=300)
g<- plot_cell_trajectory(MonRun,color_by = "State2",show_branch_points=F,use_color_gradient = F,cell_size = 0.5)
g <- g + ggplot2::scale_color_viridis_d()
g <- g + ggplot2::labs(color="Lineage\nState")
g
dev.off()

foo <- synapser::synTableQuery("SELECT * FROM syn18409904")$asDataFrame()
foo <- foo[,-c(1:2)]
rownames(foo) <- foo$col
foo <- foo[,-c(1)]
colnames(foo) <- c('State 4','State 5','State 6','State 1','State 3','State 2')
tiff(file='~/Desktop/MANUSCRIPT/figure4b.tiff',height=85,width=100,units='mm',res=300)
pheatmap::pheatmap(foo,show_rownames=F,color = viridis::viridis(100))
dev.off()


apoeDf <- data.frame(apoe=MonRun$Tissue.APOE4=='TCX.1' | MonRun$Tissue.APOE4=='TCX.2',
                     State=MonRun$State2,
                     stringsAsFactors=F)

summary(glm(apoe ~ State,apoeDf,family='binomial'))


