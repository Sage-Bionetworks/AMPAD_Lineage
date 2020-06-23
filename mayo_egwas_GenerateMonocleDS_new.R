source('MiscPreprocessing.R')

#Loading data
synapser::synLogin()

mayoObj <- synapser::synGet('syn3617054')
Dat <- data.table::fread(mayoObj$path,data.table=F)
Dat <- dplyr::select(Dat,-FID)

#synapse id of dat2 file: syn8466814
mayoCovObj <- synapser::synGet('syn3617056')
Dat2 <- data.table::fread(mayoCovObj$path,data.table=F)
Dat2 <- dplyr::select(Dat2,-FID)

fullDat <- dplyr::left_join(Dat,Dat2)
fullDatFemale <- dplyr::filter(fullDat,Sex==1)
DatExpr <- dplyr::select(fullDatFemale,dplyr::starts_with('ILMN'))
DatCov <- dplyr::select(fullDatFemale,-dplyr::starts_with('ILMN'))
X <- dplyr::select(DatCov,plate2,plate3,plate4,plate5)
X$interc <- 1
X <- data.matrix(X)
betaPlate <- solve(t(X)%*%X)%*%t(X)%*%data.matrix(DatExpr)
DatExprAdj <- DatExpr-X%*%betaPlate
# advec <- sapply(DatFemale$DiseaseStatus,function(x) if(x=='Control'){return(0)}else{return(1)})

system.time(tvalues <- apply(DatExprAdj,2,function(x,y){utilityFunctions::fastlm(y,x)},DatCov$Dxn))
pvalues <- pt(abs(tvalues),184,lower.tail=F)*2

pvalDf <- data.frame(probeID=names(pvalues),
                     pval=pvalues,
                     adj.pval=p.adjust(pvalues,method='fdr'),
                     stringsAsFactors=F)
keep_probes <- pvalDf$probeID[pvalDf$adj.pval<=0.1]

DatExprFilt <- DatExprAdj[,keep_probes]

temp <- t(DatExprFilt)
temp2 <- DatCov

colnames(temp) <- NULL
rownames(temp) <- NULL
source('LineageFunctions.R')

ba<-ecdf(MonRun$Pseudotime)
MonRun <- RunMonocleTobit((temp), temp2, C_by = 'Pseudotime',gene_short_name = keep_probes)
MonRun$Dxn<- as.factor(MonRun$Dxn)
MonRun$Dxn2 <- sapply(MonRun$Dxn,function(x) if(x==0){return('Control')}else{return('AD')})
MonRun$Dxn2 <- as.factor(MonRun$Dxn2)
MonRun$resistant<-!(MonRun$Pseudotime>quantile(ba,.8) & MonRun$Dxn==0)
MonRun$resistant2<-sapply(MonRun$resistant,function(x) if(x==0){return('Resistant')}else{return('Not Resistant')})
MonRun$resistant2 <- factor(MonRun$resistant2,levels=c('Resistant','Not Resistant'))
MonRun$E4dose <- as.factor(MonRun$E4dose)

tiff(file='~/Desktop/pseuodtime_mayoegwas.tiff',height=85,width=100,units='mm',res=300)
g<- plot_cell_trajectory(MonRun,color_by = "Dxn2",show_branch_points=F,use_color_gradient = F,cell_size = 0.5)
g <- g + ggplot2::scale_color_viridis_d()
g <- g + ggplot2::labs(color="Disease\nState")
g
dev.off()

tiff(file='~/Desktop/pseuodtime_mayoegwas_resist.tiff',height=85,width=100,units='mm',res=300)
g<- plot_cell_trajectory(MonRun,color_by = "resistant2",show_branch_points=F,use_color_gradient = F,cell_size = 0.5)
g <- g + ggplot2::scale_color_viridis_d()
g <- g + ggplot2::labs(color="Resistant\nState")
g
dev.off()

tiff(file='~/Desktop/pseuodtime_mayoegwas_apoe.tiff',height=85,width=100,units='mm',res=300)
g<- plot_cell_trajectory(MonRun,color_by = "E4dose",show_branch_points=F,use_color_gradient = F,cell_size = 0.5)
g <- g + ggplot2::scale_color_viridis_d()
g <- g + ggplot2::labs(color="APOE E4\nDose")
g
dev.off()


mrdf <- data.frame(Disease_State=MonRun$Dxn2,
                   Pseudotime=MonRun$Pseudotime,
                   Diagnosis=MonRun$Dxn,
                   apoe=MonRun$E4dose)
g <- ggplot2::ggplot(mrdf,ggplot2::aes(x=Disease_State,
                                             y=Pseudotime,
                                             color=Disease_State))
g <- g + ggplot2::geom_boxplot()
g <- g + ggplot2::geom_point(position=ggplot2::position_jitterdodge())
#g <- g + ggplot2::scale_color_viridis_d()
g <- g + ggplot2::scale_color_manual(values=viridis::viridis(3)[1:2])
g <- g + ggplot2::xlab('Disease State')
g <- g + ggplot2::labs(color="Disease\nState")

tiff(file='~/Desktop/mayo_egwas_boxplot.tiff',height=85,width=100,units='mm',res=300)
g
dev.off()

mrdf$Pseudotimescaled<- scale(mrdf$Pseudotime)
summary(glm(Diagnosis ~ Pseudotimescaled,mrdf,family='binomial'))
mrdf2 <- dplyr::filter(mrdf,apoe!=-9)
mrdf2$apoe <- factor(mrdf2$apoe,levels=c('0','1','2'))
res4<-(MASS::polr(apoe ~ Pseudotimescaled,mrdf2))
cat('p-value: ',pt(abs(summary(res4)$coef[1,3]),181,lower.tail=F)*2,'\n')


df1 <- data.frame(pt=MonRun$Pseudotime,dxn=MonRun$Dxn,stringsAsFactors=F)
View(df1)
ba<-ecdf(MonRun$Pseudotime)


df1$resistant<-as.numeric(MonRun$Pseudotime>quantile(ba,.8) & MonRun$Dxn==0)
system.time(tvalues2 <- apply(DatExprAdj,2,function(x,y){utilityFunctions::fastlm(y,x)},df1$resistant))
pvalues2 <- pt(abs(tvalues2),184,lower.tail=F)*2


convert_illumina<-function(ensemblIds){
  
  library(biomaRt)
  ensembl=biomaRt::useMart('ENSEMBL_MART_ENSEMBL',
                           dataset = 'hsapiens_gene_ensembl',
                           host='uswest.ensembl.org')
  
  genes<-biomaRt::getBM(attributes = c('illumina_humanwg_6_v3','external_gene_name'),
                        filters='illumina_humanwg_6_v3',
                        values=ensemblIds,
                        mart=ensembl)
  return(genes)
}

mappingTab <- convert_illumina(colnames(DatExprAdj))

resistPvalDf <- data.frame(illumina_humanwg_6_v3=names(pvalues2),tval=tvalues2,pval=pvalues2,p.adj=p.adjust(pvalues2,method='fdr'),stringsAsFactors = F)

resistPvalDf <- dplyr::left_join(resistPvalDf,mappingTab)

foo <- synapser::synTableQuery("SELECT * FROM syn18409904")$asDataFrame()
foo <- foo[,-c(1:2)]
rownames(foo) <- foo$col
foo <- foo[,-c(1)]
colnames(foo) <- c('State 4','State 5','State 6','State 1','State 3','State 2')
respheat<-pheatmap::pheatmap(foo,show_rownames=F,color = viridis::viridis(100))

geneClusters <- cutree(respheat$tree_row,h = 2.0)
names(geneClusters) <- rownames(foo)

dFun <- function(i,gca){
  return(names(gca[gca==i]))
}

geneClList <- lapply(1:4,dFun,geneClusters)
names(geneClList) <- paste0('Cluster',1:4)
#fob <- data.table::fread('../lineage/go_linclust.tsv',data.table=F)

resistPvalDfup <- dplyr::filter(resistPvalDf,tval>0 &p.adj<=0.05)
resistPvalDfdown <- dplyr::filter(resistPvalDf,tval<0 &p.adj<=0.05)

geneClList$upInMayoeGWASResistant <- unique(resistPvalDfup$external_gene_name)
geneClList$downInMayoeGWASResistant <- unique(resistPvalDfdown$external_gene_name)

geneClDf<-utilityFunctions::list2df(geneClList)
geneClDf2 <- utilityFunctions::list2df(geneClList[-c(5,6)])
row.names(geneClDf2) <- geneClDf2$value
pheatmap::pheatmap(foo,show_rownames=F,color = viridis::viridis(100),annotation_row = dplyr::select(geneClDf2,key))
#foo2 <- dplyr::select(comb_anova2,gene_names,summary)
geneClDf$Presence <- 1
foo3 <- tidyr::pivot_wider(geneClDf,
                           id_cols = "value",
                           names_from = "key",
                           values_from = "Presence")
foo3[is.na(foo3)] <- 0
foo3 <- data.frame(foo3,stringsAsFactors=F)
colnames(foo3)[c(6,7)] <- c('Mayo eGWAS resistant up','Mayo eGWAS resistant down')

tiff(filename = '~/Desktop/upset_2.tiff', height = 6, width = 7,units='in',pointsize=1,res=300)
UpSetR::upset(foo3,nsets = 6)
dev.off()

utilityFunctions::fisherWrapper(geneClList$upInMayoeGWASResistant,geneClList$Cluster1,intersect(mappingTab$external_gene_name,unlist(geneClList[-c(5,6)])))

utilityFunctions::fisherWrapper(geneClList$upInMayoeGWASResistant,geneClList$Cluster2,intersect(mappingTab$external_gene_name,unlist(geneClList[-c(5,6)])))

utilityFunctions::fisherWrapper(geneClList$upInMayoeGWASResistant,geneClList$Cluster3,intersect(mappingTab$external_gene_name,unlist(geneClList[-c(5,6)])))

utilityFunctions::fisherWrapper(geneClList$upInMayoeGWASResistant,geneClList$Cluster4,intersect(mappingTab$external_gene_name,unlist(geneClList[-c(5,6)])))

utilityFunctions::fisherWrapper(geneClList$downInMayoeGWASResistant,geneClList$Cluster1,intersect(mappingTab$external_gene_name,unlist(geneClList[-c(5,6)])))

utilityFunctions::fisherWrapper(geneClList$upInMayoeGWASResistant,geneClList$Cluster1,intersect(mappingTab$external_gene_name,unlist(geneClList[-c(5,6)])))

utilityFunctions::fisherWrapper(geneClList$downInMayoeGWASResistant,geneClList$Cluster2,intersect(mappingTab$external_gene_name,unlist(geneClList[-c(5,6)])))

utilityFunctions::fisherWrapper(geneClList$downInMayoeGWASResistant,geneClList$Cluster3,intersect(mappingTab$external_gene_name,unlist(geneClList[-c(5,6)])))

utilityFunctions::fisherWrapper(geneClList$downInMayoeGWASResistant,geneClList$Cluster4,intersect(mappingTab$external_gene_name,unlist(geneClList[-c(5,6)])))

