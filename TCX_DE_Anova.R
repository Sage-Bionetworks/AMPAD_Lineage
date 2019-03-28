#Loading pre-computed lineage
load("~/Documents/DriverPrediction/LineageMisc/Data/TCX_Monocle.RData")
source('MonoclePlottingCustom.R')
source('MiscPreprocessing.R')
source('LineageFunctions.R')


#getting gene symbol
gene_short_name <- Make.Gene.Symb(GeneNamesAD)

MonRun2 = MonRun
ScaledDat = ScaleMatrix(MonRun2@assayData$exprs)
MonRun$State2 <- MonRun$State
MonRun$State2[MonRun$State == 6] <- 4
MonRun$State2[MonRun$State == 1] <- 6
MonRun$State2[MonRun$State == 7] <- 1

#pre-process data for ANOVA test
l2 <- list()
l2$gene_names <- gene_short_name
l2$p_2 <- rep(0,length(gene_short_name))
l2$p_3 <- rep(0,length(gene_short_name))
l2$p_4 <- rep(0,length(gene_short_name))
l2$p_5 <- rep(0,length(gene_short_name))
l2$p_6 <- rep(0,length(gene_short_name))

l2$d_2 <- rep(0,length(gene_short_name))
l2$d_3 <- rep(0,length(gene_short_name))
l2$d_4 <- rep(0,length(gene_short_name))
l2$d_5 <- rep(0,length(gene_short_name))
l2$d_6 <- rep(0,length(gene_short_name))


for (i in 1:length(gene_short_name)){
  l <- list()
  l$x <- as.vector(t(temp[i,]))
  l$s <- as.character(MonRun$State2)
  
  df <- as.data.frame(l)
  res.aov <- aov(x ~ s, data = df)
  tk <- TukeyHSD(res.aov)
  
  l2$p_2[i] <- tk$s[1,4]
  l2$p_3[i] <- tk$s[2,4]
  l2$p_4[i] <- tk$s[3,4]
  l2$p_5[i] <- tk$s[4,4]
  l2$p_6[i] <- tk$s[5,4]
  
  l2$d_2[i] <- tk$s[1,1]
  l2$d_3[i] <- tk$s[2,1]
  l2$d_4[i] <- tk$s[3,1]
  l2$d_5[i] <- tk$s[4,1]
  l2$d_6[i] <- tk$s[5,1]
  
}

#Save data
df2 <- as.data.frame(l2)
saveRDS(df2,'Data/TCX_F_pv1_Anova_DE_1_6_switched.rds')