#Loading pre-computed lineage
load("~/Documents/DriverPrediction/LineageMisc/Data/TCX_Monocle.RData")
source('MonoclePlottingCustom.R')

dev.off()

#Diagnosis
MonRun$Diagnosis <- MonRun$Tissue.SourceDiagnosis
#Changing category names to match ROSMAP
MonRun$Diagnosis[MonRun$Diagnosis == 'TCX.AD'] <- 'AD' 
MonRun$Diagnosis[MonRun$Diagnosis == 'TCX.CONTROL'] <- 'CONTROL' 
MonRun$Diagnosis[MonRun$Diagnosis == 'TCX.PSP'] <- 'PSP' 
MonRun$Diagnosis[MonRun$Diagnosis == 'TCX.PATH_AGE'] <- 'PATH_AGE' 

jpeg('./PaperFigs/TCX_diag_new.jpg',width = 600, height = 350)
PlotMonocleDiscrete(MonRun@reducedDimS, MonRun@reducedDimK, './Data/TCX_graph.txt', 
                    clr = MonRun$Diagnosis, cp = NULL, sz1=0, sz2=0)
dev.off()

#State
MonRun$State2 <- MonRun$State
MonRun$State2[MonRun$State == 6] <- 4
MonRun$State2[MonRun$State == 7] <- 6
jpeg('./PaperFigs/TCX_state_new.jpg',width = 600, height = 350)
PlotMonocleDiscrete(MonRun@reducedDimS, MonRun@reducedDimK, './Data/TCX_graph.txt', 
                    clr = MonRun$State2, cp = NULL, sz1=0, sz2=0)
dev.off()



########################################################################################

#diagnosis boxplot
jpeg('./PaperFigs/TCX_diag_box_new.jpg',width = 600, height = 350)
p <- ggplot(MonRun@phenoData@data, aes(x=Diagnosis, y=Pseudotime,fill=Diagnosis)) + 
  geom_boxplot() + stat_summary(fun.y=mean, geom="point", shape=23, size=2)  +
  theme(axis.text=element_text(size=15), axis.title=element_text(size=15,face="bold"),
        legend.text=element_text(size=15))  
p
dev.off()

