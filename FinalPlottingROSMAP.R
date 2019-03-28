#Loading pre-computed lineage
load("~/Documents/DriverPrediction/LineageMisc/Data/DLPFC_Monocle.RData")
source('MonoclePlottingCustom.R')

dev.off()

#Diagnosis
jpeg('./PaperFigs/DLPFC_diag_new.jpg',width = 600, height = 350)
PlotMonocleDiscrete(MonRun@reducedDimS, MonRun@reducedDimK, './Data/DLPFC_graph.txt', 
                      clr = MonRun$Diagnosis, cp = NULL, sz1=0, sz2=0)
dev.off()


#Braak

jpeg('./PaperFigs/DLPFC_braaksc_new.jpg',width = 600, height = 350)
PlotMonocleDiscrete(MonRun@reducedDimS, MonRun@reducedDimK, './Data/DLPFC_graph.txt', 
                      clr = MonRun$braaksc, cp = 'RdBu', sz1=0, sz2=0)
dev.off()

#cerad

jpeg('./PaperFigs/DLPFC_ceradsc_new.jpg',width = 600, height = 350)
PlotMonocleDiscrete(MonRun@reducedDimS, MonRun@reducedDimK, './Data/DLPFC_graph.txt', 
                      clr = MonRun$ceradsc, cp = 'RdBu', sz1=0, sz2=0)
dev.off()

#cogdx
jpeg('./PaperFigs/DLPFC_cogdx1_new.jpg',width = 600, height = 350)
PlotMonocleDiscrete(MonRun@reducedDimS, MonRun@reducedDimK, './Data/DLPFC_graph.txt', 
                    clr = MonRun$cogdx.1, cp = 'RdBu', sz1=0, sz2=0)
dev.off()

########################################################################################

#diagnosis boxplot
jpeg('./PaperFigs/DLPFC_diag_box_new.jpg',width = 600, height = 350)
p <- ggplot(MonRun@phenoData@data, aes(x=Diagnosis, y=Pseudotime,fill=Diagnosis)) + 
  geom_boxplot() + stat_summary(fun.y=mean, geom="point", shape=23, size=2)  +
  theme(axis.text=element_text(size=15), axis.title=element_text(size=15,face="bold"),
        legend.text=element_text(size=15))  
p
dev.off()

#braak boxplot
jpeg('./PaperFigs/DLPFC_braaksc_box_new.jpg',width = 600, height = 350)
p <- ggplot(MonRun@phenoData@data, aes(x=braaksc, y=Pseudotime,fill=braaksc)) + 
  geom_boxplot() + stat_summary(fun.y=mean, geom="point", shape=23, size=2)  +
  theme(axis.text=element_text(size=15), axis.title=element_text(size=15,face="bold"),
        legend.text=element_text(size=15)) + scale_fill_brewer(palette="RdBu") 
p
dev.off()

#cerad boxplot
jpeg('./PaperFigs/DLPFC_ceradsc_box_new.jpg',width = 600, height = 350)
p <- ggplot(MonRun@phenoData@data, aes(x=ceradsc, y=Pseudotime,fill=ceradsc)) + 
  geom_boxplot() + stat_summary(fun.y=mean, geom="point", shape=23, size=2)  +
  theme(axis.text=element_text(size=15), axis.title=element_text(size=15,face="bold"),
        legend.text=element_text(size=15)) + scale_fill_brewer(palette="RdBu")
p
dev.off()

#cogdx boxplot
jpeg('./PaperFigs/DLPFC_cogdx_box_new.jpg',width = 600, height = 350)
p <- ggplot(MonRun@phenoData@data, aes(x=cogdx.1, y=Pseudotime,fill=cogdx.1)) + 
  geom_boxplot() + stat_summary(fun.y=mean, geom="point", shape=23, size=2)  +
  theme(axis.text=element_text(size=15), axis.title=element_text(size=15,face="bold"),
        legend.text=element_text(size=15))  + scale_fill_brewer(palette="RdBu")
p
dev.off()



