#Loading pre-computed lineage
load("~/Documents/DriverPrediction/LineageMisc/Data/TCX_Monocle.RData")
source('MonoclePlottingCustom.R')
source('MiscPreprocessing.R')
source('LineageFunctions.R')

#Get diagnositic status variable and fit smooth spline 
Y <- (MonRun$Tissue.SourceDiagnosis=='TCX.AD') + 0.0 
ft <- smooth.spline(MonRun@phenoData@data$Pseudotime, Y, df = 3)
Y_pred <- predict(ft,MonRun@phenoData@data$Pseudotime)$y

Srtd <- sort(MonRun@phenoData@data$Pseudotime,index.return = T)
Srtd <- Srtd$ix

#plotting
plot(MonRun@phenoData@data$Pseudotime[Srtd], Y[Srtd], col = 'Red', 
     xlab = 'Inferred stage', ylab = 'Fraction of AD', pch=16,cex=1.5)
lines(MonRun@phenoData@data$Pseudotime[Srtd], Y_pred[Srtd], lwd = 3, col = 'Blue')