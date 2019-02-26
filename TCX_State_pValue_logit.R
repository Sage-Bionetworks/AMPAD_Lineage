load('Data/TCX_Monocle.RData')
MonRun$State2 <- MonRun$State
MonRun$State2[MonRun$State == 6] <- 4
MonRun$State2[MonRun$State == 7] <- 6


l <- list()
l$S1 <- (MonRun$State2==1) + 0.0 
l$S2 <- (MonRun$State2==2) + 0.0 
l$S3 <- (MonRun$State2==3) + 0.0 
l$S4 <- (MonRun$State2==4) + 0.0 
l$S5 <- (MonRun$State2==5) + 0.0 
l$S6 <- (MonRun$State2==6) + 0.0 
l$P <- MonRun$Pseudotime


l$Y <- (MonRun$Tissue.APOE4 != 'TCX.0') + 0.0 

df <- as.data.frame(l)

mylogit <- glm(Y ~ S1 + S2 + S3 + S4 + S5 + S6 + P -1, data = df, family = binomial(link = "logit"))
summary(mylogit)
predicted <- predict(mylogit, df, type="response")
MonRun$Predicted <- predicted
plot_cell_trajectory(MonRun, color_by = 'Predicted')
library(ROCR)
ROCRpred <- prediction(predicted, l$Y)
ROCRperf <- performance(ROCRpred, 'tpr','fpr')
plot(ROCRperf, colorize = TRUE, text.adj = c(-0.2,1.7))