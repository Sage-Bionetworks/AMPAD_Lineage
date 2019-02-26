load('Data/DLPFC_Monocle.RData')
MonRun$State2 <- MonRun$State
MonRun$State2[MonRun$State == 4] <- 3
MonRun$State2[MonRun$State == 7] <- 3
MonRun$State2[MonRun$State == 6] <- 4
MonRun$State2[MonRun$State == 5] <- 4
MonRun$State2[MonRun$State == 8] <- 5
MonRun$State2[MonRun$State == 9] <- 6


l <- list()
l$S1 <- (MonRun$State2==1) + 0.0 
l$S2 <- (MonRun$State2==2) + 0.0 
l$S3 <- (MonRun$State2==3) + 0.0 
l$S4 <- (MonRun$State2==4) + 0.0 
l$S5 <- (MonRun$State2==5) + 0.0 
l$S6 <- (MonRun$State2==6) + 0.0 
l$P <- MonRun$Pseudotime


l$Y <- (MonRun$APOE4 != '0') + 0.0 

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

library(MASS)
l2 <- list()
l2$Y <- as.factor(MonRun$cogdx.1)
l2$P <- MonRun$Pseudotime 
df2 <- as.data.frame(l2)

m <- polr(Y ~ P , data = df2, Hess=TRUE)

## view a summary of the model
summary(m)

## store table
(ctable <- coef(summary(m)))

## calculate and store p values
p <- pnorm(abs(ctable[, "t value"]), lower.tail = FALSE) * 2

## combined table
(ctable <- cbind(ctable, "p value" = p))
