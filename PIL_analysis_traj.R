#Loading ROSMAP monocle file
load("~/Documents/DriverPrediction/LineageMisc/Data/DLPFC_Monocle.RData")

#loading CSV file containing PiL data 
df <- read.csv('./Data/ROSMAP_PiL_dep_cog.csv')
#df$study <- NULL
#df$fu_year <- NULL
#df$scaled_to <- NULL
#df$dcfdx <- NULL
#df$cesdsum <- NULL


#loading PROJid to SampleID mapping file
Map <- read.csv('./Data/ROSMAP_IDkey.csv')

#Get values for different clinical covariates at last valid visit 
l <- length(MonRun$SampleID)
cogn_ep <- rep(NaN,l)
cogn_po <- rep(NaN,l)
cogn_ps <- rep(NaN,l)
cogn_se <- rep(NaN,l)
cogn_wo <- rep(NaN,l)
cogn_global <- rep(NaN,l)
purpose_total <- rep(NaN,l)
r_depres <- rep(NaN,l)


for(i in 1:l){
  In <- which(Map$rnaseq_id == MonRun$SampleID[i])
  pID = Map$projid[In]
  pID = pID[!is.na(pID)][1]
  
  In2 <- which(df$projid==pID)
  df2 <- df[In2,]
  df2 <- df2[!is.na(df2$cogn_ep)&!is.na(df2$cogn_po) & !is.na(df2$cogn_ps) & !is.na(df2$cogn_se) 
             & !is.na(df2$cogn_wo)&!is.na(df2$cogn_global)&!is.na(df2$r_depres),]
  
  if (dim(df2)[1]>0){
    m <- max(df2$age_at_visit)
    In3 <- which(df2$age_at_visit == m)
    cogn_ep[i] <- df2$cogn_ep[In3]
    cogn_po[i] <- df2$cogn_po[In3]
    cogn_ps[i] <- df2$cogn_ps[In3]
    cogn_se[i] <- df2$cogn_se[In3]
    cogn_wo[i] <- df2$cogn_wo[In3]
    cogn_global[i] <- df2$cogn_global[In3]
    purpose_total[i] <- df2$purpose_total[In3]
    r_depres[i] <- df2$r_depres[In3]
    
  } 
    
}
  

MonRun$cogn_ep <- cogn_ep
MonRun$cogn_po <- cogn_po
MonRun$cogn_ps <- cogn_ps
MonRun$cogn_se <- cogn_se
MonRun$cogn_wo <- cogn_wo
MonRun$global <- cogn_global
MonRun$purpose_total <- purpose_total
MonRun$r_depres <- r_depres


GeneratePlots <-  function(x,y,z,ObjName){
  par(mfrow=c(1,2))
  
  fit <- lm(y ~ z)
  y_res <- y - (fit$coefficients[1] + fit$coefficients[2]*z)
  

  fit1 <- lm(y~x)
  pv <- formatC(as.numeric(summary(fit1)$coefficients[,4][2]), format = "e", digits = 2)
  rs <- formatC(summary(fit1)$r.squared, digits = 3)
  txt <- paste(c('R^2 = ',rs,', p-V = ',pv),sep='',collapse = '')
    
  #making plot1
  plot(x,y,pch = 16, cex = 1.3, col = "blue",main=txt,xlab = 'Pseudotime',ylab=ObjName)
  abline(fit1,lwd = 3, col = 'red')
  
  
  fit2 <- lm(y_res ~ x)
  pv <- formatC(as.numeric(summary(fit2)$coefficients[,4][2]), format = "e", digits = 2)
  rs <- formatC(summary(fit2)$r.squared, digits = 3)
  txt <- paste(c('R^2 = ',rs,', p-V = ',pv),sep='',collapse = '')
  txt2 <- paste(c(ObjName,'-corrected'),sep='',collapse = '')
  
  #making plot2  
  plot(x,y_res,pch = 16, cex = 1.3, col = "blue",main=txt,xlab = 'Pseudotime',ylab=txt2)
  abline(fit2,lwd = 3, col = 'red')
  
  #making monocle plot
  #plot_cell_trajectory(MonRun, color_by = ObjName)
  
  return(0)
    
}

P <- MonRun$Pseudotime
B <- as.numeric(MonRun$braaksc)

dev.off()
jpeg('./PaperFigs/cogn_ep.jpg')
GeneratePlots(P,MonRun$cogn_ep,B,'cogn_ep')
dev.off()

jpeg('./PaperFigs/cogn_do.jpg')
GeneratePlots(P,MonRun$cogn_po,B,'cogn_po')
dev.off()

jpeg('./PaperFigs/cogn_ps.jpg')
GeneratePlots(P,MonRun$cogn_ps,B,'cogn_ps')
dev.off()

jpeg('./PaperFigs/cogn_se.jpg')
GeneratePlots(P,MonRun$cogn_se,B,'cogn_se')
dev.off()

jpeg('./PaperFigs/cogn_wo.jpg')
GeneratePlots(P,MonRun$cogn_wo,B,'cogn_wo')
dev.off()

jpeg('./PaperFigs/cogn_global.jpg')
GeneratePlots(P,MonRun$global,B,'cogn_global')
dev.off()

#trajectory plotting

jpeg('./PaperFigs/cogn_ep_traj.jpg')
plot_cell_trajectory(MonRun, color_by = 'cogn_ep')
dev.off()

jpeg('./PaperFigs/cogn_po_traj.jpg')
plot_cell_trajectory(MonRun, color_by = 'cogn_po')
dev.off()

jpeg('./PaperFigs/cogn_ps_traj.jpg')
plot_cell_trajectory(MonRun, color_by = 'cogn_ps')
dev.off()

jpeg('./PaperFigs/cogn_se_traj.jpg')
plot_cell_trajectory(MonRun, color_by = 'cogn_se')
dev.off()

jpeg('./PaperFigs/cogn_wo_traj.jpg')
plot_cell_trajectory(MonRun, color_by = 'cogn_wo')
dev.off()

jpeg('./PaperFigs/cogn_global_traj.jpg')
plot_cell_trajectory(MonRun, color_by = 'cogn_global')
dev.off()


MonRun$State2 <- MonRun$State
MonRun$State2[MonRun$State == 4] <- 3
MonRun$State2[MonRun$State == 7] <- 3
MonRun$State2[MonRun$State == 6] <- 4
MonRun$State2[MonRun$State == 5] <- 4
MonRun$State2[MonRun$State == 8] <- 5
MonRun$State2[MonRun$State == 9] <- 6

#Making boxplots 
p <- ggplot(MonRun@phenoData@data, aes(x=State2, y=cogn_ep,fill=State2)) + geom_boxplot()
p + stat_summary(fun.y=mean, geom="point", shape=23, size=2)

p <- ggplot(MonRun@phenoData@data, aes(x=State2, y=cogn_global,fill=State2)) + geom_boxplot()
p + stat_summary(fun.y=mean, geom="point", shape=23, size=2)

p <- ggplot(MonRun@phenoData@data, aes(x=State2, y=cogn_ps,fill=State2)) + geom_boxplot()
p + stat_summary(fun.y=mean, geom="point", shape=23, size=2)

p <- ggplot(MonRun@phenoData@data, aes(x=State2, y=cogn_wo,fill=State2)) + geom_boxplot()
p + stat_summary(fun.y=mean, geom="point", shape=23, size=2)

p <- ggplot(MonRun@phenoData@data, aes(x=State2, y=cogn_wo,fill=State2)) + geom_boxplot()
p + stat_summary(fun.y=mean, geom="point", shape=23, size=2)

p <- ggplot(MonRun@phenoData@data, aes(x=State2, y=purpose_total,fill=State2)) + geom_boxplot()
p + stat_summary(fun.y=mean, geom="point", shape=23, size=2)
