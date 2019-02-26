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
  
  if (dim(df2)[1]>1){
    m <- max(df2$age_at_visit)
    M <- min(df2$age_at_visit)
    In3 <- which(df2$age_at_visit == m)
    In4 <- which(df2$age_at_visit == M)
    cogn_ep[i] <- (df2$cogn_ep[In3] - df2$cogn_ep[In4])/(m-M)
    cogn_po[i] <- (df2$cogn_po[In3] - df2$cogn_po[In4])/(m-M)
    cogn_ps[i] <- (df2$cogn_ps[In3] - df2$cogn_ps[In4])/(m-M)
    cogn_se[i] <- (df2$cogn_se[In3] - df2$cogn_se[In4])/(m-M)
    cogn_wo[i] <- (df2$cogn_wo[In3] - df2$cogn_wo[In4])/(m-M)
    cogn_global[i] <- (df2$cogn_global[In3] - df2$cogn_global[In3])/(m-M)
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
MonRun$r_depres <- r_depres


dev.off()

#trajectory plotting

jpeg('./PaperFigs/PIL_traits/cogn_ep_traj_del.jpg')
plot_cell_trajectory(MonRun, color_by = 'cogn_ep')
dev.off()

jpeg('./PaperFigs/PIL_traits/cogn_po_traj_del.jpg')
plot_cell_trajectory(MonRun, color_by = 'cogn_po')
dev.off()

jpeg('./PaperFigs/PIL_traits/cogn_ps_traj_del.jpg')
plot_cell_trajectory(MonRun, color_by = 'cogn_ps')
dev.off()

jpeg('./PaperFigs/PIL_traits/cogn_se_traj_del.jpg')
plot_cell_trajectory(MonRun, color_by = 'cogn_se')
dev.off()

jpeg('./PaperFigs/PIL_traits/cogn_wo_traj_del.jpg')
plot_cell_trajectory(MonRun, color_by = 'cogn_wo')
dev.off()

jpeg('./PaperFigs/PIL_traits/cogn_global_traj_del.jpg')
plot_cell_trajectory(MonRun, color_by = 'cogn_global')
dev.off()