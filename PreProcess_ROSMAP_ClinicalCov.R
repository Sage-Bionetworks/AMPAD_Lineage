Cl_dat <- read.csv('ROSMAP_clinical.csv', stringsAsFactors = F)
Mp <- read.csv('ROSMAP_IDkey.csv', stringsAsFactors = F)

Mp <- Mp[which(Mp$rnaseq_id!=''),]
In_cl <- c()

for( i in 1:length(Mp$rnaseq_id)){
  
  In_cl <- c(In_cl,which(Cl_dat$projid==Mp$projid[i]))
  
}

Mp$projid <- NULL
Mp_new <- cbind(Mp,Cl_dat[In_cl,])

In_mp <- c()

for(i in 1:length(Dat2$SampleID)){
  In_mp <- c(In_mp,which(Mp_new$rnaseq_id==Dat2$SampleID[i])[1])
  
}

Dat2_new <- cbind(Dat2,Mp_new[In_mp,])

Dat2_new$age_first_ad_dx[Dat2_new$age_first_ad_dx==''] = 0
Dat2_new$age_first_ad_dx[Dat2_new$age_first_ad_dx=='90+'] = 0
Dat2_new$age_first_ad_dx <- as.numeric(Dat2_new$age_first_ad_dx)
Dat2_new$age_death[Dat2_new$age_death=='90+'] = 0
Dat2_new$age_death[Dat2_new$age_first_ad_dx==0] = 0 
Dat2_new$DeltaD <- abs(Dat2_new$age_death - Dat2_new$age_first_ad_dx)
write.csv(Dat2_new,'ROSMAP_DLPFC_Covariates3.csv')





