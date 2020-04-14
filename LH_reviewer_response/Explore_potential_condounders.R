#Explore possible confounders (RIN, PMI, genetic principal components) by looking at 
#associations between the unadjusted pseudotime estimates (currently in synapse) with 
#each of these factors

#get the previously calculated pseudotimes & covariate files:
synapser::synLogin()
tcxLineageTimes <- synapser::synTableQuery("select * from syn17023721")$asDataFrame()
dlpfcLineageTimes <- synapser::synTableQuery("select * from syn17023795")$asDataFrame()
tcxCovObj <- synapser::synGet("syn8466814")
tcxCov <- data.table::fread(tcxCovObj$path,data.table=F)
dlpfcCovObj <- synapser::synGet("syn11024258")
dlpfcCov <- data.table::fread(dlpfcCovObj$path,data.table=F)
tcxLineageTimes <- tcxLineageTimes[,-c(1:3)]
dlpfcLineageTimes <- dlpfcLineageTimes[,-c(1:3)]
dlpfc <- dplyr::left_join(dlpfcLineageTimes,dlpfcCov,by='SampleID')
tcx <- dplyr::left_join(tcxLineageTimes,tcxCov,by='SampleID')


#Look at PMI & RIN distributions in dlpfc
summary(dlpfc$pmi)
hist(dlpfc$pmi)
summary(tcx$PMI)
hist(tcx$PMI)

summary(dlpfc$RINcontinuous)
hist(dlpfc$RINcontinuous)
summary(dlpfc$RINcontinuous2)
hist(dlpfc$RINcontinuous2) #use RINcontinuous
hist(tcx$RIN)
hist(tcx$RIN2)

#PMI is very right skewed. log-transform to correct for skew (add a 1 to avoid zeros)
dlpfc$log_pmi <- log(dlpfc$pmi + 1)
hist(dlpfc$log_pmi)  #this looks better
tcx$log_pmi <- log(tcx$PMI + 1)
hist(tcx$log_pmi)  #this looks better


#scale pseudotime:
dlpfc$pseudotime_sc <- scale(dlpfc$Pseudotime, center=F)
hist(dlpfc$pseudotime_sc)
tcx$pseudotime_sc <- scale(tcx$Pseudotime, center=F)
hist(tcx$pseudotime_sc)


#Run a linear regression between pmi & pseudotime
fit <- lm(pseudotime_sc~log_pmi, data=dlpfc)
summary(fit)
plot(fit$residuals)
library(ggplot2)
g <- ggplot(dlpfc, aes(x=pseudotime_sc, y=log_pmi)) + geom_point() + geom_smooth(method=lm) + labs(title="Pseudotime vs. PMI", x="Pseudotime (scaled)", y = "PMI (log)") + theme_classic()
#tiff(file='Laura_plots/scatter_pmiVSpseudo.tiff',height=85,width=100,units='mm',res=300)
g
#dev.off()
fit <- lm(pseudotime_sc~log_pmi, data=tcx)
summary(fit)
plot(fit$residuals)
g <- ggplot(tcx, aes(x=pseudotime_sc, y=log_pmi)) + geom_point() + geom_smooth(method=lm) + labs(title="Pseudotime vs. PMI", x="Pseudotime (scaled)", y = "PMI (log)") + theme_classic()
#tiff(file='Laura_plots/scatter_pmiVSpseudo.tiff',height=85,width=100,units='mm',res=300)
g
#dev.off()

#associations between pseudotime & rin: dlpfc
fit <- lm(pseudotime_sc~RINcontinuous, data=dlpfc)
summary(fit)
plot(fit$residuals)
g <- ggplot(dlpfc, aes(x=pseudotime_sc, y=RINcontinuous)) + geom_point() + geom_smooth(method=lm) + labs(title="Pseudotime vs. RIN", x="Pseudotime (scaled)", y = "RIN continuous") + theme_classic()
#tiff(file='Laura_plots/scatter_rinVSpseudo.tiff',height=85,width=100,units='mm',res=300)
g
#dev.off()

#RIN is highly negatively associated with pseudotime. is rin associated with case status?
ggplot(dlpfc, aes(x=Diagnosis, y=RINcontinuous)) + geom_boxplot()
dlpfc$diagnosis_casecontrol <- ifelse(dlpfc$Diagnosis=='AD', 1,
                                      ifelse(dlpfc$Diagnosis=='CONTROL',0,NA))
table(dlpfc$Diagnosis)
table(dlpfc$diagnosis_casecontrol)
fit <- lm(pseudotime_sc~diagnosis_casecontrol, data=dlpfc)
summary(fit)
fit <- lm(pseudotime_sc~diagnosis_casecontrol, data=dlpfc)
summary(fit)
summary(glm(diagnosis_casecontrol ~ RINcontinuous,dlpfc,family='binomial'))

#associations between pseudotime & rin: tcx
fit <- lm(pseudotime_sc~RINcontinuous, data=tcx)
summary(fit)
plot(fit$residuals)
g <- ggplot(tcx, aes(x=pseudotime_sc, y=RIN)) + geom_point() + geom_smooth(method=lm) + labs(title="Pseudotime vs. RIN", x="Pseudotime (scaled)", y = "RIN continuous") + theme_classic()
#tiff(file='Laura_plots/scatter_rinVSpseudo.tiff',height=85,width=100,units='mm',res=300)
g
#dev.off()


ggplot(tcx, aes(x=Tissue.SourceDiagnosis, y=RIN)) + geom_boxplot()
tcx$diagnosis_casecontrol <- ifelse(tcx$Tissue.Diagnosis=='TCX.AD', 1,
                                      ifelse(tcx$Tissue.Diagnosis=='TCX.CONTROL',0,NA))
table(tcx$Tissue.Diagnosis)
table(tcx$diagnosis_casecontrol)
fit <- lm(pseudotime_sc~diagnosis_casecontrol, data=tcx)
summary(fit)
summary(glm(diagnosis_casecontrol ~ RIN,tcx,family='binomial'))




######check for associations with PCs for genetic ancestry######
eigenvectors <- synapser::synGet("syn20820117")
class(eigenvectors)
eigenvectors2 <- data.table::fread(eigenvectors$path,data.table=F)
class(eigenvectors2)
View(eigenvectors2)
names(eigenvectors2) <- c("assay", "projid", "PC1", "PC2", "PC3", "PC4", "PC5", "PC6", "PC7", 
                          "PC8", "PC9", "PC10", "PC11", "PC12", "PC13", "PC14", "PC15", "PC16", "PC17", "PC18", "PC19", "PC20")
#get metadata and biospecimen ids from rosmap:
metadata <- synapser::synGet("syn3191087")
metadata_IDs <- read.csv(metadata$path)
metadata_IDs <- subset(metadata_IDs, select=c(projid, Study,individualID ))

#combine Study and projid
metadata_IDs$projid2 <- metadata_IDs$projid
metadata_IDs$projid<-NULL

metadata_IDs$projid = paste(metadata_IDs$Study,metadata_IDs$projid2)
metadata_IDs$projid <- gsub('\\s+', '', metadata_IDs$projid)

biospec1s <- synapser::synGet("syn21323366")
biospecs <- read.csv(biospec1s$path)
biospecsID <- subset(biospecs, select=c(individualID, specimenID))
allIDs <- merge(metadata_IDs, biospecsID)

names(allIDs)[names(allIDs) == "specimenID"] <- "SampleID"
rosmap_eigens2 <- merge(allIDs, dlpfc, by="SampleID")
rosmap_eigens2 <- merge(rosmap_eigens2, eigenvectors2, by="projid")
#save just the PCs for later analysis:
dlpfc_PCs <- subset(rosmap_eigens2, select=c(SampleID,PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10,PC11,PC12,PC13,PC14,PC15,PC16,PC17,PC18,PC19,PC20))
write.csv(dlpfc_PCs, file="DLPFC_principal_components_F.csv")
dlpfc2 <- merge(dlpfc, rosmap_eigens2)
dim(dlpfc2)

#SamplID in tcx matches projid in eigenvectors (minus the added '_TCX' in tcx df) 
tcx$projid <- gsub('_TCX', '', tcx$SampleID)
tcx2 <- merge(tcx, eigenvectors2, by="projid")
tcx_PCs <- subset(tcx2, select=c(SampleID,PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10,PC11,PC12,PC13,PC14,PC15,PC16,PC17,PC18,PC19,PC20))
write.csv(tcx_PCs, file="TCX_principal_components_F.csv")


summary(lm(pseudotime_sc~PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10, data=dlpfc2))
summary(lm(pseudotime_sc~PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10, data=tcx2))


