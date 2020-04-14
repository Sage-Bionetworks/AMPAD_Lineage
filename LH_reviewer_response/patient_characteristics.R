#collect covariate info for TCX and DLPFC data for patient characteristics tables

tcxCovObj <- synapser::synGet('syn8466814')
mayo <- read.delim(tcxCovObj$path,stringsAsFactors = F)
dlpfcCovObj <- synapser::synGet('syn11024258')
rosmap1 <- read.delim(dlpfcCovObj$path,stringsAsFactors = F)
rosmapObj <- synapser::synGet('syn3191087')
rosmap2 <- data.table::fread(rosmapObj$path,data.table=F)

#Subsetting mayo data based on brain region
In_BR <- grep('TCX',mayo$Tissue.Diagnosis)
mayo <- mayo[In_BR,]

#removing bad batches in rosmap
rosmap1 <- subset(rosmap1, rosmap1$Batch<7)


#add in braak score & cerad score to rosmap data (need to synchronize Sample IDs & join)
rosmapIdObj <- synapser::synGet('syn3382527')
rosmapId <- data.table::fread(rosmapIdObj$path,data.table=F)
rosmapId <- dplyr::select(rosmapId,projid,rnaseq_id)
rosmapRNAid<-dplyr::left_join(rosmapId,rosmap2)
#remove duplicate rows
rosmapRNAid <- unique(rosmapRNAid)
#rename apoe_genotype column 
names(rosmapRNAid)[names(rosmapRNAid) == "apoe_genotype"] <- "apoe_complete"
rosmapRNAid2 <- subset(rosmapRNAid, select=c(rnaseq_id,braaksc,ceradsc,apoe_complete,cts_mmse30_first_ad_dx, cts_mmse30_lv))
names(rosmapRNAid2)[names(rosmapRNAid2) == "rnaseq_id"] <- "SampleID"

rosmap<-dplyr::left_join(rosmap1,rosmapRNAid2, by="SampleID")

#counts of males and females
table(mayo$Sex)
table(rosmap$msex)
#diagnosis data
table(mayo$Tissue.SourceDiagnosis, mayo$Sex)
table(rosmap$Diagnosis, rosmap$msex)
#braak score
table(rosmap$braaksc, rosmap$msex)
table(rosmap$ceradsc, rosmap$msex)
table(rosmap$cogdx, rosmap$msex)

table(mayo$Tissue.APOE4, mayo$Sex)
table(rosmap$apoe_genotype, rosmap$msex)

tapply(mayo$RIN, mayo$Sex, summary)
tapply(rosmap$RINcontinuous, rosmap$msex, summary)
