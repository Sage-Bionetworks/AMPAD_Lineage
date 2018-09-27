#synapseClient::synapseLogin()
#rosmapClinicalObj <- synapseClient::synGet('syn3191087')
#rosmapUncensoredAgesObj <- synapseClient::synGet('syn7116000')
#rosmapIdMapObj <- synapseClient::synGet('syn3382527')
#rosmapCogDecline1Obj <- synapseClient::synGet('syn6182375')
#rosmapCogDecline2Obj <- synapseClient::synGet('syn6182376')
#rosmapMetanetworkObj <- synapseClient::synGet('syn8268669')

synapser::synLogin(email = 'mukhes3@uw.edu', password = 'Nirvana99')
rosmapClinicalObj <- synapser::synGet('syn3191087')
rosmapUncensoredAgesObj <- synapser::synGet('syn7116000')
rosmapIdMapObj <- synapser::synGet('syn3382527')
rosmapCogDecline1Obj <- synapser::synGet('syn6182375')
rosmapCogDecline2Obj <- synapser::synGet('syn6182376')
rosmapMetanetworkObj <- synapser::synGet('syn8268669')


####thanneer's code for fixing ids
#https://github.com/th1vairam/Brain_Reg_Net/blob/ad12b544d6d6c23be2bf94eed53a6b4f75d154d1/code/Rmd/ROSMAP_REPROCESSED.Rmd

####read in clinical file
rosmapClinical <- data.table::fread(rosmapClinicalObj$path,
                                    data.table=F)

rosmapClinical <- dplyr::select(rosmapClinical,-V1,-age_death,-age_first_ad_dx,-age_at_visit_max)
#View(rosmapClinical)


####read in uncensored age file
rosmapUncensoredAges <- data.table::fread(rosmapUncensoredAgesObj$path,
                                          data.table=F)
#View(rosmapUncensoredAges)

####merge uncensored ages with clinical file
rosmapClinical <- dplyr::left_join(rosmapClinical,rosmapUncensoredAges,by='projid')
#View(rosmapClinical)

###cognitive decline
rosmapCogDec1 <- data.table::fread(rosmapCogDecline1Obj$path,
                                   data.table=F)

rosmapCogDec1 <- dplyr::select(rosmapCogDec1,-Sample)

rosmapCogDec2 <- data.table::fread(rosmapCogDecline2Obj$path,
                                   data.table=F)

mungeIds <- function(x){
  #ros or map
  isROS<-grep('ROS',x)
  isMAP<-grep('MAP',x)
  if(length(isROS)>0){
    rosId<-strsplit(x,'ROS')[[1]][2]
    return(as.numeric(rosId))
  }else if(length(isMAP)>0){
    mapId <- strsplit(x,'MAP')[[1]][2]
    return(as.numeric(mapId))
  }else{
    return(as.numeric(x))
  }
}
projids2 <- sapply(rosmapCogDec2$CollaboratorParticipantId,mungeIds)
rosmapCogDec2$ProjectID <- projids2
rosmapCogDec2 <- dplyr::select(rosmapCogDec2,ProjectID,cogn_global_slope)
rosmapCogDec <- rbind(rosmapCogDec1,rosmapCogDec2)

#View(rosmapCogDec1)
rosmapClinical <- dplyr::left_join(rosmapClinical,
                                   rosmapCogDec,
                                   by=c("projid"="ProjectID"))
rosmapClinical$apoe_genotype <- as.factor(rosmapClinical$apoe_genotype)
rosmapClinical$ceradsc <- as.factor(rosmapClinical$ceradsc)
rosmapClinical$braaksc <- as.factor(rosmapClinical$braaksc)
rosmapClinical$cogdx <- as.factor(rosmapClinical$cogdx)



