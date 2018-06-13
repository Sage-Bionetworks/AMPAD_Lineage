library(shiny)
source('E:/SageDocs/PredictingDriverGenes/LineageMisc/LineageFunctions.R')
#library(datasets)
setwd('E:/SageDocs/PredictingDriverGenes/LineageMisc/')
source('MiscPreprocessing.R')
Dat <- read.delim('E:/SageDocs/PredictingDriverGenes/LineageMisc/MAYO_CBE_TCX_logCPM.tsv',stringsAsFactors = F)
Dat2 <- read.delim('E:/SageDocs/PredictingDriverGenes/LineageMisc/MAYO_CBE_TCX_Covariates.tsv',stringsAsFactors = F)
Cov <- read.csv('E:/SageDocs/PredictingDriverGenes/LineageMisc/mayo_igap_snps.csv',stringsAsFactors = F)
Cov[,2:22] <- round(Cov[,2:22])
AMP_mods <-  read.csv('E:/SageDocs/PredictingDriverGenes/LineageMisc/TCX_DE.csv')


#head(TempEnr)

# Define server logic required to plot various variables against mpg
shinyServer(function(input, output) {
  
  # Compute the forumla text in a reactive expression since it is 
  # shared by the output$caption and output$mpgPlot expressions
  Sex_I <- reactive(input$Sex)
  pv_I <- reactive(input$pv)
  lab_I <- reactive(input$lab)
  

  tmp <- reactive({
    if (is.null(Sex_I())){
      return(NULL)
    }
    
    if (is.null(pv_I())){
      return(NULL)
    }
    

    Sex <- Sex_I()
    pv <- pv_I()

    
    In <- which(AMP_mods$logPV >= pv)
    AMP_mods <- AMP_mods[In,]
    
    
    #Normalize all columns 
    
    GeneNames <- Dat$ensembl_gene_id
    GeneNamesAD <- AMP_mods$GeneID
    
    Names <- colnames(Dat)
    
    for (i in 1:length(Names)){
      
      Names[i] <- substring(Names[i],2)
      
    }
    
    
    colnames(Dat) <- Names
    cNames <- Dat2$SampleID
    l <- length(Names)
    
    #deleting columns not in the covariate list
    temp <- rep(T,l)
    for (i in 1:l){
      if (!(Names[i] %in% cNames)){
        temp[i] <- F
      }
    }
    
    In <- which(temp)
    #print(temp)
    Dat <- Dat[,In]
    
    #deleting extra rows in covariate list
    Names <- Names[In]
    l <- length(cNames)
    temp <- rep(T,l)
    for (i in 1:l){
      if (!(cNames[i] %in% Names)){
        temp[i] <- F
      }
    }
    In <- which(temp)
    Dat2 <- Dat2[In,]
    
    
    DatNorm <- ColNorm(Dat)
    In <- which(GeneNames %in% GeneNamesAD)
    DatNorm2 <- DatNorm[In,]
    
    In_BR <- grep('TCX',Dat2$Tissue.Diagnosis)
    DatNorm3 <- DatNorm2[,In_BR]
    Dat3 <- Dat2[In_BR,]
    
    #Keeping only female data 
    In_S <- which(Dat3$Sex == Sex)
    DatNorm4 <- DatNorm3[,In_S]
    Dat4 <- Dat3[In_S,]
    
    In_cov <- which(Cov$ID %in% Dat4$Donor_ID)
    Cov <- Cov[In_cov,]
    In_cov <- c()
    for(i in 1:length(Cov$ID)){
      temp <- which(Dat4$Donor_ID==Cov$ID[i])
      In_cov <- c(In_cov,temp[1])
    }
    Cov <- Cov[In_cov,]
    
    for (i in 23:26){
      Cov[,i] <- (Cov[,i] - min(Cov[,i]))/(max(Cov[,i])-min(Cov[,i]))
    }
    
    temp <- DatNorm4
    temp2 <- cbind(Dat4,Cov)
    rownames(temp) <- NULL
    colnames(temp) <- NULL
    MonRun <- RunMonocleTobit(temp, temp2, C_by = 'endoScore')


    tmp2 <- list()
    tmp2$MonRun <- MonRun

    tmp2
  })
  
  tmp3 <- reactive({

    if (is.null(lab_I())){
      return(NULL)
    }
    
    lab <- lab_I()
    tmp2 <- list()
    tmp2$lab <- lab
    
    tmp2
  })
  
  

  output$mpgPlot <- renderPlot({
    if (is.null(tmp())){
      return(NULL)
    }
    
    if (is.null(tmp3())){
      return(NULL)
    }
    
    t<- tmp()
    t2 <- tmp3()
    
    plot_cell_trajectory(t$MonRun, color_by = t2$lab)
    
    #if(!is.null(tmp())){
    #  t <- tmp()
    #  plot_cell_trajectory(t$MonRun, color_by = t$lab)

    #}
  })
})