library(shiny)

Sex <- c('MALE','FEMALE')
pv <- c(1:6)
lab <- c('Source','Sex','Tissue.SourceDiagnosis','Tissue.Diagnosis','Tissue.APOE4','AgeAtDeath','Tissue.Diagnosis.Sex','ID','rs6656401','rs6733839','rs10948363','rs11771145','rs9331896','rs983392','rs10792832','rs4147929','rs3865444','rs9271192','rs28834970','rs11218343','rs10498633','rs8093731','rs35349669','rs190982','rs2718058','rs1476679','rs10838725','rs17125944','rs7274581','riskScore','endoScore','immuneScore','neuroScore')


# Define UI for miles per gallon application
shinyUI(pageWithSidebar(
  
  # Application title
  headerPanel('Patient lineage visualization'),
  
  # Sidebar with controls to select the variable to plot against mpg
  # and to specify whether outliers should be included
  sidebarPanel(
    selectInput('pv', 'Choose a DE p-value cut-off:',
                choices = pv),
    selectInput('Sex', 'Choose a gender:',
                choices = Sex),
    selectInput('lab', 'Color lineage by:',
                choices = lab)
  ),
  
  # Show the caption and plot of the requested variable against mpg
  mainPanel(
    plotOutput('mpgPlot')
  )
))