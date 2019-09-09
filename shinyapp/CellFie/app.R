#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#
library(R.matlab)
library(ggplot2)
library(shiny)
library(openxlsx)
library(ggplot2)
library(gplots)
library(ggbiplot)

# Define UI for application that draws a histogram
ui <- fluidPage(
   
   # Application title
   titlePanel("Cellfie"),
   
   # Sidebar with a slider input for number of bins 
   sidebarLayout(
     sidebarPanel(
       h3("Upload data"),
       fileInput('expressionCSV','Upload Expression',accept = c('text/csv','text/comma-separated-values,text/plain','.csv')),
       fileInput('annotCSV','Upload Annotation',accept = c('text/csv','text/comma-separated-values,text/plain','.csv')),
       selectInput('model','Select Metabolic Model',rev(c('MT_iCHOv1_final.mat','MT_iHsa.mat','MT_iMM1415.mat','MT_inesMouseModel.mat',
                                                      'MT_iRno.mat','MT_quek14.mat','MT_recon_1.mat','MT_recon_2.mat','MT_recon_2_2.mat','MT_recon_2_2_entrez.mat'))),
       
       h3("Select Threshold Types"),
       selectInput('intype','Threshold Type',c('local','global')),
       selectInput('valuetype','Threshold Value Type',c('value','percentile')),
       selectInput('localtype','Threshold Value Type (local only)',c('minmaxmean','mean')),
       
       h3("Value Threshold"),
       sliderInput("valMin","Minimum Value:",min = 0,max = 1000,value = 25),
       sliderInput("valMax","Minimum Value:",min = 10,max = 1000,value = 75),
       
       h3("Percent Threshold"),
       sliderInput("percMin","Minimum Percent:",min = 0,max = 1,value = .25),
       sliderInput("percMax","Minimum Percent:",min = 0,max = 1,value = .75),
       
       actionButton('run',"Run CellFie")
       
     ),
      
      # Show a plot of the generated distribution
      mainPanel(
         tabsetPanel(
           tabPanel('Inputs',tableOutput('table'),tableOutput('table2'),plotOutput("expsshist")),
           tabPanel('Run',textOutput("runOut0"),textOutput("runOut")),
           tabPanel('Results',
                    sidebarPanel(
                      downloadButton("downloadData1","Download Score"),
                      downloadButton("downloadData2","Download Binary Score"),
                      downloadButton("downloadData3","Download Task Information")
                      ),
                    plotOutput('heatmap'),plotOutput('biplot')
                    )
         )
      )
   )
)

# Define server logic required to draw a histogram
server <- function(input, output) {
  print(getwd())

  #### load reactive variables
  
  # load gene/model
  model_list = read.xlsx('../../input/List_GenesInModels.xlsx')  
  model_genes = reactive({ unique(na.omit(model_list[[input$model]])) })
  
  # load expression & annotation data data
  expss = reactive({ req(input$expressionCSV); return( read.csv(input$expressionCSV$datapath,header = TRUE) ) })
  annot = reactive({ req(input$annotCSV); return( read.csv(input$annotCSV$datapath,header = TRUE) ) })
  
  #### generate reactive outputs for tab1: "Input"
  output$table <- renderTable({ return(head(expss()))  })
  output$table2 <- renderTable({ return(head(annot()))  })
  
  output$expsshist <- renderPlot({
    count <- sum(idx<-expss()[,1] %in% model_genes())
    hist(expssDistr<-log10(as.numeric(unlist(expss()[idx,-1]))), col = 'darkgray', border = 'white',
         main=paste(count,'of uploaded genes found in',input$model))
    
    #values = log10( ifelse(input$valuetype=='value',  c(input$valMin,input$valMax),
    #                quantile(expssDistr,c(input$percMin,input$percMax))) )
    #minmax = ifelse(input$intype=='global',  c(values[1],NA),
    #                ifelse(input$localtype=='mean',  c(log10(mean(expssDistr,na.rm=T)),NA),
    #                       c(values[1],values[2])  ))

    abline(v=log10(median(expssDistr,na.rm = T)),col='red')
    #abline(v=values[1],col='blue')
    #abline(v=values[2],col='blue')
  })
  
  ######## 
  # run cellfie
  unique_index = make.names(Sys.time()) ### for stable file saves, append to score*.csv
  dir.create(unique_index)
  
  output$runOut0 <- renderText({
    print('get some coffee, this will take a few minutes...')
  })
 
  textout = eventReactive(input$run,{
       cmd = paste('./application/run_execCellfie.sh v94/',
                   input$expressionCSV$datapath, ncol(expss())-1 , input$model , input$intype , input$valuetype ,
                   input$localtype, input$valMin,input$valMax,paste(unique_index),'.') 
      return( system( cmd , intern=T ))
   })
  
  output$runOut <- renderText({
    textout()
  })
  
  ########
  # results tab
  
  # load results
  #cellfieout = reactive({ req(f<-paste('cellfieout.mat')); return( read.csv(input$expressionCSV$datapath,header = TRUE) ) })
  score = reactive({ req(f<-file.path(unique_index,'score.csv')); return( read.csv(f,header = FALSE) ) })
  score_bin = reactive({ req(f<-file.path(unique_index,'score_binary.csv')); return( read.csv(f,header = FALSE) ) })
  tasks = reactive({ req(f<-file.path(unique_index,'taskInfo.csv')); return( read.csv(f,header = TRUE) ) })
  
  # plot results
  #output$biplot = renderPlot({
  #  dat = score()
  #  colnames(dat) = annot()[,1]
  #  rownames(dat) = tasks()[,2]
  #  print(head(dat))
  #  pr = prcomp(log(t(data.matrix(dat[rowMeans(dat)>2,]))+1e-3),center = TRUE,scale = TRUE)
  #  ggbiplot(pr,labels=colnames(dat),groups=annot()[,2],ellipse = TRUE)
    #heatmap.2(dat[rowMeans(dat)>0,],trace='none')
  #})
  # plot results
  output$heatmap = renderPlot({
    dat = score()
    colnames(dat) = annot()[,1]
    rownames(dat) = tasks()[,2]
    print(head(dat))
    heatmap.2(t(log(na.omit(data.matrix(dat[rowMeans(dat)>2,]))+1e-3)),mar=c(10,5),trace='none')
  })
  
  # download results
  output$downloadData1 = downloadHandler(filename = function(){paste('cellfie_score-',Sys.Date(),'.csv',sep='')}, content = function(file){ write.csv(score() ,file)})
  output$downloadData2 = downloadHandler(filename = function(){paste('cellfie_scorebinary-',Sys.Date(),'.csv',sep='')}, content = function(file){ write.csv(score_bin() ,file)})
  output$downloadData3 = downloadHandler(filename = function(){paste('cellfie_tasks-',Sys.Date(),'.csv',sep='')}, content = function(file){ write.csv(tasks() ,file)})
}

# Run the application 
shinyApp(ui = ui, server = server)

