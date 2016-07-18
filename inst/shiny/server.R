library(shiny)
package.phenology <- require('phenology')

# Define server logic required to draw a histogram
shinyServer(function(input, output) {
  
  # Expression that generates a histogram. The expression is
  # wrapped in a call to renderPlot to indicate that:
  #
  #  1) It is "reactive" and therefore should re-execute automatically
  #     when inputs change
  #  2) Its output type is a plot
  
  
  output$distPlot <- renderPlot({
    
    if (!package.phenology) {
      par(mar=c(0, 0, 0, 0))
      plot(x=c(0, 1), y=c(0, 1), axes=FALSE, 
           xaxt="n", yaxt="n", main="", 
           xlab = "", ylab = "", 
           xaxs="i", yaxs="i", type="n")
      text(x = 0.5, y=0.5, labels = "The phenology package is not installed!", 
           col="red", cex = 2)
      text(x = 0.5, y=0.3, labels = "Contact your site administrator", 
           col="green", cex = 2)
    }
    
    if (is.null(input$file1) & package.phenology) {
      par(mar=c(0, 0, 0, 0))
      plot(x=c(0, 1), y=c(0, 1), axes=FALSE, 
           xaxt="n", yaxt="n", main="", 
           xlab = "", ylab = "", 
           xaxs="i", yaxs="i", type="n")
      text(x = 0.5, y=0.6, labels = "First, choose the model with the buttons down and then", 
           col="red", cex = 2)
      text(x = 0.5, y=0.5, labels = "load dataset using the button at the left.", 
           col="red", cex = 2)
      text(x = 0.5, y=0.4, labels = "When the model will be fitted, results will be shown here.", 
           col="green", cex = 2)
      text(x = 0.5, y=0.3, labels = "It can take one minute or more. Be patient!", 
           col="green", cex = 2)
      text(x = 0.5, y=0.2, labels = "It is better (more rapid) to load the dataset after choosing the model and the header option.", 
           col="black", cex = 1)
      
    }
    if (!is.null(input$file1) & package.phenology) {
      
      inFile <- input$file1
      dest <- inFile$datapath
      
      # dest <- "/Users/marc/Documents/Espace_de_travail_R/package_phenology/LUTH\ 2002/Yalimapo_2002_GS.txt"
      
      table <- read.delim(dest, header = input$Header, stringsAsFactors = FALSE)
      
      lby <- ncol(table)
      
      if (lby == 3) {
        nm <- levels(as.factor(table[, 3]))
        if (length(nm)>1) nm <- NULL
      } else {
        nm <- inFile$name
      }
      
      df2 <- table
  
      Formated <- add_phenology(add=df2, silent=TRUE, name=nm)

      if (input$Min == 2) {
        pfixed <- c(Flat=0)
      } else {
        pfixed <- c(Flat = 0, Min = 0)
      }
      Par <- par_init(data=Formated, parametersfixed=pfixed)
      if (input$Length == 1) Par <- LBLE_to_L(Par)
      if (input$Min == 3) {
        Par <- MinBMinE_to_Min(Par)
        Par <- toggle_Min_PMin(Par)
      }
      
      zLogic <- ifelse(input$Zero==1, TRUE, FALSE)

      result <- fit_phenology(data=Formated, parametersfit = Par,
                              parametersfixed=pfixed, trace=0, silent=TRUE, 
                              zero_counts = zLogic)
      
      # x <- installed.packages()
      
      x <- plot(result, series="all", moon=FALSE, progressbar = FALSE,
                growlnotify = FALSE)
      
      output$resultsInfo <- renderPrint({print(x)})
      output$resultsRaw <- renderPrint({print(result)})
    }
  })
})