library(shiny)
package.phenology <- require('phenology')

load(file="babel.Rdata")

# Define server logic required to draw a histogram
shinyServer(function(input, output, session) {
  
  session$onSessionEnded(stopApp)
  
  output$Titre_UI <- renderUI({
    h1(babel["Titre" , input$language], align = "center")
  })
  output$Objet_UI <- renderUI({
    p(strong(babel["Objet" , input$language]))
  })
  output$WebVersion_UI <- renderUI({
    p(babel["WebVersion" , input$language], 
        a(babel["Ici" , input$language]
          , href="https://cran.r-project.org/package=phenology"
          , target="_blank"))
  })
  output$Developpeur_UI <- renderUI({
    p(babel["Developpeur" , input$language]
      , a("Marc Girondot"
          , href="https://max2.ese.u-psud.fr/epc/conservation/index.html"
          , target="_blank"))
  })
  
  output$Publication_UI <- renderUI({
    p(babel["Publication" , input$language],
      a(
        img(src="PDF.png", height=35, width=35), 
        "Girondot, M. 2010. Estimating density of animals during migratory waves: 
        application to marine turtles at nesting site. Endangered Species Research, 
        12, 85-105.", 
        href="http://www.int-res.com/articles/esr2010/12/n012p095.pdf",
        target="_blank"))
  })

  output$Liste_UI <- renderUI({
    h3(babel["Liste" , input$language])
  })
  
  output$RecommandeSWOT_UI <- renderUI({
    p(babel["RecommandeSWOT" , input$language]
      , a(babel["Norme" , input$language], 
          href=babel["NormePDF" , input$language],
          target="_blank")
      , babel["SWOT" , input$language])
  })
  output$RecommandeRASTOMA_UI <- renderUI({
    p(babel["RecommandeRASTOMA" , input$language])
  })
  
  output$ChargeData_UI <- renderUI({
    h4(babel["ChargeData" , input$language])
  })
  output$TitreCol_UI <- renderUI({
    checkboxInput("Header", babel["TitreCol" , input$language], value = FALSE, width = NULL)
  })
  output$NoteFormat_UI <- renderUI({
    helpText(babel["NoteFormat" , input$language])
  })
  
  output$ChargeExemple_UI <- renderUI({
    helpText(babel["ChargeExemple" , input$language])
  })
  
  
  output$Options_UI <- renderUI({
    h4(babel["Options" , input$language])
  })
  
  output$Min0_UI <- renderUI({  
    radioButtons("Min", babel["Min0" , input$language], list(yes=1, no=2), selected=2, inline = TRUE)
  })
  output$Sym_UI <- renderUI({  
    radioButtons("Length", babel["Sym" , input$language], list(yes=1, no=2), selected=2, inline = TRUE)
  })
  output$Report0_UI <- renderUI({  
    radioButtons("Zero", babel["Report0" , input$language], list(yes=1, no=2), selected=1, inline = TRUE)
  })
  output$FormatDate_UI <- renderUI({  
    textInput("FormatDate", babel["FormatDate" , input$language], value = "", width = NULL, placeholder = "%d/%m/%Y")
  })
  output$NoteDate_UI <- renderUI({  
    helpText(babel["NoteDate" , input$language])
  })
  output$MoisRef_UI <- renderUI({  
    numericInput("MonthRef", babel["MoisRef" , input$language], value = NULL, min=1, max=12, step=1, width = NULL)
  })
  output$NoteMoisRef_UI <- renderUI({  
    helpText(babel["NoteMoisRef" , input$language])
  })
  output$SerieGraphique_UI <- renderUI({  
    numericInput("Plot", babel["SerieGraphique" , input$language], value = 1, min=1, step=1, width = NULL)
  })
  output$NoteGraph_UI <- renderUI({  
    helpText(babel["NoteGraph" , input$language])
  })
  
  output$Synthese1_UI <- renderUI({
    h4(babel["Synthese1" , input$language])
  })
  # output$Synthese2_UI <- renderUI({
  #   h4(babel["Synthese2" , input$language])
  # })
  output$NoteStockage_UI <- renderUI({  
    helpText(babel["NoteStockage" , input$language])
  })  
  

  outplotfn <- eventReactive(eventExpr=input$goButton, ignoreNULL = FALSE, valueExpr={
    iMonthRef <- input$MonthRef
    ilanguage <- input$language
    if (is.null(iMonthRef)) {
      
      par(mar=c(0, 0, 0, 0))
      plot(x=c(0, 1), y=c(0, 1), axes=FALSE,
           xaxt="n", yaxt="n", main="",
           xlab = "", ylab = "",
           xaxs="i", yaxs="i", type="n")
      text(x = 0.5, y=0.6, labels = babel["Info1" , ilanguage],
           col="red", cex = 1.6)
      text(x = 0.5, y=0.5, labels = babel["Info2" , ilanguage],
           col="red", cex = 1.6)
      text(x = 0.5, y=0.4, labels = babel["Info3" , ilanguage],
           col="green", cex = 1.6)
      text(x = 0.5, y=0.3, labels = babel["Info4" , ilanguage],
           col="green", cex = 1.6)
      
    } else {
    

    ifile1 <- input$file1
    iHeader <- input$Header
    iMonthRef <- input$MonthRef
    iFormatDate <- input$FormatDate
    iMin <- input$Min
    iLength <- input$Length
    iZero <- input$Zero
    
    
    
#    ilanguage <- 1
#    ifile1 <- NULL
#    iHeader <- 1
#    iMonthRef <- 1
#    iFormatDate <- NULL
#    iMin <- 1
#    iLength <- 1
#    iZero <- 1
    


    if (!is.null(ifile1)) {
      inFile <- ifile1
      dest <- inFile$datapath
      
      # dest <- "/Users/marc/Documents/Espace_de_travail_R/package_phenology/LUTH\ 2002/Yalimapo_2002_GS.txt"
      
      table <- read.delim(dest, header = iHeader, stringsAsFactors = FALSE)
      if (ncol(table)==1) table <- read.csv(dest, header = iHeader, stringsAsFactors = FALSE)
      if (ncol(table)==1) table <- read.csv2(dest, header = iHeader, stringsAsFactors = FALSE)
      
    } else {
      table <- Gratiot
      table[, 1] <- as.character(table[, 1])
      
      inFile <- list(name = "Cayenne 2001 Leatherbacks")
    }
    
    lby <- ncol(table)
    
    
    
    if (lby == 1) {
      output$resultsInfo <- renderText({"Error in the data file"})
    } else {
      
      if (lby == 3) {
        nm <- levels(as.factor(table[, 3]))
        if (length(nm)>1) {
          nm <- NULL
          table[, 3] <- enc2native(table[, 3])
        } else {
          nm <- enc2native(nm)
        }
      } else {
        nm <- inFile$name
        nm <- enc2native(nm)
      }
    
    
    
#     df2 <- table
      zmonthRef <- 1
    if (is.numeric(iMonthRef)) zmonthRef <- iMonthRef
    if (iFormatDate == "") zFormatDate <- "%d/%m/%Y" else zFormatDate <- iFormatDate
    
    output$table <- renderTable({
      table
    })
    
    zLogic <- ifelse(iZero==1, TRUE, FALSE)
    
    
    Formated <- add_phenology(add=table, silent=TRUE, name=nm,
                              format=zFormatDate,
                              month_ref = zmonthRef, 
                              ZeroCounts.default = zLogic)
    
    pfixed <- NULL
    if (max(sapply(Formated, function(xxx) max(xxx[, "nombre"])))<10) pfixed <- c(Theta=+Inf)
    
    if (iMin == 2) {
      pfixed <- c(pfixed, Flat=0)
    } else {
      pfixed <- c(pfixed, Flat = 0, Min = 0)
    }
    Par <- par_init(data=Formated, fixed.parameters=pfixed)
    if (iLength == 1) Par <- LBLE_to_L(Par)
    if (iMin == 3) {
      Par <- MinBMinE_to_Min(Par)
      Par <- toggle_Min_PMin(Par)
    }
    
    result <- fit_phenology(data=Formated, fitted.parameters = Par,
                            fixed.parameters=pfixed,
                            silent=TRUE,
                            control = list(trace = 0, REPORT =
                                             100, maxit = 500), 
                            hessian = TRUE)
    
    
    x <- summary(result, series="all", print=FALSE)
    
    plot(result, series= max(as.numeric(input$Plot), length(Formated)), 
         plot.objects=c("observations", "ML", "ML.SD", "ML.quantiles"))
    
    output$resultsInfo <- renderText({
      
    rtxt <- paste("\n", babel["NumeroSerie" , ilanguage], nrow(x$synthesis), "\n", sep = "")
      
      for (i in 1:nrow(x$synthesis)) {
         nmser <- rownames(x$synthesis)[i]
         tirets <- paste0(rep("-", nchar(nmser)), sep="", collapse = "")
         pp <- ""
         if (i != 1)
           pp <- "\n"
         rtxt <- paste0(rtxt, pp, tirets, "\n", sep = "", collapse = "")
         rtxt <- paste0(rtxt, nmser, "\n", collapse = "")
         rtxt <- paste0(rtxt, tirets, "\n", sep = "", collapse = "")
         mnponte <- x$synthesis[i, "with_obs_Mean"]
         # sdponte <- x$synthesis[[i]]$estimates["SD.ML.with.obs"]
        rtxt <- paste0(rtxt, babel["Total" , ilanguage], format(mnponte, digits = floor(log10(mnponte) +
                                                                                          5)), "\n", sep = "", collapse = "")
         mnponte1 <- x$synthesis[i, "with_obs_Low_ML"]
         mnponte2 <- x$synthesis[i, "with_obs_High_ML"]
        if ((mnponte1 == 0) | (mnponte2 == 0)) {
          rtxt <- paste0(rtxt, babel["IC95ND" , ilanguage])
        }
        else {
          rtxt <- paste0(rtxt, babel["IC95" , ilanguage], format(mnponte1,
                                                                 digits = floor(log10(mnponte1) + 5)), " - ",
                         format(mnponte2, digits = floor(log10(mnponte2) +
                                                           5)), sep="", collapse = "")
        }
      }
      rtxt
    })
    
    }
    
    }

    })
  

  output$distPlot <- renderPlot({

    outplotfn()

  })

  
  
})