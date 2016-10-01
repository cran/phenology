library(shiny)
package.phenology <- require('phenology')

load(file="babel.Rdata")

# Define server logic required to draw a histogram
shinyServer(function(input, output) {
  
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
          , href="http://max2.ese.u-psud.fr/epc/conservation/Girondot/Publications/Marc.html"
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
  
  output$GSYalimapo_UI0 <- renderUI({
    h4(babel["ChargeExemple" , input$language])
  })
  output$GSYalimapo_UI <- renderUI({
    actionButton("GSYalimapo", babel["ChargeData" , input$language])
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
  

  output$distPlot <- renderPlot({
    
    ilanguage <- input$language
    igo <- input$goButton
    iGS <- input$GSYalimapo
    
    
    isolate({
    ifile1 <- input$file1
    iHeader <- input$Header
    iMonthRef <- input$MonthRef
    iFormatDate <- input$FormatDate
    iMin <- input$Min
    iLength <- input$Length
    iZero <- input$Zero
    
    
    
    if (!package.phenology) {
      par(mar=c(0, 0, 0, 0))
      plot(x=c(0, 1), y=c(0, 1), axes=FALSE,
           xaxt="n", yaxt="n", main="",
           xlab = "", ylab = "",
           xaxs="i", yaxs="i", type="n")
      text(x = 0.5, y=0.5, labels = babel["Install1" , ilanguage]
           , col="red", cex = 2)
      text(x = 0.5, y=0.3, labels = babel["Install2" , ilanguage]
           , col="green", cex = 2)
    } else {
      
      iGS2 <- ifelse(length(iGS)==0, TRUE, iGS == 0)
      
      if (identical(is.null(ifile1) & iGS2, TRUE)) {
    
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
        
        y2002 <- !iGS2
        go <- (igo != 0)
        
        if (!is.null(ifile1)) {
          inFile <- ifile1
          dest <- inFile$datapath
          
          # dest <- "/Users/marc/Documents/Espace_de_travail_R/package_phenology/LUTH\ 2002/Yalimapo_2002_GS.txt"
          
          table <- read.delim(dest, header = iHeader, stringsAsFactors = FALSE)
        } else {
          table <-   structure(list(Date = c("01/03/02", "02/03/02", "03/03/02", "04/03/02", 
                                           "05/03/02", "06/03/02", "07/03/02", "08/03/02", "09/03/02", "10/03/02", 
                                           "11/03/02", "12/03/02", "13/03/02", "14/03/02", "15/03/02", "16/03/02", 
                                           "17/03/02", "18/03/02", "19/03/02", "20/03/02", "21/03/02", "22/03/02", 
                                           "23/03/02", "24/03/02", "25/03/02", "26/03/02", "27/03/02", "28/03/02", 
                                           "29/03/02", "30/03/02", "31/03/02", "01/04/02", "02/04/02", "03/04/02", 
                                           "04/04/02", "05/04/02", "06/04/02", "07/04/02", "08/04/02", "09/04/02", 
                                           "10/04/02", "11/04/02", "12/04/02", "13/04/02", "14/04/02", "15/04/02", 
                                           "16/04/02", "17/04/02", "18/04/02", "19/04/02", "20/04/02", "21/04/02", 
                                           "22/04/02", "23/04/02", "24/04/02", "25/04/02", "26/04/02", "27/04/02", 
                                           "28/04/02", "29/04/02", "30/04/02", "01/05/02", "02/05/02", "03/05/02", 
                                           "04/05/02", "05/05/02", "06/05/02", "07/05/02", "08/05/02", "09/05/02", 
                                           "10/05/02", "11/05/02", "12/05/02", "13/05/02", "14/05/02", "15/05/02", 
                                           "16/05/02", "17/05/02", "18/05/02", "19/05/02", "20/05/02", "21/05/02", 
                                           "22/05/02", "23/05/02", "24/05/02", "25/05/02", "26/05/02", "27/05/02", 
                                           "28/05/02", "29/05/02", "30/05/02", "31/05/02", "01/06/02", "02/06/02", 
                                           "03/06/02", "04/06/02", "05/06/02", "06/06/02", "07/06/02", "08/06/02", 
                                           "09/06/02", "10/06/02", "11/06/02", "12/06/02", "13/06/02", "14/06/02", 
                                           "15/06/02", "16/06/02", "17/06/02", "18/06/02", "19/06/02", "20/06/02", 
                                           "21/06/02", "22/06/02", "23/06/02", "24/06/02", "25/06/02", "26/06/02", 
                                           "27/06/02", "28/06/02", "29/06/02", "30/06/02", "01/07/02", "02/07/02", 
                                           "03/07/02", "04/07/02", "05/07/02", "06/07/02", "07/07/02", "08/07/02", 
                                           "09/07/02", "10/07/02", "11/07/02", "12/07/02", "13/07/02", "14/07/02", 
                                           "15/07/02", "16/07/02", "17/07/02", "18/07/02", "19/07/02", "20/07/02", 
                                           "21/07/02", "22/07/02", "23/07/02", "24/07/02", "25/07/02", "26/07/02", 
                                           "27/07/02", "28/07/02", "29/07/02", "30/07/02", "31/07/02", "01/08/02", 
                                           "02/08/02", "03/08/02", "04/08/02", "05/08/02", "06/08/02", "07/08/02", 
                                           "08/08/02", "09/08/02", "10/08/02", "11/08/02", "12/08/02", "13/08/02", 
                                           "14/08/02", "15/08/02", "16/08/02", "17/08/02", "18/08/02", "19/08/02", 
                                           "20/08/02", "21/08/02", "22/08/02", "23/08/02", "24/08/02", "25/08/02", 
                                           "26/08/02", "27/08/02", "28/08/02", "29/08/02", "30/08/02", "31/08/02", 
                                           "01/09/02", "02/09/02", "03/09/02", "04/09/02", "05/09/02", "06/09/02", 
                                           "07/09/02", "08/09/02", "09/09/02", "10/09/02", "11/09/02", "12/09/02", 
                                           "13/09/02", "14/09/02", "15/09/02", "16/09/02", "17/09/02", "18/09/02", 
                                           "19/09/02", "20/09/02", "21/09/02", "22/09/02", "23/09/02", "24/09/02", 
                                           "25/09/02", "26/09/02", "27/09/02", "28/09/02", "29/09/02", "30/09/02", 
                                           "01/10/02", "02/10/02", "03/10/02", "04/10/02", "05/10/02", "06/10/02", 
                                           "07/10/02", "08/10/02", "09/10/02", "10/10/02", "11/10/02", "12/10/02", 
                                           "13/10/02", "14/10/02", "15/10/02", "16/10/02", "17/10/02", "18/10/02", 
                                           "19/10/02", "20/10/02", "21/10/02", "22/10/02", "23/10/02", "24/10/02", 
                                           "25/10/02", "26/10/02", "27/10/02", "28/10/02", "29/10/02", "30/10/02", 
                                           "31/10/02", "01/11/02", "02/11/02", "03/11/02", "04/11/02", "05/11/02", 
                                           "06/11/02", "07/11/02", "08/11/02", "09/11/02", "10/11/02", "11/11/02", 
                                           "12/11/02", "13/11/02", "14/11/02", "15/11/02", "16/11/02", "17/11/02", 
                                           "18/11/02", "19/11/02", "20/11/02", "21/11/02", "22/11/02"), 
                                    Number = c(3L, 2L, 1L, 1L, 6L, 2L, 2L, 2L, 1L, NA, 4L, NA, 1L, 
                                           NA, 1L, NA, 3L, 1L, 4L, 4L, 2L, 3L, 1L, 7L, 7L, 5L, 4L, 7L, 
                                           4L, 4L, NA, NA, NA, 10L, 21L, 9L, 9L, 17L, 11L, 16L, 18L, 
                                           32L, NA, NA, 30L, NA, 40L, 20L, 16L, 18L, 39L, 44L, 32L, 
                                           57L, 33L, 43L, 24L, 25L, 15L, 20L, 51L, 40L, 50L, 46L, 51L, 
                                           35L, 14L, 50L, 56L, 68L, 69L, 72L, 66L, NA, NA, 63L, 34L, 
                                           48L, 31L, NA, NA, NA, NA, 74L, 81L, 58L, NA, NA, 56L, 36L, 
                                           37L, 44L, 61L, NA, NA, 35L, NA, 62L, 44L, 69L, NA, NA, NA, 
                                           56L, 46L, 27L, 42L, 56L, NA, 69L, 59L, 55L, 35L, 47L, 38L, 
                                           NA, 30L, 28L, 38L, 24L, NA, NA, NA, 27L, 35L, 30L, 36L, 41L, 
                                           NA, NA, 27L, 21L, 27L, 15L, 22L, NA, 22L, 25L, 26L, 23L, 
                                           23L, NA, 14L, NA, NA, NA, NA, NA, 10L, 10L, 7L, 4L, 4L, 5L, 
                                           2L, NA, NA, 15L, 7L, 3L, NA, NA, NA, NA, 2L, NA, NA, NA, 
                                           NA, NA, 3L, 3L, 2L, 4L, 1L, 0L, NA, NA, 4L, 3L, 3L, 0L, NA, 
                                           NA, NA, 0L, 0L, 0L, 0L, 0L, NA, NA, 0L, 0L, 1L, 0L, 0L, NA, 
                                           NA, 1L, NA, 2L, NA, 0L, NA, NA, 2L, NA, 0L, NA, 0L, NA, NA, 
                                           0L, 1L, 0L, 1L, 1L, NA, NA, NA, 0L, 0L, 0L, 0L, 0L, NA, 3L, 
                                           NA, 0L, NA, 0L, NA, NA, 0L, NA, 1L, NA, 0L, NA, NA, 1L, NA, 
                                           NA, NA, NA, NA, NA, 0L, NA, 1L, NA, NA, NA, NA, NA, 1L, 1L, 
                                           0L, 0L, NA, NA, NA, NA, NA, NA, 1L)), .Names = c("Date", "Number"
                                           ), class = "data.frame", row.names = c(NA, -267L))
          
          inFile <- list(name = "Yalimapo 2002")
        }
        
        lby <- ncol(table)
        
        if (lby == 3) {
          nm <- levels(as.factor(table[, 3]))
          if (length(nm)>1) nm <- NULL
        } else {
          nm <- inFile$name
        }
        
        df2 <- table
        if (is.na(iMonthRef)) zmonthRef <- NULL else zmonthRef <- iMonthRef
        if (iFormatDate=="") zFormatDate <- NULL else zFormatDate <- iFormatDate
        
        output$table <- renderTable({
          df2
        })
        
        Formated <- add_phenology(add=df2, silent=TRUE, name=nm,
                                  format=zFormatDate,
                                  month_ref = zmonthRef)
        
        if (iMin == 2) {
          pfixed <- c(Flat=0)
        } else {
          pfixed <- c(Flat = 0, Min = 0)
        }
        Par <- par_init(data=Formated, parametersfixed=pfixed)
        if (iLength == 1) Par <- LBLE_to_L(Par)
        if (iMin == 3) {
          Par <- MinBMinE_to_Min(Par)
          Par <- toggle_Min_PMin(Par)
        }
        
        zLogic <- ifelse(iZero==1, TRUE, FALSE)
        
       result <- fit_phenology(data=Formated, parametersfit = Par,
                                 parametersfixed=pfixed,
                                 trace=0, silent=TRUE,
                                zero_counts = zLogic)
 
        
       x <- plot(result, series="all", moon=FALSE, progressbar = FALSE,
                 growlnotify = FALSE, Plot=FALSE)
       
       plot(result, series= as.numeric(input$Plot), moon=FALSE, progressbar = FALSE,
            growlnotify = FALSE)
       
       
       
       output$resultsInfo <- renderText({
         
         rtxt <- paste("\n", babel["NumeroSerie" , ilanguage], length(x), "\n", sep = "")
         for (i in 1:length(x)) {
           nmser <- names(x)[i]
           tirets <- paste0(rep("-", nchar(nmser)), sep="", collapse = "")
           pp <- ""
           if (i != 1) 
             pp <- "\n"
           rtxt <- paste0(rtxt, pp, tirets, "\n", sep = "", collapse = "")
           rtxt <- paste0(rtxt, nmser, "\n", collapse = "")
           rtxt <- paste0(rtxt, tirets, "\n", sep = "", collapse = "")
           mnponte <- x[[i]]$estimates[1]
           sdponte <- x[[i]]$estimates[2]
           if (sdponte >= 1e-04) {
             rtxt <- paste0(rtxt, babel["Total" , ilanguage], format(mnponte, digits = floor(log10(mnponte) + 
                                                                                               5)), " ; SD ", format(sdponte, digits = floor(log10(sdponte) + 
                                                                                                                                               5)), "\n", sep = "", collapse = "")
           }
           else {
             rtxt <- paste0(rtxt, babel["Total" , ilanguage], format(mnponte, digits = 3, 
                                                                     nsmall = 3), " ; SD NA\n", sep = "", collapse = "")
           }
           mnponte1 <- x[[i]]$estimates[5]
           mnponte2 <- x[[i]]$estimates[6]
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
       # output$resultsRaw <- renderPrint({print(result)})
        
        
        
      }
    }
    
    })

  })
})