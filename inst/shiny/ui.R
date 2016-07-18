# library(shiny); library("phenology"); runApp("/Users/marc/Documents/Espace_de_travail_R/shiny/Phenology")
# library(shiny); runApp("http://max3.ese.u-psud.fr:3838/phenology/")


library(shiny)

# Define UI for application that draws a histogram
shinyUI(fluidPage(
  
  # Application title
  
  wellPanel(
    h1("Phenology analysis", align = "center")
    , p(img(src="logo.png", height=126, width=626), align="center")
    , p(strong("The phenology R package is a set of tools to easily estimate 
               density of animal present in a site when animals are present 
               in a site during only part of the year."))
    , p("This web server version is a simplified version of the complete tools available ", 
        a("here."
          , href="https://cran.r-project.org/package=phenology"
          , target="_blank"))
    , p("Phenology package is developped by "
        , a("Marc Girondot"
            , href="http://max2.ese.u-psud.fr/epc/conservation/Girondot/Publications/Marc.html"
            , target="_blank"))
    , tags$hr()
    , p(" The method has been published in ",
        a(
          img(src="PDF.png", height=35, width=35), 
          "Girondot, M. 2010. Estimating density of animals during migratory waves: 
          application to marine turtles at nesting site. Endangered Species Research, 
          12, 85-105.", 
          href="http://www.int-res.com/articles/esr2010/12/n012p095.pdf",
          target="_blank"))
  )
  
  # Sidebar with a slider input for the number of bins
  , sidebarLayout(
    sidebarPanel(
      wellPanel(
        h4("Load dataset with counts")
        , fileInput("file1", label=NULL, accept=c('text/csv', 
                                          'text/comma-separated-values',
                                          'text/plain', '.csv'))
        , checkboxInput("Header", "Have the data column headers ?", value = FALSE, width = NULL)
        , helpText("Note: The input file must be a .csv or .txt. 
                   The first column must be the dates and the second be the counts. 
                   An optional third column can be used to indicate the origin 
                   of data (for example using name of beaches) to aggregate several datasets 
                   in a single analysis.")
        )
      ,  wellPanel(h4("Options to evaluate phenology")
                   , radioButtons("Min", "Should the Min parameter be fixed to 0 ?", list(yes=1, no=2), selected=2, inline = TRUE)
                   , radioButtons("Length", "Should the season be symmetrical around peak ?", list(yes=1, no=2), selected=2, inline = TRUE)
                   , radioButtons("Zero", "Have the 0 counts been reported ?", list(yes=1, no=2), selected=1, inline = TRUE)
      )
        )
    ,
    
    # Show a plot of the generated distribution
    mainPanel(
      plotOutput("distPlot")
      , h4("Synthesis for nest counts")
      , verbatimTextOutput(outputId="resultsInfo")
      , h4("Synthesis for phenology model")
      , verbatimTextOutput(outputId="resultsRaw")
    )
    )
))