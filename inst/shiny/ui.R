# library(shiny); library("phenology"); runApp("/Users/marc/Documents/Espace_de_travail_R/shiny/Phenology")
# library(shiny); runApp("http://max3.ese.u-psud.fr:3838/phenology/")


library(shiny)

load(file="babel.Rdata")

radioButtons_withHTML <- function (inputId, label, choices, selected = NULL, inline = FALSE, 
                                   width = NULL) 
{
  generateOptions_withHTML <- function (inputId, choices, selected, inline, type = "checkbox") 
  {
    options <- mapply(choices, names(choices), FUN = function(value, 
                                                              name) {
      inputTag <- tags$input(type = type, name = inputId, value = value)
      if (value %in% selected) 
        inputTag$attribs$checked <- "checked"
      if (inline) {
        tags$label(class = paste0(type, "-inline"), inputTag, 
                   tags$span(HTML(name)))
      }
      else {
        tags$div(class = type, tags$label(inputTag, tags$span(HTML(name))))
      }
    }, SIMPLIFY = FALSE, USE.NAMES = FALSE)
    div(class = "shiny-options-group", options)
  }
  choices <- shiny:::choicesWithNames(choices)
  selected <- if (is.null(selected)) 
    choices[[1]]
  else {
    shiny:::validateSelected(selected, choices, inputId)
  }
  if (length(selected) > 1) 
    stop("The 'selected' argument must be of length 1")
  options <- generateOptions_withHTML(inputId, choices, selected, inline, 
                                      type = "radio")
  divClass <- "form-group shiny-input-radiogroup shiny-input-container"
  if (inline) 
    divClass <- paste(divClass, "shiny-input-container-inline")
  tags$div(id = inputId, style = if (!is.null(width)) 
    paste0("width: ", validateCssUnit(width), ";"), class = divClass, 
    shiny:::controlLabel(inputId, label), options)
}



shinyUI(fluidPage(
  
  wellPanel(
    p(radioButtons_withHTML('language', label=NULL, choices = list('<img src="gb.png">'="gb"
                                                                            , '<img src="pt.png">'="pt"
                                                                            , '<img src="fr.png">'="fr"
                                                                   , '<img src="sp.png">'="sp"
                                                                   )
                            , selected="gb"
                            , inline = TRUE) 
    #   radioButtons(inputId="language"
    #                , label=NULL
    #                , choices=list(English="gb"
    #                               , português="pt"
    #                               , français="fr"
    #                ) 
    #                , selected="gb"
    #                , inline = TRUE
    # )
    , align="right")
  )
  
  # Application title
  
  , wellPanel(
    uiOutput("Titre_UI")
    , p(img(src="logo.png", height=60, width=626), align="center")
    , uiOutput("Objet_UI")
    , uiOutput("WebVersion_UI")
    , uiOutput("output$Developpeur_UI")
    , tags$hr()
    , uiOutput("Publication_UI")
    , tags$hr()
    , uiOutput("Liste_UI")
    
    , p(a("Girard A., Godgenger M.-C., Gibudi A., Fretey J., Billes A., Roumet D., Bal G., Bréheret N., Bitsindou A., Leeuwe H.V., Verhage B., Ricois S., Bayé J.-P., Carvalho J., Lima H., Neto E., Angoni H., Ayissi I., Bebeya C., Folack J., Nguegim J.R. and Girondot M. (2016) Marine turtles nesting activity assessment and trend along the Central African Atlantic coast for the period of 1999-2008. International Journal of Marine Science and Ocean Technology, 3(3), 21-32."
          , href="https://www.researchgate.net/publication/304572110_Marine_turtles_nesting_activity_assessment_and_trend_along_the_Central_African_Atlantic_coast_for_the_period_of_1999-2008"
          , target="_blank"))
    , p(a("Delcroix E., Bédel S., Santelli G. and Girondot M. (2013) Monitoring design for quantification of marine turtle nesting with limited human effort: a test case in the Guadeloupe Archipelago. Oryx, 48(1), 95-105."
          , href="https://www.researchgate.net/publication/259527763_Monitoring_design_for_quantification_of_marine_turtle_nesting_with_limited_effort_A_test_case_in_the_Guadeloupe_archipelago"
          , target="_blank"))    
    , p(a("Girondot M., Bédel S., Delmoitiez L., Russo M., Chevalier J., Guéry L., Hassine S.B., Féon H. and Jribi I. (2015) Spatio-temporal distribution of ", em("Manta birostris"), " in French Guiana waters. Journal of the Marine Biological Association of the United Kingdom, 95(1), 153-160."
        , href="https://www.researchgate.net/publication/270898326_Spatio-temporal_distribution_of_Manta_birostris_in_French_Guiana_waters"
        , target="_blank"))
    , p(a("Girondot M. and Rizzo A. (2015) Bayesian framework to integrate traditional ecological knowledge into ecological modeling: A case study. Journal of Ethnobiology, 35(2), 337-353."
          , href="https://www.researchgate.net/publication/279848741_Bayesian_Framework_to_Integrate_Traditional_Ecological_Knowledge_into_Ecological_Modeling_A_Case_Study"
          , target="_blank"))
    , p("Weber S.B., Weber N., Ellick J., Avery A., Frauenstein R., Godley B.J., Sim J., Williams N. and Broderick A.C. (2014) Recovery of the South Atlantic’s largest green turtle nesting population. Biodiversity and Conservation, 23(12), 3005-3018.")
    , p("Whiting A.U., Chaloupka M., Pilcher N., Basintal P. and Limpus C.J. (2014) Comparison and review of models describing sea turtle nesting abundance. Marine Ecology Progress Series, 508, 233-246.")
    , p("Stringell T.B., Clerveaux W.V., Godley B.J., Phillips Q., Ranger S., Richardson P.B., Sanghera A. and Broderick A.C. (2015) Protecting the breeders: research informs legislative change in a marine turtle fishery. Biodiversity and Conservation, 24(7), 1775-1796.")
    , p("Bellini C., Santos A.J.B., Grossman A., Marcovaldi M.A. and Barata P.C.R. (2013) Green turtle (", em("Chelonia mydas"), ") nesting on Atol das Rocas, north-eastern Brazil, 1990–2008. Journal of the Marine Biological Association of the United Kingdom, 93(04), 1117-1132.")
    
    )
  , wellPanel(    
    a(img(src="logoSWOT.png", height=72, width=215)
      , href="http://www.seaturtlestatus.org"
      , target="_blank")
    , uiOutput("RecommandeSWOT_UI")
  )
  , wellPanel(
    a(img(src="logoRASTOMA.png", height=79, width=532)
      , href="http://www.rastoma.org"
      , target="_blank")
    , uiOutput("RecommandeRASTOMA_UI")
  )
  
  # Sidebar with a slider input for the number of bins
  , sidebarLayout(
    sidebarPanel(
      wellPanel(
        uiOutput("GSYalimapo_UI0")
        , uiOutput("GSYalimapo_UI")
      )
      , wellPanel(
        uiOutput("ChargeData_UI")
        , fileInput("file1", label=NULL, accept=c('text/csv', 
                                                  'text/comma-separated-values',
                                                  'text/plain', '.csv'))
        , uiOutput("TitreCol_UI")
        , actionButton("goButton", "Go!")
        , uiOutput("NoteFormat_UI")
        )
      ,  wellPanel(
        uiOutput("Options_UI")
        , uiOutput("Min0_UI")
        , uiOutput("Sym_UI")
        , uiOutput("Report0_UI")
        , uiOutput("FormatDate_UI")
        , uiOutput("NoteDate_UI")
        , uiOutput("MoisRef_UI")
        , uiOutput("NoteMoisRef_UI")
        , uiOutput("SerieGraphique_UI")
      )
      )
    ,
    
    # Show a plot of the generated distribution
    mainPanel(
      plotOutput("distPlot")
      , uiOutput("NoteGraph_UI")
      , uiOutput("Synthese1_UI")
      , verbatimTextOutput(outputId="resultsInfo")
#      , uiOutput("Synthese2_UI")
#      , verbatimTextOutput(outputId="resultsRaw")
      , tableOutput("table")
      , uiOutput("NoteStockage_UI")
      )
    )
      ))