#' Parameter_Global_Year transforms a set of parameters from Year to global effect, or reverse
#' @title Transform a set of parameters from Year to global effect, or reverse
#' @author Marc Girondot
#' @return Return a set of modified parameters
#' @param parameters Set of current parameters
#' @param parname Name of parameter to transform
#' @param series Set of series (see description)
#' @param sep_year Character used to separate the year
#' @param perYear TRUE if year-specific values must be setup
#' @description This function is used to transform a set of parameters 
#' that uses Peak, LengthB, LengthE, B, E, or Length to the same parameter 
#' with Year effect, or reverse.\cr
#' The parameter series can be or a result from add_phenology() or from 
#' fit_phenology() or simply a vector of names.
#' @examples 
#' \dontrun{
#' Parameter_Global_Year(parameters=c("Peak_Beach1-2018"=151, "Peak_Beach1-2019"=161), 
#'                      parname="Peak", perYear=FALSE)
#' Parameter_Global_Year(parameters=c("Peak"=151), 
#'                       series = c("beach1", "beach2"), 
#'                      parname="Peak", perYear=TRUE)
#' }
#' @export


Parameter_Global_Year <-
  function(parameters=stop("A set of parameters must be indicated"), 
           parname=c("Peak", "LengthB", "LengthE", "B", "E", "Length"), 
           series=NULL, 
           sep_year="-", 
           perYear=TRUE) {
    
    parname <- match.arg(parname, c("Peak", "LengthB", "LengthE", "B", "E", "Length"), 
                         several.ok = FALSE)
    
    if (!is.null(class(series))) {
      if (class(series)=="phenologydata") series <- names(series)
      if (class(series)=="phenology") series <- names(series$data)
    }
    px <- which(substr(names(parameters), 1, nchar(parname)) == parname)
    if (perYear) {
      # J'en ai un seul; je dois créer
      if (is.null(series)) stop("A list of series must be provided")
      sseries <- gsub(paste0(sep_year, "([^", sep_year, "]*)$"), "_\\1", series)
      
      cn <- strsplit(sseries, "_")
      rg <- unique(sapply(cn, FUN = function(x) x[[length(x)]]))
      
      ppname <- paste0(parname, "_", rg)
      pp <- rep(parameters[px], length(ppname))
      names(pp) <- ppname
      parameters <- c(parameters[-px], pp)
    } else {
      
      pp <- mean(parameters[px])
      names(pp) <- parname
      parameters <- c(parameters[-px],  pp)
      
    }
    
    
    return(parameters)
  }
