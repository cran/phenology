#' fixed.parameters0 generates a set of fixed parameters for series with only 0 counts
#' @title Generate a set of fixed parameters for series with only 0 counts
#' @author Marc Girondot
#' @return Return a set of parameters
#' @param series Set of series generated with add_phenology()
#' @description This function generates a set of fixed parameters for series with only 0 counts.\cr
#' The parameter series must be a result from add_phenology().
#' @examples 
#' \dontrun{
#' refdate <- as.Date("2001-01-01")
#' data_Gratiot <- add_phenology(Gratiot, name="Complete", 
#' 	reference=refdate, format="%d/%m/%Y")
#' pfixed <- fixed.parameters0(data_Gratiot)
#' }
#' @export


fixed.parameters0 <-
  function(series=stop("A result from add_phenology() must be provided.")) {
    
    pfixed <- NULL
    sumSeries <- sapply(series, FUN = function(x) sum(x$nombre))
    sumSeries0 <- names(sumSeries[sumSeries == 0])
    if (!identical(sumSeries0, character(0))) {
      p0 <- rep(0, length(sumSeries0))
      names(p0) <- paste0("Max_", sumSeries0)
      pfixed <- c(pfixed, p0)
    }
    
    
    
    return(pfixed)
  }
