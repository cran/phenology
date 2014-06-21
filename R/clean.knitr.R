#' clean.knitr deletes temporary files created during knitr compile
#' @title Delete temporary files created during knitr compile
#' @author Marc Girondot
#' @return Nothing
#' @description Delete temporary files created during knitr compile in working directory
#' @examples
#' \dontrun{
#' clean.knitr()
#' }
#' @export


clean.knitr <- function() system(paste0("cd '", getwd(), "';rm *.gz;rm *.toc;rm *.tex;rm *.log"))