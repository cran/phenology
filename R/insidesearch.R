#' insidesearch Search a string within files of a folder
#' @title Search a string within files of a folder
#' #' @author Marc Girondot
#' @return Return a vector with filenames
#' @param path Path of the folder to search
#' @param pattern Pattern for file names to searcg
#' @param search String to search in files
#' @param ... Options for readLines(), example warn = FALSE
#' @description Search for a string inside the files of a folder \cr
#' @examples
#' \dontrun{
#' library(phenology)
#' # Read a file with data
#' }
#' @export


insidesearch <- function(path=".", pattern="*.R", ..., 
                         search=stop("A text to be searched for is necessary")) {
  ls <- list.files(path, pattern)
  if (identical(ls, character(0))) {
    warning("No file match this pattern at this path")
    return(invisible(NULL))
  }
  
  filesok <- NULL
  
  for (f in ls) {
    fc <- readLines(con = f, ...)
    x <- gsub(search, "", fc)
    if (any(x != fc)) {
      cat("file ", f, "\n")
      cat("lines ", c(1:length(fc))[x != fc], "\n")
      filesok <- c(filesok, f)
    }
  }
  return(invisible(filesok))
}