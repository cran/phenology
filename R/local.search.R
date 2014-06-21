#' local.search() returns path of file serached in local disk based on its file name
#' @title Return path of file searched for in local disk based on its file name
#' @author Marc Girondot
#' @return A vector with paths
#' @param pattern The name of file to be searched for. Can use wildcards *
#' @param directory The path of directory to be explored in for Windows
#' @param folder The path of folder to be explored in for Unix based systems
#' @description Return path of file searched for in local disk based on its file name.
#' It has been tested only with Windows XP and MacOSX.
#' @examples
#' \dontrun{
#' RnwFiles <- local.search("*.Rnw")
#' }
#' @export


local.search <- function(pattern, directory="", folder="$HOME") {
  if (.Platform$OS.type=="unix") {
    command <- paste0("find ", folder," -type f -name '", pattern, "'")
    dest <- system(command, intern=TRUE, ignore.stdout=FALSE, ignore.stderr=TRUE)
    } else {
    dest <- shell(paste0("dir ", directory, pattern, " /b/s "), 
                          intern=TRUE, ignore.stdout=FALSE, ignore.stderr=TRUE)
  }
  if (length(dest) == 0) {return(NULL)} else {return(dest)}
}