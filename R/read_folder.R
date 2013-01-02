#' read_folder opens all files from a folder
#' @title Open all file present in a folder and create a list with these files
#' @author Marc Girondot
#' @return Return a list of the data in the files of the folder (directory for windows users)
#' @param folder Where to search for files
#' @param read Function used to read file. Ex: read.delim
#' @param ... Parameters send to the read function
#' @description To create a new list, the syntaxe is \cr
#' datalist<-read_folder(folder=".", read=read.delim, header=FALSE)\cr
#' @examples 
#' \dontrun{
#' library(phenology)
#' # Read all the files from a folder/directory
#' Gratiot<-read_folder(".")
#' }
#' @export


read_folder <- function(folder=".", read=read.delim, ...) {

	lf <- list.files(folder)
	
	ladd <- list()
	
	if (length(lf)!=0) {
	
	previous<-getwd()
	
	setwd(folder)
	
	for (i in 1:length(lf)) {
	
		linter <- read(lf[i], ...)
		ladd <- c(ladd, list(linter))
	
	}
	
	names(ladd) <- lf
	
	setwd(previous)
	
	return(ladd)
	
	} else {
		print("No file in folder/directory or folder/directory does not exist")
		return(NULL)
	}

}
