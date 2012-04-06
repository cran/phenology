#' plot_map plots a likelihood map with Delta and Phi varying.
#' @title Plot a likelihood map with Delta and Phi varying.
#' @author Marc Girondot
#' @return Return None
#' @param map A map generated with map_phenology.
#' @param pdf TRUE or FALSE, indicates if a pdf file is generated.
#' @param pdfname Name of file if pdf=TRUE
#' @param col Colors could be heat.colors(128) or rainbow(64) or col=gray(c(seq(0, 1, length.out=128)))
#' @param help If TRUE, an help is displayed
#' @description This function plots a likelihood map obtained after map_phenology.
#' @export



plot_map <-
function(map=NULL, pdf=FALSE, pdfname="Map.pdf", col=heat.colors(128), help=FALSE) {

if(help || is.null(map)) {
	cat("This function plots a likelihood map obtained after map_phenology.\n")
	cat("The syntaxe is:\n")
	cat("plot_map(map=mapx, pdf=TRUE, pdfname='NameMap.pdf')\n")
	cat("The color can be changed using the optional parameter col, for example\n")
	cat("col=heat.colors(128) or rainbow(64) or col=gray(c(seq(0, 1, length.out=128)))\n")

} else {
	
	
require(grDevices)
require(graphics)

x <- map$Phi
y <- map$Delta

input<-map$input

if (is.element('fields', installed.packages()[,1]) == FALSE) {install.packages('fields') }
library("fields")


image.plot(x, y, input, col=col, axes=TRUE, xlab="Phi", ylab="Delta", nlevel = length(col))
image.plot(x, y, input, zlim=c(min(input, na.rm=TRUE), max(input, na.rm=TRUE)), col=col, axes=TRUE, xlab="Phi", ylab="Delta", nlevel = length(col))
  
if (pdf) {

pdf(pdfname)

image.plot(x, y, input, col=col, axes=TRUE, xlab="Phi", ylab="Delta", nlevel = length(col))
image.plot(x, y, input, zlim=c(min(input, na.rm=TRUE), max(input, na.rm=TRUE)), col=col, axes=TRUE, xlab="Phi", ylab="Delta", nlevel = length(col))
  
dev.off()

}
	
}
}
