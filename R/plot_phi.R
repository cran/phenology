#' plot_phi plots the best likelihood for fixed Phi value.
#' @title Plot the best likelihood for fixed Phi value.
#' @author Marc Girondot
#' @return Return None
#' @param map A map generated with map_phenology
#' @param pdf TRUE or FALSE, indicates if a pdf file is generated.
#' @param pdfname Name of pdf file
#' @param help If TRUE, an help is displayed
#' @description The function "plot_phi" plots the best likelihood for fixed Phi value.
#' @export


plot_phi <-
function(map=NULL, pdf=FALSE, pdfname="Map.pdf", help=FALSE) {

if(help || is.null(map)) {
	cat("This function plots a likelihood lineplot obtained after map_phenology.\n")
	cat("The syntaxe is:\n")
	cat("plot_phi(map=mapx, pdf=TRUE, pdfname='NameMap.pdf')\n")

} else {

effetphi<-NULL
for(i in 1:length(map$Phi)) effetphi<-c(effetphi, min(map$input[i,], na.rm=TRUE))
plot(map$Phi, effetphi, type="l", xlab="Phi", ylab="-Ln L", bty="n", main=map$Data)

if (pdf) {

pdf(pdfname)

plot(map$Phi, effetphi, type="l", xlab="Phi", ylab="-Ln L", bty="n", main=map$Data)

dev.off()

}
}
}
