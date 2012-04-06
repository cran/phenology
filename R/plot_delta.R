#' plot_delta plots the likelihood delta for fixed Phi value.
#' @title Plot a likelihood lineplot obtained after map_phenology.
#' @author Marc Girondot
#' @return Return None
#' @param map A map generated with map_phenology
#' @param Phi Phi value or NULL
#' @param pdf TRUE or FALSE, indicates if a pdf file is generated.
#' @param pdfname Name of pdf file
#' @param help If TRUE, an help is displayed
#' @description This function plots a likelihood lineplot obtained after map_phenology.
#' @export

plot_delta <-
function(map=NULL, Phi=NULL, pdf=FALSE, pdfname="Map.pdf", help=FALSE) {

if(help || is.null(map)) {
	cat("This function plots a likelihood lineplot obtained after map_phenology.\n")
	cat("The syntaxe is:\n")
	cat("plot_delta(map=mapx, pdf=TRUE, pdfname='NameMap.pdf')\n")

} else {

input<-map$input

if (is.null(Phi)) {
Dv=as.vector(input)
pos=which.min(Dv)
j0=floor(pos/nrow(input))+1
i0=pos%%nrow(input)
} else {
i0<-which(map$Phi>=Phi)[1]
}

effetDelta<-map$input[i0,]

plot(map$Delta, effetDelta, type="l", xlab="Delta", ylab="-Ln L", bty="n", main=paste(map$Data, " - Phi=", map$Phi[i0], sep=""))

if (pdf) {

pdf(pdfname)

plot(map$Delta, effetDelta, type="l", xlab="Delta", ylab="-Ln L", bty="n", main=paste(map$Data, " - Phi=", map$Phi[i0], sep=""))

dev.off()

}
}
}
