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
#' @examples
#' library("phenology")
#' # Read a file with data
#' # Gratiot<-read.delim("http://max2.ese.u-psud.fr/epc/conservation/BI/Complete.txt", , header=FALSE)
#' data(Gratiot)
#' # Generate a formatted list nammed data_Gratiot 
#' data_Gratiot<-add_format(origin=NULL, add=Gratiot, name="Complete", reference=as.Date("2001-01-01"), format="%d/%m/%Y")
#' # Generate initial points for the optimisation
#' parg<-par_init(data_Gratiot, parametersfixed=NULL)
#' # Run the optimisation
#' \dontrun{
#' result_Gratiot<-fit_phenology(data=data_Gratiot, parametersfit=parg, parametersfixed=NULL, trace=1)
#' }
#' data(result_Gratiot)
#' # Extract the fitted parameters
#' parg1<-extract_result(result_Gratiot)
#' # Add constant Alpha and Tau values [day d amplitude=(Alpha+Nd*Beta)^Tau with Nd being the number of counts for day d]
#' pfixed<-c(parg1, Alpha=0, Tau=1)
#' pfixed<-pfixed[-which(names(pfixed)=="Theta")]
#' # The only fitted parameter will be Beta
#' parg2<-c(Beta=0.5, parg1["Theta"])
#' # Generate a likelihood map [default Phi=seq(from=0.1, to=20, length.out=100) but it is very long]
#' # Take care, it takes 20 hours ! The data map_Gratiot has the result
#' \dontrun{
#' map_Gratiot<-map_phenology(data=data_Gratiot, Phi=seq(from=0.1, to=20, length.out=100), parametersfit=parg2, parametersfixed=pfixed)
#' }
#' data(map_Gratiot)
#' # Plot the map
#' plot_map(map=map_Gratiot, pdf=FALSE, col=heat.colors(128))
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
