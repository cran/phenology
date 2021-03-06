#' plot_phi plots the best likelihood for fixed Phi value.
#' @title Plot the best likelihood for fixed Phi value.
#' @author Marc Girondot
#' @return Return None
#' @param map A map generated with map_phenology
#' @param ... Parameters for plot
#' @description The function "plot_phi" plots the best likelihood for each Phi value.
#' @family Phenology model
#' @examples
#' library("phenology")
#' # Read a file with data
#' data(Gratiot)
#' # Generate a formatted list nammed data_Gratiot 
#' data_Gratiot<-add_phenology(Gratiot, name="Complete", 
#' 		reference=as.Date("2001-01-01"), format="%d/%m/%Y")
#' # Generate initial points for the optimisation
#' parg<-par_init(data_Gratiot, fixed.parameters=NULL)
#' # Run the optimisation
#' \dontrun{
#' result_Gratiot<-fit_phenology(data=data_Gratiot, 
#' 		fitted.parameters=parg, fixed.parameters=NULL)
#' }
#' data(result_Gratiot)
#' # Extract the fitted parameters
#' parg1<-extract_result(result_Gratiot)
#' # Add constant Alpha and Tau values 
#' # [day d amplitude=(Alpha+Nd*Beta)^Tau with Nd being the number of counts for day d]
#' pfixed<-c(parg1, Alpha=0, Tau=1)
#' pfixed<-pfixed[-which(names(pfixed)=="Theta")]
#' # The only fitted parameter will be Beta
#' parg2<-c(Beta=0.5, parg1["Theta"])
#' # Generate a likelihood map [default Phi=seq(from=0.1, to=20, length.out=100) but it is very long]
#' # Take care, it takes 20 hours ! The data map_Gratiot has the result
#' \dontrun{
#' map_Gratiot<-map_phenology(data=data_Gratiot, 
#' 		Phi=seq(from=0.1, to=20, length.out=100), 
#' 		fitted.parameters=parg2, fixed.parameters=pfixed)
#' }
#' data(map_Gratiot)
#' # Plot the min(-Ln L) for Phi varying at any delta value
#' plot_phi(map=map_Gratiot)
#' @export


plot_phi <-
function(map=NULL, ...) {

effetphi<-NULL
for(i in 1:length(map$Phi)) effetphi<-c(effetphi, min(map$input[i,], na.rm=TRUE))
do.call(plot, modifyList(list(x=map$Phi, y=effetphi, type="l", xlab="Phi", ylab="-Ln L", bty="n", main=""), list(...)))

}
