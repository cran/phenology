#' plot_delta plots the likelihood delta for fixed Phi value.
#' @title Plot a likelihood lineplot obtained after map_phenology.
#' @author Marc Girondot
#' @return Return None
#' @param map A map generated with map_phenology
#' @param Phi Phi value or NULL
#' @param ... Parameters for plot
#' @description This function plots a likelihood lineplot obtained after map_phenology.
#' @family Phenology model
#' @examples
#' \dontrun{
#' library("phenology")
#' # Read a file with data
#' data(Gratiot)
#' # Generate a formatted list nammed data_Gratiot 
#' data_Gratiot<-add_phenology(Gratiot, name="Complete", 
#' 		reference=as.Date("2001-01-01"), format="%d/%m/%Y")
#' # Generate initial points for the optimisation
#' parg<-par_init(data_Gratiot, fixed.parameters=NULL)
#' # Run the optimisation
#' result_Gratiot<-fit_phenology(data=data_Gratiot, 
#' 		fitted.parameters=parg, fixed.parameters=NULL)
#' data(result_Gratiot)
#' # Extract the fitted parameters
#' parg1<-extract_result(result_Gratiot)
#' # Add constant Alpha and Tau values 
#' # [day d amplitude=(Alpha+Nd*Beta)^Tau with Nd being the number of counts for day d]
#' pfixed<-c(parg1, Alpha=0, Tau=1)
#' pfixed<-pfixed[-which(names(pfixed)=="Theta")]
#' # The only fitted parameter will be Beta
#' parg2<-c(Beta=0.5, parg1["Theta"])
#' # Generate a likelihood map 
#' # [default Phi=seq(from=0.1, to=20, length.out=100) but it is very long]
#' # Take care, it takes 20 hours ! The data map_Gratiot has the result
#' map_Gratiot<-map_phenology(data=data_Gratiot, 
#' 		Phi=seq(from=0.1, to=20, length.out=100), 
#' 		fitted.parameters=parg2, fixed.parameters=pfixed)
#' data(map_Gratiot)
#' # Plot the min(-Ln L) for Delta varying with Phi equal to the value for maximum likelihood
#' plot_delta(map=map_Gratiot)
#' # Plot the min(-Ln L) for Delta varying with Phi the nearest to 15
#' plot_delta(map=map_Gratiot, Phi=15)
#' }
#' @export

plot_delta <- function(map=NULL, Phi=NULL, ...) {
  
  input<-map$input
  
  if (is.null(Phi)) {
    Dv=as.vector(input)
    pos=which.min(Dv)
    j0=floor(pos/nrow(input))+1
    i0=pos%%nrow(input)
  } else {
    i0 <- which(map$Phi>=Phi)[1]
  }
  
  xregle <- c(0, 1+map$Phi[i0]/2)
  
  effetDelta <- map$input[i0,]
  
  do.call(plot, modifyList(list(x=map$Delta, y=effetDelta, type="l", 
                                xlab="Delta", ylab="-Ln L", xlim=xregle, 
                                bty="n", main=paste0("Phi=", map$Phi[i0])), 
                           list(...)))
  
}
