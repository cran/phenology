#' compare_AIC compares the AIC of several outputs obtained with the same data.
#' @title Compares the AIC of several outputs
#' @author Marc Girondot
#' @return A list with DeltaAIC and Akaike weight for the models.
#' @param result A list with the results to be compared
#' @param help If TRUE, an help is displayed
#' @description This function is used to compares the AIC of several outputs obtained with the same data but with different set of parameters.
#' @examples
#' # Read a file with data
#' library("phenology")
#' #' Gratiot<-read.delim("http://max2.ese.u-psud.fr/epc/conservation/BI/Complete.txt", , header=FALSE)
#' data(Gratiot)
#' # Generate a formated list nammed data_Gratiot 
#' data_Gratiot<-add_format(origin=NULL, add=Gratiot, name="Complete", reference=as.Date("2001-01-01"), format="%d/%m/%Y")
#' # Fix parameter FLat to 0
#' pfixed=c(Flat=0)
#' # Generate initial points for the optimisation
#' parg<-par_init(data_Gratiot, parametersfixed=pfixed)
#' # Fit is done
#' # result_Gratiot_Flat<-fit_phenology(data=data_Gratiot, parametersfit=parg2, parametersfixed=pfixed, trace=1)
#' data(result_Gratiot_Flat)
#' # Generate initial points for the optimisation
#' parg<-par_init(data_Gratiot, parametersfixed=NULL)
#' # Run the optimisation
#' # result_Gratiot<-fit_phenology(data=data_Gratiot, parametersfit=parg, parametersfixed=NULL, trace=1)
#' data(result_Gratiot)
#' # Compare both models
#' outputAIC<-compare_AIC(list(full=result_Gratiot, Flat=result_Gratiot_Flat))
#' @export


compare_AIC <-
function(result=NULL, help=FALSE) {
if(help) {
	cat("This function is used to compares the AIC of several outputs obtained\n")
	cat("with the same data but with different set of parameters.\n")
	cat("The syntax is outputAIC<-compare_AIC(list(test1=result_test1, test2=result_test2))\n")

} else {
if (!is.null(result)) {
	if (!is.list(result) || (is.null(names(result))) || (any(names(result)==""))) {
		print("The results must be included within a list with names; see example.")
	} else {
	out<-NULL
	l<-length(result)
		if (l<2) {
			print("A least two results must be provided.")
		} else {
			aic<-NULL
			name<-names(result)
			for (i in 1:l) {
				if (inherits(result[[i]], "NestsResult")) {
					aic<-c(aic, result[[i]]$AIC)
				} else {
					aic<-c(aic, 2*result[[i]]$value+2*(length(result[[i]]$par)))
				}
			}
			
			bestaic<-min(aic)
			ser<-which(aic==bestaic)
			deltaaic<-aic-bestaic
			aw<-exp(-0.5*deltaaic)
			saw=sum(aw)
			aw<-aw/saw
			
			out<-data.frame(cbind(AIC=aic, DeltaAIC=deltaaic, Akaike_weigth=aw), row.names=name)
			print(paste("The lowest AIC (",bestaic ,") is for series ", name[ser], " with Akaike weight=", aw[ser], sep=""))
			
			return(out)
		}
	}
}
}
}
