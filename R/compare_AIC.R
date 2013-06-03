#' compare_AIC compares the AIC of several outputs obtained with the same data.
#' @title Compares the AIC of several outputs
#' @author Marc Girondot
#' @return A list with DeltaAIC and Akaike weight for the models.
#' @param ... Successive results to be compared as lists
#' @description This function is used to compares the AIC of several outputs obtained with the same data but with different set of parameters.\cr
#' The parameters must be lists with $aic or $AIC or $value and $par elements.\cr
#' If several objects are within the same list, there AIC is sumed.\cr
#' For example, compare_AIC(g1=list(group), g2=list(separe1, separe2)) can be used to compare a single model onto two different sets of data against each set of data fited with its own set of parameters.
#' @examples
#' \dontrun{
#' # Read a file with data
#' library("phenology")
#' Gratiot<-read.delim("http://max2.ese.u-psud.fr/epc/conservation/BI/Complete.txt", header=FALSE)
#' data(Gratiot)
#' # Generate a formated list nammed data_Gratiot 
#' refdate <- as.Date("2001-01-01")
#' data_Gratiot<-add_phenology(Gratiot, name="Complete", reference=refdate, format="%d/%m/%Y")
#' # Fix parameter FLat to 0
#' pfixed=c(Flat=0)
#' # Generate initial points for the optimisation
#' parg<-par_init(data_Gratiot, parametersfixed=pfixed)
#' # Fit is done
#' result_Gratiot_Flat<-fit_phenology(data=data_Gratiot, parametersfit=parg2, 
#' 		parametersfixed=pfixed, trace=1)
#' data(result_Gratiot_Flat)
#' # Generate initial points for the optimisation
#' parg<-par_init(data_Gratiot, parametersfixed=NULL)
#' # Run the optimisation
#' result_Gratiot<-fit_phenology(data=data_Gratiot, 
#' 		parametersfit=parg, parametersfixed=NULL, trace=1)
#' data(result_Gratiot)
#' # Compare both models
#' outputAIC<-compare_AIC(full=result_Gratiot, Flat=result_Gratiot_Flat)
#' }
#' @export


compare_AIC <-
function(...) {
result <- list(...)

if (is.list(result) & length(result)==1) result <- unlist(result, recursive=FALSE)

if (!is.null(result)) {
	if (!is.list(result) || (is.null(names(result))) || (any(names(result)==""))) {
		print("The results must be included within a list with names; see example.")
		return(invisible())
	} else {
	out<-NULL
	l<-length(result)
		if (l<2) {
			print("A least two results must be provided.")
			return(invisible())
		} else {
			aic<-NULL
			name<-names(result)
			for (i in 1:l) {

        encours <- result[i]
        
        if (is.null(encours[[1]]$aic) & is.null(encours[[1]]$AIC) & is.null(encours[[1]]$value))
          encours <- encours[[1]]
        
 
# je n'ai pas de AIC ou rien pour le calculer
          sumAIC <- 0
            AICencours <- NULL
            for (j in 1:length(encours)) {
              encours2 <- encours[[j]]
              if (!is.null(encours2$AIC)) AICencours <- encours2$AIC
              if (!is.null(encours2$aic)) AICencours <- encours2$aic
              if (!is.null(encours2$value) & !is.null(encours2$par))
                AICencours <- 2*encours2$value+2*(length(encours2$par))
            if (is.null(AICencours))
              {
                print(paste("Object", name[i], "has not the required format"))
                return(invisible())
              }
            sumAIC <- sumAIC + AICencours
            }
			aic <- c(aic, sumAIC)	
			
			}
			
			bestaic<-min(aic)
			ser<-which(aic==bestaic)
			deltaaic<-aic-bestaic
			aw<-exp(-0.5*deltaaic)
			saw=sum(aw)
			aw<-aw/saw
			
			out<-data.frame(cbind(AIC=aic, DeltaAIC=deltaaic, Akaike_weight=aw), row.names=name)
			print(paste("The lowest AIC (",sprintf("%.3f", bestaic) ,") is for series ", name[ser], " with Akaike weight=", sprintf("%.3f", aw[ser]), sep=""))
			
			return(out)
		}
	}
}
}
