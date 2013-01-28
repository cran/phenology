#' phenology_swot is a simplified function for phenology package developped as part of SWOT project
#' @title Fit a nesting season of marine turtles
#' @author Marc Girondot \email{marc.girondot@@u-psud.fr}
#' @return Return a list of with data and result
#' @param header Does the timeseries has header.
#' @param reference Date used as reference. Is the day 1.
#' @param month_ref Reference month. Generally will be 1 (January) or 7 (July).
#' @param format Format for dates.
#' @description This function tries to make very simple the use of this package.\cr
#' In most of the cases, values for header, reference and format need not to be set\cr
#' because they are detected automatically.
#' @export


phenology_swot <-
function(header=NULL, reference=NULL, month_ref= NULL, format=NULL) {

#  header=NULL; reference=NULL; month_ref= NULL; format=NULL
	premier <- TRUE
	Formated <- NULL
	
		repeat {

			if (premier) {
		 		cat("Choose the file to be read:\n")
		 	} else {
		 		cat("Choose another file to be read or press Cancel to begin the analysis\n")
			}
			
			nm <- try(file.choose(), silent=TRUE)
			if(class(nm)=="try-error") {
				if (premier) return(invisible())
				break()
			}
			
			premier <- FALSE
			
				obj_list <-  lapply(nm,readLines, warn=FALSE)
			  obj_listec <- obj_list[[1]]
			  obj_list[[1]] <- obj_listec[obj_listec!=""]
			
				rp <- .read_phenology(obj_list, header, reference, month_ref, format, nm)

					
				Formated <- add_phenology(previous=Formated, add=rp$DATA, reference=rp$reference, format=rp$format)
				

	}

				pfixed <- c(Min=0, Flat=0)
				Par <- par_init(data=Formated, parametersfixed=pfixed)
				print("Test 1/4: Please be patient...")
				result1 <- fit_phenology(data=Formated, parametersfit = Par, parametersfixed=pfixed, trace=0, silent=TRUE)
				x <- result1$par
				pfixed <- c(Flat=0)
				Par <- c(x, Min=0.1)
				print("Test 2/4: Please be patient...")
				result2 <- fit_phenology(data=Formated, parametersfit = Par, parametersfixed=pfixed, trace=0, silent=TRUE)
				pfixed <- c(Min=0)
				Par <- c(x, Flat=1)
				print("Test 3/4: Please be patient...")
				result3 <- fit_phenology(data=Formated, parametersfit = Par, parametersfixed=pfixed, trace=0, silent=TRUE)
				pfixed <- NULL
				Par <- c(x, Flat=1, Min=0.1)
				print("Test 4/4: Please be patient...")
				result4 <- fit_phenology(data=Formated, parametersfit = Par, parametersfixed=pfixed, trace=0, silent=TRUE)

				result <- list(result1, result2, result3, result4)[which.min(c(2*(result1$value+length(result1$par)), 2*(result2$value+length(result2$par)), 2*(result3$value+length(result3$par)), 2*(result4$value+length(result4$par))))][[1]]

				out <- plot(result, series="all", moon=FALSE)
				
				out
			

}
