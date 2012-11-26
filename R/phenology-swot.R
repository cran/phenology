#' phenology_swot is a simplified function for phenology package.
#' @title Fit a nesting season of marine turtles
#' @author Marc Girondot \email{marc.girondot@@u-psud.fr}
#' @return None
#' @param header Does the timeseries has header.
#' @param reference Date used as reference - day 1.
#' @param format Format for dates.
#' @description This function tries to make very simple the use of this package.\cr
#' In most of the cases, values for header, reference and format need not to be set\cr
#' because they are detected automically.
#' @export


phenology_swot <-
function(header=NULL, reference=NULL, format=NULL) {

	premier <- TRUE
	Formated <- NULL
	
			repeat {

			if (premier) {
		 		cat("Choose the file to be read:\n")
		 	} else {
		 		cat("Choose another file to be read: Cancel to begin analysis\n")
			}
			
			nm <- try(file.choose(), silent=TRUE)
			if(class(nm)=="try-error") {
				if (premier) return(invisible())
				break()
			}
			
			premier <- FALSE
			
				if (is.null(header)) header <- FALSE
				

				obj_list <-  lapply(nm,readLines, warn=FALSE)
				

				analyse <- unlist(obj_list)[2]
				if (length(grep("\t", analyse))!=0) {
					spl <- "\t"
				
				
				} else {
					if (length(grep(";", analyse))!=0) {
						spl <- ";"
					
					
					} else {
						spl <- ","
					}
				}
				dta <- unlist(strsplit(analyse, spl))[1]
				
# Je teste / et -
				dtaspl <- "/"
				if (length(grep("/", dta))!=0) {dtaspl <- "/"}
				if (length(grep("-", dta))!=0) {dtaspl <- "-"}
				
				analyse1 <- unlist(obj_list)[1]
				dta1 <- unlist(strsplit(analyse1, spl))[1]

				if (length(grep(dtaspl, dta1))!=0) {header <- FALSE} else {header <- TRUE}
				DATA <- read.table(nm, header = header, sep = spl, stringsAsFactors = FALSE)
			# DATA[,1] j'ai les dates
			
				es <- as.data.frame(strsplit(DATA[,1], dtaspl), stringsAsFactors=FALSE)
				
				y <- "%y"
				an <- max(c(max(nchar(es[1,])), max(nchar(es[2,])), max(nchar(es[3,]))))
				if (an==4) y <- "%Y"
				
				nf <- length(levels(as.factor(as.numeric(es[1,]))))
				nf <- c(nf, length(levels(as.factor(as.numeric(es[2,])))))
				nf <- c(nf, length(levels(as.factor(as.numeric(es[3,])))))
				
				ypos <- which.min(nf)
				dpos <- which.max(nf)
				mpos <- c(1:3)[-c(ypos, dpos)]
				
				if (is.null(format)) {
				
				if (ypos == 1) {format <- y} else {
					if (dpos==1) {format <- "%d"} else {format <- "%m"}}
				format <- paste(format, dtaspl, sep="")
				if (ypos == 2) {format <- paste(format, y, sep="") } else {
					if (dpos==2) {format <- paste(format, "%d", sep="")} else {format <- paste(format, "%m", sep="")}}
				format <- paste(format, dtaspl, sep="")
				if (ypos == 3) {format <- paste(format, y, sep="") } else {
					if (dpos==3) {format <- paste(format, "%d", sep="")} else {format <- paste(format, "%m", sep="")}}
				
				}
				
				if (is.null(reference)) {
				
				year <- min(as.numeric(es[ypos,]))
				if (year == 0 & max(as.numeric(es[ypos,])) == 99) year=1999
				if (year<100) {
					if (year<50) {year <- year + 2000} else {year <- year + 1900}
				}
						
				year <- as.character(year)
				
				if (nf[ypos] ==2) {
					# 2 années donc juin
					reference <- as.Date(paste(year, "-06-01", sep=""))
				} else {
					reference <- as.Date(paste(year, "-01-01", sep=""))
				}
				
				}
					
				Formated <- add_format(origin=Formated, add=DATA, reference=reference, format=format, name=basename(nm))

}

				pfixed <- c(Min=0, Flat=0)
				Par <- par_init(data=Formated, parametersfixed=pfixed)
				print("Test 1/4: Please be patient...")
				result1 <- fit_phenology(data=Formated, parametersfit = Par, parametersfixed=pfixed, trace=0, silent=TRUE)
				pfixed <- c(Flat=0)
				Par <- par_init(data=Formated, parametersfixed=pfixed)
				print("Test 2/4: Please be patient...")
				result2 <- fit_phenology(data=Formated, parametersfit = Par, parametersfixed=pfixed, trace=0, silent=TRUE)
				pfixed <- c(Min=0)
				Par <- par_init(data=Formated, parametersfixed=pfixed)
				print("Test 3/4: Please be patient...")
				result3 <- fit_phenology(data=Formated, parametersfit = Par, parametersfixed=pfixed, trace=0, silent=TRUE)
				pfixed <- NULL
				Par <- par_init(data=Formated, parametersfixed=pfixed)
				print("Test 4/4: Please be patient...")
				result4 <- fit_phenology(data=Formated, parametersfit = Par, parametersfixed=pfixed, trace=0, silent=TRUE)

				result <- list(result1, result2, result3, result4)[which.min(c(2*(result1$value+length(result1$par)), 2*(result2$value+length(result2$par)), 2*(result3$value+length(result3$par)), 2*(result4$value+length(result4$par))))][[1]]

				plot_phenology(result=result, series="all", moon=FALSE)
			

}
