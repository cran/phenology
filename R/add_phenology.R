#' add_phenology creates a new dataset.
#' @title Create a new dataset or add a timeserie to a previous dataset.
#' @author Marc Girondot
#' @return Return a list of formated data that can be used ith fit_phenology()
#' @param previous Previous data formated with add_phenology or NULL [default] if no previous data exist
#' @param add The data to be added. It can be a set of several entities that uses the same reference and date format
#' @param name The name of the monitored site
#' @param reference as.Date('2001-12-31') The date used as 1st date
#' @param sep.dates Separator used to separate dates when incertitude is included
#' @param month_ref If no reference date is given, use this month as a reference
#' @param format The format of the dates
#' @param colname.Date Name or number of column with dates
#' @param colname.Number Name or number of column with numbers
#' @param colname.Rookery Name or number of column with rookery names
#' @param include0 Does timeseries with only 0 should be included?
#' @param datepeakfor0 If series with no observation are included, where add a 1 value in ordinal date (see description)
#' @param expandRange0Observation If TRUE, the range of date with 0 observations are expanded into individual dates
#' @param silent Does information about added timeseries is shown
#' @description To create a new dataset, the syntaxe is :\cr
#' data <- add_phenology(add=newdata, name="Site", reference=as.Date('2001-12-31'), 
#' format='\%d/\%m/\%y')\cr\cr
#' To add a dataset to a previous one, the syntaxe is :\cr
#' data <- add_phenology(previous=previousdata, add=newdata, name='Site', \cr
#' reference=as.Date('2001-01-01'), format="\%Y-\%m-\%d") \cr\cr
#' The dataset to be added must include 2 or 3 columns.\cr
#' The colname.Date included the dates in the format specified by 
#' the parameter format. If the number of nests is known  
#' for an exact date, then only one date must be indicated.\cr  
#' If the number of nests is known for a range of date, the 
#' first and last dates must be separated but a sep.dates character.\cr
#' For example: 1/2/2000-10/2/2000\cr
#' Note that date in the colname.Date column can be already formated and in this case 
#' the parameter format is ignored.\cr\cr
#' The colname.Number includes the number of nests observed for 
#' this date or this range of dates.\cr
#' The colname.Rookery is optional and includes the name of the rookeries.\cr\cr
#' If only two columns are indicated, the name can be indicated as  
#' a parameter of the function with the parameter name. If no name is indicated,  
#' the default name Site will be used, but take care, only one 
#' rookery of this name can be used.\cr\cr
#' Several rookeries can be included in the same file but in this case 
#' the rookery name is obligatory at the colname.Rookery column.\cr\cr
#' #' The model cannot be fitted if a timeseries has no observation because the trivial 
#' solution is of course with max=0. The solution is to include a fake false observation at the closest 
#' position of the peak, and then the estimated number of nests/tracks will be the estimated number - 1.\cr
#' If include0 is TRUE, then the series with no observation are included and one observation is added 
#' at the monitored date the closest of datepeakfor0.\cr
#' The normal way to manage such a situation is as followed:\cr
#' 1- Format data with include0 being FALSE\cr
#' 2- Fit parameters using fdf <- fit_phenology()\cr
#' 3- Format data with include0 being TRUE and datepeakfor0=fdf$par["Peak"]\cr
#' 4- Fix previsouly fitted parameters using pfixed <- fdf$par\cr
#' 5- Generate new set of parameters with par_init(data, fixed.parameters=pfixed)\cr
#' 6- Run again fit_phenology()\cr\cr
#' Some problems that can occur:\cr
#' If a name is defined as a third column of a data.frame and a name is 
#' defined also with name parameter, the third column has priority.\cr
#' Two different timeseries MUST have different name and characters _ and 
#' space are forbiden in timeseries names. They are automatically changed if they are present.
#' @family Phenology model
#' @examples
#' \dontrun{
#' library(phenology)
#' # Read a file with data
#' data(Gratiot)
#' # Generate a formatted list nammed data_Gratiot 
#' refdate <- as.Date("2001-01-01")
#' data_Gratiot <- add_phenology(Gratiot, name="Complete", 
#' 	reference=refdate, format="%d/%m/%Y")
#' 	
#' # Generate initial points for the optimisation
#' parg <- par_init(data_Gratiot, fixed.parameters=NULL)
#' # Run the optimisation
#' result_Gratiot <- fit_phenology(data=data_Gratiot, fitted.parameters=parg, 
#' 	fixed.parameters=NULL)
#' data(result_Gratiot)
#' # Plot the phenology and get some stats
#' output <- plot(result_Gratiot)
#' 
#' #############################################
#' # Example of use of include0 and datepeakfor0
#' #############################################
#' # Let create a times series with only 0
#' data0 <- data.frame(Date=c("11/3/2015", "12/3/2015", "13/3/2015-18/3/2015", "25/3/2015"), 
#'                     Number=c(0, 0, 0, 0), 
#'                     Beach=rep("Site", 4), stringsAsFactors=FALSE)
#' # Here I don't include beach with no observation: error message
#' try1 <- add_phenology(data0, format="%d/%m/%Y", month_ref=1, include0=FALSE)
#' # Here I include timeseries with no observation
#' try1 <- add_phenology(data0, format="%d/%m/%Y", month_ref=1, include0=TRUE, datepeakfor0=100)
#' try1 <- add_phenology(data0, format="%d/%m/%Y", month_ref=1, include0=TRUE, datepeakfor0=73)
#' try1 <- add_phenology(data0, format="%d/%m/%Y", month_ref=1, include0=TRUE, datepeakfor0=70)
#' # It can be done in two steps
#' try1 <- add_phenology(data0, format="%d/%m/%Y", month_ref=1, include0=TRUE)
#' try2 <- add_phenology(previous=try1, include0=TRUE, datepeakfor0=100)
#' # Here I include the series without observation
#' try1 <- add_phenology(add=data0, format="%d/%m/%Y", month_ref=1, 
#'                       include0=TRUE, expandRange0Observation=TRUE)
#' }
#' @export


add_phenology <-
  function(add=NULL, name="Site", reference=NULL, 
           month_ref= NULL, sep.dates="-", 
           colname.Date=1, colname.Number=2, colname.Rookery=3, 
           format="%d/%m/%Y", previous=NULL, include0=FALSE, datepeakfor0=NULL, 
           expandRange0Observation=TRUE, 
           silent=FALSE) {
    
    
    # name=NULL; reference=NULL; month_ref= NULL; sep.dates="-"; format="%d/%m/%Y"; previous=NULL; colname.Date=1; colname.Number=2; colname.Rookery=3; silent=FALSE; include0=FALSE; datepeakfor0=NULL; expandRange0Observation=TRUE
    
    if (class(previous) != "phenologydata" && !is.null(previous)) {
      stop("The previous dataset must be already formated using add_phenology()!")
    }
    
    if (!is.null(add)) {
      
      if (!is.null(format)) {
        if (any(grepl(sep.dates, format))) {
          stop("Separator between day, month and year cannot be the same as the separator between two dates")
        }
      }
      
      if (class(try(add[, colname.Date], silent = TRUE)) == "try-error") {
        stop("The columns with dates does not exist")
      }
      
      if (class(try(add[, colname.Number], silent = TRUE)) == "try-error") {
        stop("The columns with numbers does not exist")
      }
      
      if (class(try(add[, colname.Rookery], silent = TRUE)) == "try-error") {
        if (!silent) message("The columns with rookery name does not exist; I create one")
        add <- cbind(add, Site=rep(name, nrow(add)))
      }
      
      if (is.null(reference) & is.null(month_ref)) {
        stop("reference or month_ref must be supplied")
      }
      
      if (is.factor(add[,colname.Date])) add[,colname.Date] <- as.character(add[,colname.Date])
      if (is.factor(add[,colname.Rookery])) add[,colname.Rookery] <- as.character(add[,colname.Rookery])
      
      d <- add[,colname.Date]
      if (is.character(d)) {
        d2 <- strsplit(d, split=sep.dates)
        dp1 <- as.Date(unlist(lapply(d2, FUN = function(x) x[1])), format=format)
        dp2 <- as.Date(unlist(lapply(d2, FUN = function(x) x[2])), format=format)
      } else {
        dp1 <- d
        dp2 <- rep(NA, length(d))
      }
      add <- cbind(add, D18989898=dp1, D28989898=dp2)
      
      if (any(is.na(dp1))) {
        stop(paste("Date format for column colname.Date is not correct; check line(s)", which(is.na(dp1))))
      }
      
      if (is.character(reference)) reference <- as.Date(reference, format=format)
      
      intermediaire <- list()
      
      for (site in unique(add[, colname.Rookery])) {
        if (!silent) message(paste0("Site: ", site))
        dfadd <- add[add[, colname.Rookery] == site, ]
        dfadd <- dfadd[order(dfadd[, "D18989898"]), ]
        
        if (is.null(reference)) {
          # si month_ref > premier mois de la série, c'est l'année n-1
          # sinon c'est l'année n
          premieredate <- dfadd[1, "D18989898"]
          premiermois <- as.POSIXlt(premieredate)$mon+1
          
          if (premiermois < month_ref) {
            premieredate <- premieredate - ifelse(as.POSIXlt(premieredate-365)$mday==as.POSIXlt(premieredate)$mday, 365, 366)
            refencours <- as.character(premieredate)
            substr(refencours, 9, 10) <- "01"
            m <- paste0("0", as.character(month_ref))
            substr(refencours, 6, 7) <- substr(m, nchar(m)-2, nchar(m))
            refencours <- as.Date(refencours)
          } else {
            refencours <- as.character(premieredate)
            substr(refencours, 9, 10) <- "01"
            m <- paste0("0", as.character(month_ref))
            substr(refencours, 6, 7) <- substr(m, nchar(m)-2, nchar(m))
            refencours <- as.Date(refencours)
          }
        } else {
          refencours <- reference
        }
        
        if (!silent) message(paste("Reference date:", as.character(refencours)))
        
        df <- data.frame(Date=dfadd$D18989898, 
                         Date2=dfadd$D28989898, 
                         nombre=dfadd[, colname.Number], 
                         ordinal=as.numeric(dfadd$D18989898-refencours), 
                         ordinal2=as.numeric(dfadd$D28989898-refencours), 
                         Modeled=rep(NA, nrow(dfadd)), 
                         LnL=rep(NA, nrow(dfadd)), 
                         stringsAsFactors = FALSE)
        if ((expandRange0Observation) & (any(!is.na(df[, "Date2"])))) {
          dfec <- df[-(1:nrow(df)), ]
          
          for (idf in 1:nrow(df)) {
            if ((!is.na(df[idf, "Date2"])) & (df[idf, "nombre"] == 0)) {
              dt <- seq(from=df[idf, "Date"], to=df[idf, "Date2"], by="1 day")
              dfe0 <- data.frame(Date=dt, 
                                 Date2=rep(as.Date(NA), length(dt)), 
                                 nombre=rep(0, length(dt)), 
                                 ordinal=as.numeric(dt-refencours), 
                                 ordinal2=rep(NA, length(dt)), 
                                 Modeled=rep(NA, length(dt)), 
                                 LnL=rep(NA, length(dt)), 
                                 stringsAsFactors = FALSE)
              dfec <- rbind(dfec, dfe0)
            } else {
              dfec <- rbind(dfec, df[idf, ])
            }
          }
          df <- dfec
        }
        
        
        attributes(df) <- modifyList(attributes(df), list(reference=refencours))
        
        lg <- ifelse(as.POSIXlt(refencours+365)$mday == as.POSIXlt(refencours)$mday, 365, 366)
        
        # test cohérence
        if (max(c(df$ordinal, df$ordinal2)+1, na.rm = TRUE) > lg) {
          stop("More than one year for this series")
        }
        
        test <- logical(lg)
        for (i in 1:nrow(df)) {
          if (is.na(df[i, "ordinal2"])) {
            if (test[df[i, "ordinal"]+1]) {
              stop(paste("Error; some dates are duplicated. Look at", as.character(df[i, "Date"])))
            } else {
              test[df[i, "ordinal"]+1] <- TRUE
            }
          } else {
            if (any(test[(df[i, "ordinal"]:df[i, "ordinal2"])+1])) {
              stop(paste("Error; some dates are duplicated. Look at", as.character(df[i, "Date"])))
            } else {
              test[(df[i, "ordinal"]:df[i, "ordinal2"])+1] <- TRUE
            }
          }
        }
        df <- list(df)
        names(df) <- site
        intermediaire <- c(intermediaire, df)
        
      }
      
      previous <- c(previous, intermediaire)
      class(previous) <- "phenologydata"
      
    }
    
    if (length(previous) == 0) {
      warning("No timeseries is included.")
    } else {
      
      
      if (any(grepl("_", names(previous)))) {
        if (!silent) print("The character _ is forbiden in names of timesseries. It has been changed to '.'.")
        names(previous) <- gsub("_", ".", names(previous))
      }
      
      if (any(grepl("-", names(previous)))) {
        if (!silent) print("The character - in timeseries must be used to separate beach name and year.")
        # names(previous) <- gsub("-", ".", names(previous))
      }
      
      if (any(grepl(" ", names(previous)))) {
        if (!silent) print("The character ' ' is forbiden in names of timesseries. It has been changed to '.'.")
        names(previous) <- gsub(" ", ".", names(previous))
      }
      
      if (any(duplicated(names(previous)))) {
        stop("The names of timesseries must be unique.")
      }
      
      intermediaire <- list()
      for (site in names(previous)) {
        
        df <- previous[[site]]
        if (sum(df$nombre, na.rm = TRUE) == 0) {
          if (!silent) message(paste0("Site: ", site))
          
          if (include0) {
            if (is.null(datepeakfor0)) {
              if (!silent) message("This series contains no observation; datepeakfor0 is not defined.")
              attributes(df) <- modifyList(attributes(df), list(fakeobservation=NULL))
              
            } else {
              # datepeakfor0
              pos <- NULL
              dispos_min <- +Inf
              pos_min <- NULL
              
              for (i in 1:nrow(df)) {
                if (!is.na(df[i, "ordinal2"])) {
                  if ((datepeakfor0 >= df[i, "ordinal"]) & (datepeakfor0 <= df[i, "ordinal2"])) {
                    pos <- i
                  } else {
                    if (abs(datepeakfor0 - df[i, "ordinal"]) < dispos_min) {
                      dispos_min <- abs(datepeakfor0 - df[i, "ordinal"])
                      pos_min <- i
                    }
                    if (abs(datepeakfor0 - df[i, "ordinal2"]) < dispos_min) {
                      dispos_min <- abs(datepeakfor0 - df[i, "ordinal2"])
                      pos_min <- i
                    }
                    
                  }
                } else {
                  if (df[i, "ordinal"] == datepeakfor0) {
                    pos <- i
                  } else {
                    if (abs(datepeakfor0 - df[i, "ordinal"]) < dispos_min) {
                      dispos_min <- abs(datepeakfor0 - df[i, "ordinal"])
                      pos_min <- i
                    }
                  }
                }
                if (!is.null(pos)) break
              }
              
              if (is.null(pos)) pos <- pos_min
              df[pos, "nombre"] <- 1
              if (!silent) message("This series contains no observation; one fake observation was created")
              attributes(df) <- modifyList(attributes(df), list(fakeobservation=i))
            }
            
            df <- list(df)
            names(df) <- site
            intermediaire <- c(intermediaire, df)
            
            
          } else {
            if (!silent) message("This series contains no observation; it is not included")
          }
        } else {
          df <- list(df)
          names(df) <- site
          intermediaire <- c(intermediaire, df)
          
        }
        
      }
      
      previous <- intermediaire
      class(previous) <- "phenologydata"
    }
    
    return(previous)
  }
