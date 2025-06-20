% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/add_phenology.R
\name{add_phenology}
\alias{add_phenology}
\title{Create a new dataset or add a timeserie to a previous dataset.}
\usage{
add_phenology(
  add = stop("New data must be given."),
  name = "Site",
  reference = NULL,
  month_ref = NULL,
  sep.dates = "-",
  end.season.date = NULL,
  colname.Date = 1,
  colname.Number = 2,
  colname.Rookery = "Site",
  colname.CountTypes = NULL,
  CountTypes.default = "exact",
  colname.A = NULL,
  A.default = NA,
  colname.S = NULL,
  S.default = NA,
  colname.ZeroCounts = NULL,
  ZeroCounts.default = TRUE,
  format = "\%d/\%m/\%Y",
  previous = NULL,
  include0 = FALSE,
  datepeakfor0 = NULL,
  expandRange0Observation = TRUE,
  check.overlapping.dates = TRUE,
  silent = FALSE
)
}
\arguments{
\item{add}{The data to be added. It can be a set of several entities that uses the same reference and date format.}

\item{name}{The name of the monitored site.}

\item{reference}{as.Date('2001-12-31') The date used as day 0 in the ordinal date model.}

\item{month_ref}{If no reference date is given, use this month as a reference.}

\item{sep.dates}{Separator used to separate dates when incertitude is included.}

\item{end.season.date}{The date corresponding to the end of nesting season.}

\item{colname.Date}{Name or number of column with dates.}

\item{colname.Number}{Name or number of column with numbers.}

\item{colname.Rookery}{Name or number of column with rookery names.}

\item{colname.CountTypes}{Model of count type. It can be "exact" (default), "minimum" or a number to indicate the maximum possible.}

\item{CountTypes.default}{The default of CountTypes if colname.CountTypes is not provided.}

\item{colname.A}{The A parameter of the detection model}

\item{A.default}{Default A value.}

\item{colname.S}{The S parameter of the detection model}

\item{S.default}{Default S value.}

\item{colname.ZeroCounts}{The name of the column to indicate whether zero counts are included (TRUE is default).}

\item{ZeroCounts.default}{The default for ZeroCounts.}

\item{format}{The format of the dates.}

\item{previous}{Previous data formated with add_phenology or NULL [default] if no previous data exist.}

\item{include0}{Does timeseries with only 0 should be included?}

\item{datepeakfor0}{If series with no observation are included, where add a 1 value in ordinal date (see description)}

\item{expandRange0Observation}{If TRUE, the range of date with 0 observations are expanded into individual dates}

\item{check.overlapping.dates}{If TRUE, will check for date overlapping}

\item{silent}{Does information about added timeseries is shown}
}
\value{
Return a list of formated data that can be used ith fit_phenology()
}
\description{
To create a new dataset, the syntaxe is :\cr
data <- add_phenology(add=newdata, name="Site", reference=as.Date('2001-12-31'),
format='\\%d/\\%m/\\%y')\cr\cr
To add a dataset to a previous one, the syntax is :\cr
data <- add_phenology(previous=previousdata, add=newdata, name='Site', \cr
reference=as.Date('2001-01-01'), format="\\%Y-\\%m-\\%d") \cr\cr
The dataset to be added must include 2 or 3 columns.\cr
The colname.Date included the dates in the format specified by
the parameter format. If the number of nests is known
for an exact date, then only one date must be indicated.\cr
If the number of nests is known for a range of date, the
first and last dates must be separated by a sep.dates character.\cr
For example: 1/2/2000-10/2/2000\cr
Note that date in the colname.Date column can be already formated and in this case
the parameter format is ignored.\cr\cr
The colname.Number includes the number of nests observed for
this date or this range of dates.\cr
The colname.Rookery is optional and includes the name of the rookeries.\cr\cr
If only two columns are indicated, the name can be indicated as
a parameter of the function with the parameter name. If no name is indicated,
the default name Site will be used, but take care, only one
rookery of this name can be used.\cr\cr
Several rookeries can be included in the same file but in this case
the rookery name is obligatory at the colname.Rookery column.\cr\cr
The model cannot be fitted if a timeseries has no observation because the trivial
solution is of course with max=0. The solution is to include a fake false observation at the closest
position of the peak, and then the estimated number of nests/tracks will be the estimated number - 1.\cr
If include0 is TRUE, then the series with no observation are included and one observation is added
at the monitored date the closest of datepeakfor0.\cr
The normal way to manage such a situation is as followed:\cr
1- Format data with include0 being FALSE\cr
2- Fit parameters using fdf <- fit_phenology()\cr
3- Format data with include0 being TRUE and datepeakfor0=fdf$par["Peak"]\cr
4- Fix previsouly fitted parameters using pfixed <- fdf$par\cr
5- Generate new set of parameters with par_init(data, fixed.parameters=pfixed)\cr
6- Run again fit_phenology()\cr\cr
Some problems that can occur:\cr
If a name is defined as a third column of a data.frame and a name is
defined also with name parameter, the third column has priority.\cr
Two different timeseries MUST have different name and characters _ and
space are forbiden in timeseries names. They are automatically changed if they are present.
}
\details{
add_phenology creates a new dataset.
}
\examples{
\dontrun{
library(phenology)
# Read a file with data
data(Gratiot)
# Generate a formatted list nammed data_Gratiot 
refdate <- as.Date("2001-01-01")
data_Gratiot <- add_phenology(Gratiot, name="Complete", 
	reference=refdate, format="\%d/\%m/\%Y")
	
# Generate initial points for the optimisation
parg <- par_init(data_Gratiot, fixed.parameters=NULL)
# Run the optimisation
result_Gratiot <- fit_phenology(data=data_Gratiot, fitted.parameters=parg, 
	fixed.parameters=NULL)
data(result_Gratiot)
# Plot the phenology and get some stats
output <- plot(result_Gratiot)

#############################################
# Example of use of include0 and datepeakfor0
#############################################
# Let create a times series with only 0
data0 <- data.frame(Date=c("11/3/2015", "12/3/2015", "13/3/2015-18/3/2015", "25/3/2015"), 
                    Number=c(0, 0, 0, 0), 
                    Beach=rep("Site", 4), stringsAsFactors=FALSE)
# Here I don't include beach with no observation: error message
try1 <- add_phenology(data0, format="\%d/\%m/\%Y", month_ref=1, include0=FALSE)
# Here I include timeseries with no observation
try1 <- add_phenology(data0, format="\%d/\%m/\%Y", month_ref=1, include0=TRUE, datepeakfor0=100)
try1 <- add_phenology(data0, format="\%d/\%m/\%Y", month_ref=1, include0=TRUE, datepeakfor0=73)
try1 <- add_phenology(data0, format="\%d/\%m/\%Y", month_ref=1, include0=TRUE, datepeakfor0=70)
# It can be done in two steps
try1 <- add_phenology(data0, format="\%d/\%m/\%Y", month_ref=1, include0=TRUE)
try2 <- add_phenology(previous=try1, include0=TRUE, datepeakfor0=100)
# Here I include the series without observation
try1 <- add_phenology(add=data0, format="\%d/\%m/\%Y", month_ref=1, 
                      include0=TRUE, expandRange0Observation=TRUE)
                      
##############################################################################
# Example of A and S parameters to say that only half of a beach was monitored
##############################################################################

refdate <- as.Date("2001-01-01")
data_Gratiot <- add_phenology(Gratiot, name="Complete1", 
                             	reference=refdate, format="\%d/\%m/\%Y")
	      S.default = 10
	      # Let Complete1 be of length 2 and Complete0.5 be of length 1
	      length.Complete1 <- 2
	      length.Complete0.5 <- 1
	      A.default = log(1/(length.Complete0.5/length.Complete1)-1)/(-S.default*4)
	      # For day 0, the detection probability is
	      1/(1+exp(-(4*S.default)*(A.default-0)))
	data_Gratiot <- add_phenology(previous=data_Gratiot, 
	      add=Gratiot, 
	      name="Complete0.5", 
	      A.default = A.default,
	      S.default = S.default, 
	      reference=refdate, 
	      format="\%d/\%m/\%Y")
	parg <- c(par_init(data_Gratiot, fixed.parameters=c(Min=0, Flat=0)), PMin=0.1)
	result_Gratiot <- fit_phenology(data=data_Gratiot, fitted.parameters=parg, 
	                                fixed.parameters=c(Flat=0))
	result_Gratiot$par
	# it shows that Max_Complete1 is half of Max_Complete0.5; all is ok
	summary(result_Gratiot)$synthesis
	
}
}
\seealso{
Other Phenology model: 
\code{\link{AutoFitPhenology}()},
\code{\link{BE_to_LBLE}()},
\code{\link{Gratiot}},
\code{\link{LBLE_to_BE}()},
\code{\link{LBLE_to_L}()},
\code{\link{L_to_LBLE}()},
\code{\link{MarineTurtles_2002}},
\code{\link{MinBMinE_to_Min}()},
\code{\link{adapt_parameters}()},
\code{\link{add_SE}()},
\code{\link{extract_result}()},
\code{\link{fit_phenology}()},
\code{\link{likelihood_phenology}()},
\code{\link{logLik.phenology}()},
\code{\link{map_Gratiot}},
\code{\link{map_phenology}()},
\code{\link{par_init}()},
\code{\link{phenology}()},
\code{\link{phenology2fitRMU}()},
\code{\link{phenology_MHmcmc}()},
\code{\link{phenology_MHmcmc_p}()},
\code{\link{plot.phenology}()},
\code{\link{plot.phenologymap}()},
\code{\link{plot_delta}()},
\code{\link{plot_phi}()},
\code{\link{print.phenology}()},
\code{\link{print.phenologymap}()},
\code{\link{print.phenologyout}()},
\code{\link{remove_site}()},
\code{\link{result_Gratiot}},
\code{\link{result_Gratiot1}},
\code{\link{result_Gratiot2}},
\code{\link{result_Gratiot_Flat}},
\code{\link{summary.phenology}()},
\code{\link{summary.phenologymap}()},
\code{\link{summary.phenologyout}()}
}
\author{
Marc Girondot \email{marc.girondot@gmail.com}
}
\concept{Phenology model}
