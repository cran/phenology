#' add_format creates a new dataset.
#' @title Create a new dataset or add a timeserie to a previous one.
#' @author Marc Girondot
#' @return Return a list of formated data
#' @param origin previousdata or NULL if no previous data exists
#' @param add newdata The data to be added
#' @param name 'Site' The name of the monitored site
#' @param reference as.Date('2001-12-31') The date used as 1st date
#' @param format The format of the date in the file
#' @param help If TRUE, an help is displayed
#' @description To create a new dataset, the syntaxe is \cr
#' data<-add_format(add=newdata, name="Site", reference=as.Date('2001-12-31'), format='%d/%m/%y')\cr
#' To add a dataset to a previous one, the syntaxe is \cr
#' data<-add_format(origin=previousdata, add=newdata, name='Site', reference=as.Date('2001-12-31'), format='%d/%m/%y')\cr
#' \cr
#' The dataset to be added must include 2 or 3 columns.\cr
#' The first one is the date in the format specified by 
#' the parameter format=. If the number of nests is known 
#' for an exact data, then only one date must be indicated 
#' If the number of nests is known for a range of date, the 
#' first and last dates must be separated but a -.\cr
#' For example: 1/2/2000-10/2/2000\cr
#' The second column is the number of nests observed for 
#' this date or this range of dates.\cr
#' The third column is optional and is the name of the rookery.\cr
#' If only two columns are indicated, the name can be indicated as 
#' a parameter of the function with name=. If no name is indicated, 
#' the default name Site will be used, but take care, only one 
#' rookery this name can be used.\cr
#' Several rookeries can be included in the same file but in this case 
#' the rookery name is obligatory at the third column.
#' @examples 
#' library(phenology)
#' # Read a file with data
#' # Gratiot<-read.delim("http://max2.ese.u-psud.fr/epc/conservation/BI/Complete.txt", , header=FALSE)
#' data(Gratiot)
#' # Generate a formatted list nammed data_Gratiot 
#' data_Gratiot<-add_format(origin=NULL, add=Gratiot, name="Complete", reference=as.Date("2001-01-01"), format="%d/%m/%Y")
#' # Generate initial points for the optimisation
#' parg<-par_init(data_Gratiot, parametersfixed=NULL)
#' # Run the optimisation
#' # result_Gratiot<-fit_phenology(data=data_Gratiot, parametersfit=parg, parametersfixed=NULL, trace=1)
#' data(result_Gratiot)
#' # Plot the phenology and get some stats
#' plot_phenology(result=result_Gratiot, pdf=FALSE)
#' @export


add_format <-
function(origin=NULL, add=stop("A dataset must be indicated"), name=NULL, reference=NULL, format="%d/%m/%y", help=FALSE) {
if(help) {
	cat("To create a new dataset, the syntaxe is \n")
	cat("data<-add_format(add=newdata, name='Site', \n")
	cat("+   reference=as.Date('2001-12-31'), format='%d/%m/%y')\n")
	cat("To add a dataset to a previous one, the syntaxe is \n")
	cat("data<-add_format(origin=previousdata, add=newdata, name='Site', \n")
	cat("+   reference=as.Date('2001-12-31'), format='%d/%m/%y')\n")
	cat("\n")
	cat("The dataset to be added must include 2 or 3 columns.\n")
	cat("The first one is the date in the format specified by\n")
	cat("the parameter format=. If the number of nests is known \n")
	cat("for an exact data, then only one date must be indicated\n")
	cat("If the number of nests is known for a range of date, the\n")
	cat("first and last dates must be separated but a -.\n")
	cat("For example: 1/2/2000-10/2/2000\n")
	cat("The second column is the number of nests observed for\n")
	cat("this date or this range of dates.\n")
	cat("The third column is optional and is the name of the rookery.\n")
	cat("If only two columns are indicated, the name can be indicated as \n")
	cat("a parameter of the function with name=. If no name is indicated,\n")
	cat("the default name Site will be used, but take care, only one \n")
	cat("rookery this name can be used.\n")
	cat("Several rookeries can be included in the same file but in this case\n")
	cat("the rookery name is obligatory at the third column.\n")
} else {

	if (is.null(name)) {name<-deparse(substitute(add))}
	
# fichier pnd
# si 5 colonnes, j'en retire les 2 dernières
	if (dim(add)[2]==5) {
		add<-add[,-c(3:5)]
	}
#Si j'ai 3 colonnes et rien sur la 3ème, je n'en garde que deux
	if (dim(add)[2]==3) {
		if (all(is.na(add[,3]))) {
			add<-add[,-3]
		}
	}

	if (dim(add)[2]==2) {
	
# J'ai deux colonnes et le nom des dans name
	print(name)
# Je n'ai pas de nom de site. Il n'y a donc qu'une seule série
		colnames(add)=c("Date", "Nombre")	
		add$Date<-as.character(add$Date)
		add[, 2]<-as.character(add[, 2])
# je retire toutes les lignes pour lesquelles je n'ai pas d'observation de ponte - 19012012
		add<-add[(add[,1]!="") & (!is.na(add[,1])) & (gsub("[0-9 ]", "", add[,2])=="") & (!is.na(add[,2])) & (add[,2]!=""),1:2]
		
		if (dim(add)[1]==0) {
			print(paste("The timeseries", name, "is empty; check it"))
		} else {

		
		# Je le mets en numérique: 21/3/2012
		add[, 2]<-as.numeric(add[, 2])
		addT<-data.frame(Date=rep(as.Date("1900-01-01"), dim(add)[1]), Date2=rep(as.Date("1900-01-01"), dim(add)[1]), nombre=rep(NA, dim(add)[1]), ordinal=rep(NA, dim(add)[1]), ordinal2=rep(NA, dim(add)[1]), Modeled=rep(NA, dim(add)[1]), LnL=rep(NA, dim(add)[1]))
		for(i in 1:dim(add)[1]) {
			essai<-unlist(strsplit(add$Date[i], "-"))
			addT$Date[i]<-as.Date(essai[1], format)
			if (is.na(addT$Date[i])) {
				print(paste("Error date ", essai[1], sep=""))
			} else {
				addT$ordinal[i]<-as.numeric(addT$Date[i]-reference+1)
				if (length(essai)==2) {
					addT$Date2[i]<-as.Date(essai[2], format)
					addT$ordinal2[i]<-as.numeric(addT$Date2[i]-reference+1)
				} else {
					addT$Date2[i]<-NA
					addT$ordinal2[i]<-NA
				}
			}
		}
		addT$nombre<-add$Nombre
		nb<-length(origin)
		origin$kyYI876Uu<-addT
		names(origin)[nb+1]<-name
		}
		
#		return(origin)
	} else {
# J'ai un nom de site. Il n'y a donc peut-être plusieurs séries - 21012012
		colnames(add)=c("Date", "Nombre", "Site")	
		add$Date<-as.character(add$Date)
		add[, 2]<-as.character(add[, 2])
# je retire toutes les lignes pour lesquelles je n'ai pas d'observation de ponte - 19012012
		add<-add[(add[,1]!="") & (!is.na(add[,1])) & (gsub("[0-9 ]", "", add[,2])=="") & (!is.na(add[,2])) & (add[,2]!=""),1:3]

# Je le mets en numérique: 21/3/2012
		add[, 2]<-as.numeric(add[, 2])
		
# il faut que je mette les noms sur les sites sans
		siteencours=""
		for(i in 1:dim(add)[1]) {
			if (add$Site[i]=="") {add$Site[i]<-siteencours} else {siteencours<-add$Site[i]}
		}
# Maintenant j'ai tous les sites remplis
# mais je dois avoir la liste des sites
		fac <- factor(add$Site)
# je génère un nouveau data.frame avec facteur par facteur
		dtaorigin=origin
		for(i in 1:length(levels(fac))) {
			siteencours<-levels(fac)[i]
			addencours<-add[add[,3]==siteencours, 1:2]
			dtaorigin<-add_format(origin=dtaorigin, add=addencours, name=siteencours, reference=reference, format=format)
		}
		origin<-dtaorigin

	}
	problem<-FALSE
	for(i in 1:length(origin)) {
		problem<-(problem) || (any(origin[[i]]$ordinal[!is.na(origin[[i]]$ordinal)]>365)) || (any(origin[[i]]$ordinal2[!is.na(origin[[i]]$ordinal2)]>365)) || (any(origin[[i]]$ordinal[!is.na(origin[[i]]$ordinal)]<0)) || (any(origin[[i]]$ordinal2[!is.na(origin[[i]]$ordinal2)]<0))

	}
#	print(problem)
	if (problem) {
		cat("Problem in at least one date; check them. Take care about format used.\n")
		cat("Within a file, all dates must conform to the same format.\n")
	}	
	
	return(origin)
}
}
