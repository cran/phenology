#' ind_long_lat is used to manage ncdf information
#' @title Return or the index in ncdf object from lat/longitude or inverse
#' @author Marc Girondot
#' @return Or the index in ncdf object from lat/longitude or inverse
#' @param ncdf a ncdf object
#' @param long	longitude in decimal format
#' @param lat	latitude in decimal format
#' @param indice.long	Index of longitude
#' @param indice.lat	Index of latitude
#' @description Return or the index in ncdf object from lat/longitude or inverse
#' @examples
#' \dontrun{
#' url <- "ftp://ftp.cdc.noaa.gov/Datasets/noaa.oisst.v2.highres/"
#' url <- paste0(url, "sst.day.mean.2012.v2.nc")
#' dest <- paste(Sys.getenv("HOME"), "/sst.day.mean.2012.v2.nc", sep="")
#' download.file(url, dest)
#' library("ncdf4")
#' dta2012<-nc_open(paste(Sys.getenv("HOME"), "/sst.day.mean.2012.v2.nc", sep=""))
#' indices <- ind_long_lat(ncdf=dta2012, lat=5.89, long=-20.56)
#' coordinates <- ind_long_lat(ncdf=dta2012, indice.lat=20, indice.long=30)
#' }
#' @export


ind_long_lat<-function(ncdf=stop("The ncdf data must be supplied"), 
                       long=NA, lat=NA, indice.long=NA, indice.lat=NA) {
  if (!is.na(long) & !is.na(lat)) {
maxindicelg <- ncdf$dim$lon$len
maxindicelt <- ncdf$dim$lat$len
if (long==0) long <- 0.001
if (long==360) long <- 359.999
if (lat==90) lat <- 89.999
long<-long%%360

maxlg <- ncdf$dim$lon$vals[maxindicelg]
minlg <- ncdf$dim$lon$vals[1]
alg<-(maxindicelg-1)/(maxlg-minlg)
blg<-1-alg*minlg
lat<-((lat+90)%%180)-90
maxlt <- ncdf$dim$lat$vals[maxindicelt]
minlt <- ncdf$dim$lat$vals[1]
alt<-(maxindicelt-1)/(maxlt-minlt)
blt<-1-alt*minlt
return(c(indice.long=round(alg*long+blg), indice.lat=round(alt*lat+blt)))
} else {
  if (!is.na(indice.long) & !is.na(indice.lat)) {
    # Je fournis les indices et je calcule les coordonnées
    maxindicelg <- ncdf$dim$lon$len
    maxindicelt <- ncdf$dim$lat$len
    maxlg <- ncdf$dim$lon$vals[maxindicelg]
    minlg <- ncdf$dim$lon$vals[1]
    maxlt <- ncdf$dim$lat$vals[maxindicelt]
    minlt <- ncdf$dim$lat$vals[1]
    a_lg <- (maxlg-minlg)/(maxindicelg+1)
    b_lg <- minlg - a_lg
    a_lt <- (maxlt-minlt)/(maxindicelt+1)
    b_lt <- minlt - a_lt
    return(c(long=a_lg*indice.long+b_lg, lat=a_lt*indice.lat+b_lt))
    
  } else {
    warning("Check the parameters")
  }
}
}