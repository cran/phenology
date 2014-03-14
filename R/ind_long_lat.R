#' ind_long_lat is used to manage ncdf information
#' @title Return or the index in ncdf object from lat/longitude or inverse
#' @author Marc Girondot
#' @return Or the index in ncdf object from lat/longitude or inverse
#' @param ncdf an object read from package ncdf4, RnetCDF or ncdf
#' @param long longitude in decimal format
#' @param lat	latitude in decimal format
#' @param indice.long	Index of longitude
#' @param indice.lat	Index of latitude
#' @param name.lon Name of argument for longitude, default is lon
#' @param name.lat Name of argument for latitude, default is lat
#' @description Return or the index in ncdf object from lat/longitude or inverse
#' @examples
#' \dontrun{
#' url <- "ftp://ftp.cdc.noaa.gov/Datasets/noaa.oisst.v2.highres/"
#' url <- paste0(url, "sst.day.mean.2012.v2.nc")
#' dest <- paste(Sys.getenv("HOME"), "/sst.day.mean.2012.v2.nc", sep="")
#' download.file(url, dest)
#' library("ncdf4")
#' dta2012 <- nc_open(dest)
#' indices <- ind_long_lat(ncdf=dta2012, lat=5.89, long=-20.56)
#' coordinates <- ind_long_lat(ncdf=dta2012, indice.lat=20, indice.long=30)
#' library("RNetCDF")
#' dta2012 <- open.nc(dest)
#' indices <- ind_long_lat(ncdf=dta2012, lat=5.89, long=-20.56)
#' coordinates <- ind_long_lat(ncdf=dta2012, indice.lat=20, indice.long=30)
#' library("ncdf")
#' dta2012 <- open.ncdf(dest)
#' indices <- ind_long_lat(ncdf=dta2012, lat=5.89, long=-20.56)
#' coordinates <- ind_long_lat(ncdf=dta2012, indice.lat=20, indice.long=30)
#' }
#' @export


ind_long_lat<-function(ncdf=stop("The ncdf data must be supplied"), 
                       long=NA, lat=NA, indice.long=NA, indice.lat=NA,
                       name.lon="lon", name.lat="lat") {
  maxindicelt <- NULL  
  if (class(ncdf)=="ncdf4") {
    maxindicelg <- ncdf$dim[[name.lon]]$len
    maxindicelt <- ncdf$dim[[name.lat]]$len
    maxlg <- ncdf$dim[[name.lon]]$vals[maxindicelg]
    minlg <- ncdf$dim[[name.lon]]$vals[1]
    maxlt <- ncdf$dim[[name.lat]]$vals[maxindicelt]
    minlt <- ncdf$dim[[name.lat]]$vals[1]
  }
  
  if (class(ncdf)=="ncdf") {
    maxindicelt <- ncdf$dim[[name.lat]]$len
    maxlt <- ncdf$dim[[name.lat]]$vals[maxindicelt]
    minlt <- ncdf$dim[[name.lat]]$vals[1]
    maxindicelg <- ncdf$dim[[name.lon]]$len   
    maxlg <- ncdf$dim[[name.lon]]$vals[maxindicelg]
    minlg <- ncdf$dim[[name.lon]]$vals[1]
  }
  
  if (class(ncdf)=="NetCDF") {    
    # latitude: range and length
    range.1 <- att.get.nc(ncdf, name.lat, "actual_range")
    maxindicelt <- dim.inq.nc(ncdf, name.lat)$length
    maxlt <- range.1[2]
    minlt <- range.1[1]
    
    # longitude: range and length
    range.2 <- att.get.nc(ncdf, name.lon, "actual_range")
    maxindicelg <- dim.inq.nc(ncdf, name.lon)$length
    maxlg <- range.2[2]
    minlg <- range.2[1]   
  }
  
  if (is.null(maxindicelt) | is.null(maxindicelg)) {
    warning("Check the ncdf data; it is not recognized")
    return(invisible())
  }
  
  if (!is.na(long) & !is.na(lat)) {

if (long==0) long <- 0.0001
if (long==360) long <- 359.9999
if (lat==90) lat <- 89.9999

long<-long%%360
alg<-(maxindicelg-1)/(maxlg-minlg)
blg<-1-alg*minlg

lat<-((lat+90)%%180)-90
alt<-(maxindicelt-1)/(maxlt-minlt)
blt<-1-alt*minlt

return(c(indice.long=round(alg*long+blg), indice.lat=round(alt*lat+blt)))
} else {
  if (!is.na(indice.long) & !is.na(indice.lat)) {
    # Je fournis les indices et je calcule les coordonnées
    return(c(long=as.numeric(indice.long)*360/maxindicelg-minlg, 
             lat=(as.numeric(indice.lat)-1)*180/maxindicelt+minlt))
    
  } else {
    warning("Check the parameters")
    return(invisible())    
  }
}
}