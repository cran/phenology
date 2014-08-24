#' convert.ts Convert one Date-Time from one timezone to another
#' @title Convert one Date-Time from one timezone to another
#' @author Marc Girondot
#' @return A POSIXlt date converted
#' @param x The date-time in POSIXlt or POSIXct format
#' @param tz The timezone
#' @description Convert one Date-Time from one timezone to another.
#' Available timezones can be shown using OlsonNames()
#' @examples
#' d <- as.POSIXlt("2010-01-01 17:34:20", tz="UTC")
#' convert.ts(d, tz="America/Guatemala")
#' @export


convert.ts <- function(x, tz=Sys.timezone()) {
return(as.POSIXct(format(as.POSIXct(x), tz=tz, usetz=TRUE), tz=tz))
}
