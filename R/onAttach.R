.onAttach <- function(libname, pkgname) {
    actual <- utils::packageDescription(pkgname)[["Version"]]

    previouswarn <- getOption("warn")
    options(warn=2)
# put here your package CRAN description
    cranpage <- paste("http://cran.r-project.org/web/packages/", pkgname, "/index.html", sep="")
    description <- try(
        utils::read.delim(cranpage, header=FALSE, stringsAsFactors=FALSE),
    silent=TRUE)

    packageStartupMessage(paste("Welcome in package", pkgname, "!"))

    if (!inherits(description, "try-error")) {
        recent <- sub("<td>", "", description[14, 1])
        recent <- sub("</td></tr>", "", recent)

            if (as.numeric(actual) < as.numeric(recent)) {
                m <- paste("Your version is ", actual, ". Most recent is ", recent, sep="")
	      packageStartupMessage(m)
                packageStartupMessage("Use update.packages() to update...")
            }

    } else {
                m <- paste("Your version is ", actual)
                m <- paste(m, ". No internet connection is available to check for update.")
                packageStartupMessage(m)
  
    }

    options(warn=previouswarn)

}