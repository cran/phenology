#' Tagloss_daymax returns the maximum number of days an individual has been observed in a dataset.
#' @title Return the maximum number of days an individual has been observed in a dataset.
#' @author Marc Girondot \email{marc.girondot@@gmail.com}
#' @return Return the maximum number of days an individual has been observed in a dataset.
#' @param individuals Set of indivuals
#' @param what By default is max, but can be min, mean or all
#' @description This function must be used to get the value of mx in Tagloss_L.\cr
#' @family Model of Tag-loss
#' @examples
#' \dontrun{
#' library(phenology)
#' # Example
#' data_f_21 <- Tagloss_format(outLR, model="21")
#' daymax(data_f_21)
#' }
#' @export


Tagloss_daymax <- function (individuals, what = "max") {
  # Je pense que ça n'est plus nécessaire, j'essaye de le retirer
  if (FALSE)
  if (inherits(individuals, "TaglossData")) 
    class(individuals) <- unique(append("data.frame", class(individuals)))
  
  if ((!inherits(individuals,  "matrix")) & (!inherits(individuals, "data.frame"))) 
    individuals <- matrix(data = individuals, nrow = 1, 
                          dimnames = list(NULL,  names(individuals)))
  cn <- colnames(individuals)
  totnames <- c("NLR_LR", "NLR_L0", "NLR_0R", "NL0_L0", "N0R_0R", 
                "NL0_00", "N0R_00", "NLR_00", "N22", "N21", "N11", "N10", 
                "N20")
  newnames <- totnames[!(totnames %in% cn)]
  if (length(newnames) != 0) {
    individuals <- cbind(individuals, matrix(data = rep(NA, 
                                                        nrow(individuals) * length(newnames)), nrow = nrow(individuals), 
                                             dimnames = list(NULL, newnames)))
  }
  mx <- 0
  mn <- +Inf
  me <- NULL
  for (individu in 1:nrow(individuals)) {
    NLR_LR <- individuals[individu, "NLR_LR"]
    NLR_L0 <- individuals[individu, "NLR_L0"]
    NLR_0R <- individuals[individu, "NLR_0R"]
    NL0_L0 <- individuals[individu, "NL0_L0"]
    N0R_0R <- individuals[individu, "N0R_0R"]
    NL0_00 <- individuals[individu, "NL0_00"]
    N0R_00 <- individuals[individu, "N0R_00"]
    NLR_00 <- individuals[individu, "NLR_00"]
    N22 <- individuals[individu, "N22"]
    N21 <- individuals[individu, "N21"]
    N11 <- individuals[individu, "N11"]
    N20 <- individuals[individu, "N20"]
    N10 <- individuals[individu, "N10"]
    NLR_LR <- ifelse(is.na(NLR_LR), 0, NLR_LR)
    NLR_L0 <- ifelse(is.na(NLR_L0), 0, NLR_L0)
    NLR_0R <- ifelse(is.na(NLR_0R), 0, NLR_0R)
    NL0_L0 <- ifelse(is.na(NL0_L0), 0, NL0_L0)
    N0R_0R <- ifelse(is.na(N0R_0R), 0, N0R_0R)
    NL0_00 <- ifelse(is.na(NL0_00), 0, NL0_00)
    N0R_00 <- ifelse(is.na(N0R_00), 0, N0R_00)
    NLR_00 <- ifelse(is.na(NLR_00), 0, NLR_00)
    N22 <- ifelse(is.na(N22), 0, N22)
    N21 <- ifelse(is.na(N21), 0, N21)
    N11 <- ifelse(is.na(N11), 0, N11)
    N20 <- ifelse(is.na(N20), 0, N20)
    N10 <- ifelse(is.na(N10), 0, N10)
    mx <- max(c(mx, NLR_LR, NLR_LR + NLR_L0 + NL0_L0, NLR_LR + 
                  NLR_0R + N0R_0R, NLR_LR + NLR_L0 + NL0_L0 + NL0_00, 
                NLR_LR + NLR_0R + N0R_0R + N0R_00, NLR_00, N22, N22 + 
                  N21 + N11 + N10, N20), na.rm = TRUE)
    mn <- min(c(mn, NLR_LR, NLR_LR + NLR_L0 + NL0_L0, NLR_LR + 
                  NLR_0R + N0R_0R, NLR_LR + NLR_L0 + NL0_L0 + NL0_00, 
                NLR_LR + NLR_0R + N0R_0R + N0R_00, NLR_00, N22, N22 + 
                  N21 + N11 + N10, N20), na.rm = TRUE)
    me <- c(me, max(c(NLR_LR, NLR_LR + NLR_L0 + NL0_L0, NLR_LR + 
                        NLR_0R + N0R_0R, NLR_LR + NLR_L0 + NL0_L0 + NL0_00, 
                      NLR_LR + NLR_0R + N0R_0R + N0R_00, NLR_00, N22, N22 + 
                        N21 + N11 + N10, N20), na.rm = TRUE))
  }
  out <- mx
  if (what == "min") 
    out <- mn
  if (what == "mean") 
    out <- mean(me)
  if (what == "all") 
    out <- me
  return(out)
}
