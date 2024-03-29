#' Tagloss_format formats a CMR dataset into a file that Tagloss_L can use.
#' @title Format a CMR dataset into a file that Tagloss_L can use.
#' @author Marc Girondot
#' @return Return the maximum number of days an individual has been observed in a dataset.
#' @param data CMR file
#' @param model Can be "21" or "LR"
#' @param column.Date Name of the column with Date
#' @param column.ID Name of the column with ID
#' @param column.left Name of the column with left tag
#' @param column.right Name of the column with right tag
#' @param keeponly2 If TRUE, keep only individuals tagged with 2 tags
#' @param progressbar Is a progressbar been shown?
#' @description This function formats a CMR dataset to a file that Tagloss_L can use.\cr
#' The format of data is a data.frame with 4 columns:\cr
#' ID is the column with the permanent identification code\cr
#' L is the column with the non-permanent code located at left\cr
#' R is the column with the non-permanent code located at right\cr
#' Date is the column with the date of observation\cr
#' Note that R and L columns can be exchanged if 21 model is used.
#' @family Model of Tag-loss
#' @examples
#' \dontrun{
#' library(phenology)
#' # Example
#' head(outLR)
#' data_f_21 <- Tagloss_format(outLR, model="21")
#' data_f_LR <- Tagloss_format(outLR, model="LR")
#' }
#' @export

Tagloss_format <- function(data, model="21", progressbar=TRUE, 
                           column.Date="Date", column.ID="ID", 
                           column.left="L", column.right="R", 
                           keeponly2=TRUE) {
  
  # model="21"; progressbar=TRUE 
  # column.Date="Date"; column.ID="ID" 
  # column.left="L"; column.right="R"
  
  # ID
  ID <- unique(data[, column.ID])
  
  model <- match.arg(arg=model, choices = c("21", "LR"))
  
  if (progressbar) {
    if (is.element('progress', installed.packages()[,1])) {
      # library("progress")
      pb <- getFromNamespace("progress_bar", ns="progress")$new(
        format = "  completion [:bar] :percent eta: :eta",
        total = length(ID), clear = FALSE)
      libp <- TRUE
    } else {
      libp <- FALSE
      pb <- txtProgressBar(min=0, max=length(ID), style=3)
    }
  }
  cptlv <- 0
  
  # options(warn = 2)
  
  if (model=="21") {
    # model="21"
    out <- data.frame(ID=character(), N22=numeric(), N21=numeric(), N11=numeric(), 
                      N10=numeric(), N20=numeric())
    for (individu in ID) {
      
      if (progressbar) {
        if (libp) {
          pb$tick() 
        } else {
          cptlv <- cptlv + 1
          setTxtProgressBar(pb, cptlv)
        }
      }
      
      N22 <- N21 <- N11 <- N10 <- N20 <- NA
      data_i <- subset(data, subset=(data[, column.ID]==individu))
      data_i <- data_i[order(data_i[, column.Date]), ]
      # Test supplémentaire: si même date, soit idem alors je supprime, sinon je prends
      # celui avec deux bagues
      
      if (nrow(data_i) > 1) {
        date_pre <- data_i[1, column.Date]
        for (i in 2:nrow(data_i)) {
          if (data_i[i, column.Date] == date_pre) {
            BL1 <- data_i[i-1, column.left]
            BR1 <- data_i[i-1, column.right]
            if (is.na(BL1)) BL1 <- ""
            if (is.na(BR1)) BR1 <- ""
            BL2 <- data_i[i, column.left]
            BR2 <- data_i[i, column.right]
            if (is.na(BL2)) BL2 <- ""
            if (is.na(BR2)) BR2 <- ""
            if ((BL1 != "") + (BR1 != "") == (BL2 != "") + (BR2 != "")) {
              # Les mêmes
              data_i[i-1, column.Date] <- NA
            } else {
              if ((BL1 != "") + (BR1 != "") > (BL2 != "") + (BR2 != "")) {
                # Date i-1 a 2 bagues et date i moins, je garde i
                data_i[i-1, column.Date] <- NA
              } else {
                data_i[i, column.Date] <- NA
              }
            }
          }
          if (!is.na(data_i[i, column.Date])) date_pre <- data_i[i, column.Date]
        }
        data_i <- data_i[!is.na(data_i[, column.Date]), ]
      }
      
      BL <- data_i[1, column.left]
      BR <- data_i[1, column.right]
      # Je dois virer quand une nouvelle bague apparaît
      if (is.na(BL)) BL <- ""
      if (is.na(BR)) BR <- ""
      
      keep <- TRUE
      if ((BR == "") & (BL == "")) keep <- FALSE
      if (((BR == "") | (BL == "")) & keeponly2) keep <- FALSE
      
      if ((nrow(data_i) > 1) & keep) {
        
        data_i[, column.left] <- ifelse(data_i[, column.left] != BL, "", data_i[, column.left])
        data_i[, column.right] <- ifelse(data_i[, column.right] != BR, "", data_i[, column.right])
        data_i[, column.left] <- ifelse(is.na(data_i[, column.left]), "", data_i[, column.left])
        data_i[, column.right] <- ifelse(is.na(data_i[, column.right]), "", data_i[, column.right])
        
        data_i[1:max(which(data_i[, column.left] != "")), column.left] <- data_i[1, column.left]
        data_i[1:max(which(data_i[, column.right] != "")), column.right] <- data_i[1, column.right]
        
        data_i <- cbind(data_i, Ordinal=as.numeric(data_i[, column.Date]-data_i[1, column.Date]))
        N22 <- max(subset(data_i, subset=(data_i[, column.left] != "") & (data_i[, column.right] != ""))$Ordinal)
        data_i <- subset(data_i, subset=(data_i[, column.left] == "") | (data_i[, column.right] == ""))
        if (nrow(data_i) != 0) {
          if ((data_i[1, column.right] == "") & (data_i[1, column.left] == "")) {
            N20 <- data_i[1, "Ordinal"] - N22
          } else {
            if (data_i[1, column.right] != "") {
              N21 <- data_i[1, "Ordinal"] - N22
              data_i <- data_i[-1, ]
              if (nrow(data_i) != 0) {
                if (nrow(subset(data_i, subset=(data_i[, column.right] != ""))) == 0) N11 <- NA else N11 <- max(subset(data_i, subset=(data_i[, column.right] != ""))$Ordinal) - N22 - N21
                data_i <- subset(data_i, subset=(data_i[, column.right] == ""))
                if (nrow(data_i) != 0) {
                  N10 <- data_i[1, "Ordinal"] - N22 - N21 - N11
                }
              }
            } else {
              N21 <- data_i[1, "Ordinal"] - N22
              data_i <- data_i[-1, ]
              if (nrow(data_i) != 0) {
                if (nrow(subset(data_i, subset=(data_i[, column.left] != ""))) == 0) N11 <- NA else N11 <- max(subset(data_i, subset=(data_i[, column.left] != ""))$Ordinal) - N22 - N21
                data_i <- subset(data_i, subset=(data_i[, column.left] == ""))
                if (nrow(data_i) != 0) {
                  N10 <- data_i[1, "Ordinal"] - N22 - N21 - N11
                }
              }
            }
          }
        }
        out <- rbind(out, data.frame(ID=individu, N22=N22, N21=N21, N11=N11, N10=N10, N20=N20))
      }
    }
    # fin modèle 21
  } else {
    # model="LR"
    out <- data.frame(ID=character(), NLR_LR=numeric(), NLR_0R=numeric(), NLR_L0=numeric(), 
                      N0R_0R=numeric(), NL0_L0=numeric(), N0R_00=numeric(), NL0_00=numeric(), NLR_00=numeric())
    
    #  options(warn=2)
    for (individu in ID) {
      
      if (progressbar) {
        if (libp) {
          pb$tick() 
        } else {
          cptlv <- cptlv + 1
          setTxtProgressBar(pb, cptlv)
        }
      }
      
      # print(individu)
      NLR_LR <- NLR_0R <- NLR_L0 <- N0R_0R <- NL0_L0 <- N0R_00 <- NL0_00 <- NLR_00 <- NA
      data_i <- subset(data, subset=(data[, column.ID] == individu))
      data_i <- data_i[order(data_i[, column.Date]), ]
      
      # Test supplémentaire: si même date, soit idem alors je supprime, sinon je prends
      # celui avec deux bagues
      
      if (nrow(data_i) > 1) {
        date_pre <- data_i[1, column.Date]
        for (i in 2:nrow(data_i)) {
          if (data_i[i, column.Date] == date_pre) {
            BL1 <- data_i[i-1, column.left]
            BR1 <- data_i[i-1, column.right]
            if (is.na(BL1)) BL1 <- ""
            if (is.na(BR1)) BR1 <- ""
            BL2 <- data_i[i, column.left]
            BR2 <- data_i[i, column.right]
            if (is.na(BL2)) BL2 <- ""
            if (is.na(BR2)) BR2 <- ""
            if ((BL1 != "") + (BR1 != "") == (BL2 != "") + (BR2 != "")) {
              # Les mêmes
              data_i[i-1, column.Date] <- NA
            } else {
              if ((BL1 != "") + (BR1 != "") > (BL2 != "") + (BR2 != "")) {
                # Date i-1 a 2 bagues et date i moins, je garde i
                data_i[i-1, column.Date] <- NA
              } else {
                data_i[i, column.Date] <- NA
              }
            }
          }
          if (!is.na(data_i[i, column.Date])) date_pre <- data_i[i, column.Date]
        }
        data_i <- data_i[!is.na(data_i[, column.Date]), ]
      }
      
      
      # Je dois virer quand une nouvelle bague apparaît
      BL <- data_i[1, column.left]
      BR <- data_i[1, column.right]
      if (is.na(BL)) BL <- ""
      if (is.na(BR)) BR <- ""
      
      keep <- TRUE
      if ((BR == "") & (BL == "")) keep <- FALSE
      if (((BR == "") | (BL == "")) & keeponly2) keep <- FALSE
      
      if ((nrow(data_i) > 1) & keep) {
        
        data_i[, column.left] <- ifelse(data_i[, column.left] != BL, "", data_i[, column.left])
        data_i[, column.right] <- ifelse(data_i[, column.right] != BR, "", data_i[, column.right])
        data_i[, column.left] <- ifelse(is.na(data_i[, column.left]), "", data_i[, column.left])
        data_i[, column.right] <- ifelse(is.na(data_i[, column.right]), "", data_i[, column.right])
        
        
        data_i <- cbind(data_i, Ordinal=as.numeric(data_i[, column.Date]-data_i[1, column.Date]))
        LR.matrix <- matrix(rep(NA, 2*data_i[nrow(data_i), "Ordinal"]+2), ncol=2)
        LR.matrix[data_i[which(data_i[, column.left] != ""), "Ordinal"]+1, 1] <- 1
        LR.matrix[data_i[which(data_i[, column.right] != ""), "Ordinal"]+1, 2] <- 1
        LR.matrix[data_i[which(data_i[, column.left] == ""), "Ordinal"]+1, 1] <- 0
        LR.matrix[data_i[which(data_i[, column.right] == ""), "Ordinal"]+1, 2] <- 0
        dl <- paste0(LR.matrix[, 1], LR.matrix[, 2])
        nb0R <- nbL0 <- NA
        if (identical(which(dl=="11"), integer())) {
          nbLR <- 0
        } else {
          nbLR <- max(which(dl=="11"))
        }
        
        nbL0 <- which(dl=="10")
        if (!identical(nbL0, integer(0))) nbL0 <- na.omit(ifelse(nbL0<nbLR, NA, nbL0))
        nb0R <- which(dl=="01")
        if (!identical(nb0R, integer(0))) nb0R <- na.omit(ifelse(nb0R<nbLR, NA, nb0R))
        nb00 <- which(dl=="00")
        if (!identical(nb00, integer(0))) nb00 <- na.omit(ifelse(nb00<max(nbLR, nb0R, nbL0, na.rm = TRUE), NA, nb00))
        attributes(nbL0) <- NULL
        attributes(nb0R) <- NULL
        attributes(nb00) <- NULL
        
        NLR_LR_calcul <- nbLR - 1
        if (nbLR != 1) NLR_LR <- NLR_LR_calcul
        if (!identical(nbL0, integer(0)) & !identical(nbL0, logical(0))) NLR_L0 <- min(nbL0) - NLR_LR_calcul -1
        if (!identical(nb0R, integer(0)) & !identical(nb0R, logical(0))) NLR_0R <- min(nb0R) - NLR_LR_calcul -1
        if (!identical(nbL0, integer(0)) & !identical(nbL0, logical(0))) if (length(nbL0) != 1) NL0_L0 <- max(nbL0) - NLR_LR_calcul - NLR_L0 -1
        if (!identical(nb0R, integer(0)) & !identical(nb0R, logical(0))) if (length(nb0R) != 1) N0R_0R <- max(nb0R) - NLR_LR_calcul - NLR_0R -1
        if (!identical(nb00, integer(0)) & !identical(nb00, logical(0))) NL0_00 <- min(nb00) - NLR_LR_calcul - NLR_L0 - ifelse(is.na(NL0_L0), 0, NL0_L0) -1
        if (!identical(nb00, integer(0)) & !identical(nb00, logical(0))) N0R_00 <- min(nb00) - NLR_LR_calcul - NLR_0R - ifelse(is.na(N0R_0R), 0, N0R_0R) -1
        if (!identical(nb00, integer(0))  & !identical(nb00, logical(0)) & (is.na(NL0_00) & is.na(N0R_00))) NLR_00 <- min(nb00) - NLR_LR_calcul -1
        
        out <- rbind(out, data.frame(ID=individu, NLR_LR=NLR_LR, NLR_0R=NLR_0R, NLR_L0=NLR_L0, N0R_0R=N0R_0R, 
                                     NL0_L0=NL0_L0, N0R_00=N0R_00, NL0_00=NL0_00, NLR_00=NLR_00))
        
      }
    }
  }
  out <- addS3Class(out, "TaglossData")
  return(out)
}
