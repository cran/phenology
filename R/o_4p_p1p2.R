#' Model of tagloss based on Rivalan data
#' @title Model of tagloss based on Rivalan data
#' @author Marc Girondot
#' @docType data
#' @name o_4p_p1p2
#' @description Model of tagloss based on Rivalan data
#' @keywords datasets
#' @usage o_4p_p1p2
#' @family Model of Tag-loss
#' @references Rivalan, P., Godfrey, M.H., Prévot-Julliard, A.-C., Girondot, M., 2005. Maximum likelihood estimates of tag loss in leatherback sea turtles. Journal of Wildlife Management 69, 540-548.
#' @examples
#' \dontrun{
#' library(phenology)
#' # Example
#' data_f_21 <- Tagloss_format(outLR, model="21")
#' 
#' # Same data fitted with new model
#' par <- c(D1_1 = 100.15324837975547, A_1 = 5.9576927964120188, 
#'          B_1 = 8.769924225871069, B_2 = 8.2353860179664125)
#' pfixed <- c(D2D1_1 = 2568, D3D2_1 = 2568, D2D1_2 = 2568, D3D2_2 = 2568)
#' o_4p_p1p2 <- Tagloss_fit(data=data_f_21, fitted.parameters = par, 
#'                          fixed.parameters = pfixed, 
#'                          model_before = "par['C_1']=par['B_1'];
#'                          par['A_2']=par['A_1'];
#'                          par['C_2']=par['B_2'];
#'                          par['D1_2']=par['D1_1']", hessian=TRUE)
#' data(o_4p_p1p2)
#' plot(o_4p_p1p2, model="1", col="red")
#' plot(o_4p_p1p2, model="2", col="blue", add=TRUE)
#' legend("topright", legend=c("2->1", "1->0"), lty=1, col=c("blue", "red"))
#' plot(o_4p_p1p2, model="Cumul")
#' }
NULL
