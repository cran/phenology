#' Tagloss_L returns the -log likelihood of a set of individuals under a model of tagloss.
#' @title Return the -log likelihood of a set of individuals under a model of tagloss.
#' @author Marc Girondot
#' @return Return the -log likelihood of a set of individuals
#' @param individuals Set of indivuals
#' @param par Set of parameters
#' @param days.maximum Maximum number of days. Can be determined using Tagloss_daymax()
#' @param fixed.parameters Set of fixed parameters
#' @param model_before Transformation of parameters before to use Tagloss_model()
#' @param model_after Transformation of parameters after to use Tagloss_model()
#' @param names.par Name of parameters. Normally unused. 
#' @param groups Number of groups for parallel computing
#' @param cores Number of cores to use for parallel computing
#' @param progressbar Is shown a progressbar?
#' @description This function must be used within optim().\cr
#'   model_before is applied to the par parameter.\cr
#'   model_after is applied after par is separated in p1, p2, pL1, pL2, pR1 and pR2 parameters.\cr
#' progressbar is set to FALSE if cores is different from 1.
#' @family Model of Tag-loss
#' @examples
#' \dontrun{
#' library(phenology)
#' 
#' # Example with 21 format of data
#' 
#' data_f_21 <- Tagloss_format(outLR, model="21")
#' par <- structure(c(49.5658922243074, 808.136085362158, 106.283783786853, 
#' 5.22150592456511, 8.00608716525864, 8.32718202233396, 150.612916258503, 
#' 715.865805125223, 2242.06574225966, 119.212383120678, 10.1860735529433, 
#' 7.14231725937626), .Names = c("D1_2", "D2D1_2", "D3D2_2", "A_2", 
#' "B_2", "C_2", "D1_1", "D2D1_1", "D3D2_1", "A_1", "B_1", "C_1"))
#' pfixed <- NULL
#' # All the data are analyzed; the N20 are very long to compute
#' Tagloss_L(individuals=data_f_21, par=par, days.maximum=Tagloss_daymax(data_f_21), 
#'           fixed.parameters=pfixed, cores=1, progressbar=TRUE)
#' # Without the N20 the computing is much faster
#' data_f_21_fast <- subset(data_f_21, subset=(is.na(data_f_21$N20)))
#' Tagloss_L(individuals=data_f_21_fast, par=par, days.maximum=Tagloss_daymax(data_f_21_fast), 
#'           fixed.par=pfixed, cores=1, progressbar=TRUE)
#' o <- Tagloss_fit(data=data_f_21_fast, fitted.parameters=par)
#' # Here it is the result of the previous function
#' o <- structure(list(par = structure(c(49.5658922243074, 808.136085362158, 
#' 106.283783786853, 5.22150592456511, 8.00608716525864, 8.32718202233396, 
#' 150.612916258503, 715.865805125223, 2242.06574225966, 119.212383120678, 
#' 10.1860735529433, 7.14231725937626), .Names = c("D1_2", "D2D1_2", 
#' "D3D2_2", "A_2", "B_2", "C_2", "D1_1", "D2D1_1", "D3D2_1", "A_1", 
#' "B_1", "C_1")), value = 5841.93084262461, counts = structure(c(1093L, 
#' NA), .Names = c("function", "gradient")), convergence = 0L, message = NULL, 
#'     hessian = structure(c(0.0469808583147824, 0.000133240973809734, 
#'     6.68478605803102e-05, -2.53581288234273, -1.25931342154217, 
#'     -0.124650568977813, -2.46700437855907e-05, -1.11413100967184e-05, 
#'     -3.18323145620525e-06, 0, -0.0182945996130002, -0.00510601694259094, 
#'     0.000133240973809734, 1.45519152283669e-05, 7.50333128962666e-06, 
#'     -0.00452587300969753, -0.0191316757991444, -0.0255117811320815, 
#'     -1.13686837721616e-06, -1.36424205265939e-06, -2.27373675443232e-07, 
#'     0, 0.000335830918629654, -0.000448608261649497, 6.68478605803102e-05, 
#'     7.50333128962666e-06, 4.32009983342141e-06, -0.00226373231271282, 
#'     -0.00954059942159802, -0.0127809016703395, -4.54747350886464e-07, 
#'     -4.54747350886464e-07, -2.27373675443232e-07, 0, 0.000176896719494835, 
#'     -0.000224190443987027, -2.53581288234273, -0.00452587300969753, 
#'     -0.00226373231271282, 223.422489398217, 41.4073996353181, 
#'     3.77875949197914, 0.000986460690910462, 0.000398813426727429, 
#'     0.000117665877041873, 0, 0.727547330825473, 0.194675862985605, 
#'     -1.25931342154217, -0.0191316757991444, -0.00954059942159802, 
#'     41.4073996353181, 189.534394394286, 28.3386068531399, 0.00216437001654413, 
#'     0.00241834641201422, 0.000652562448522076, 0, 0.841939595375152, 
#'     1.0472297162778, -0.124650568977813, -0.0255117811320815, 
#'     -0.0127809016703395, 3.77875949197914, 28.3386068531399, 
#'     70.250493081403, -0.00022441781766247, -0.000161662683240138, 
#'     0.000257614374277182, 0, -0.578908839088399, 1.08917492980254, 
#'     -2.46700437855907e-05, -1.13686837721616e-06, -4.54747350886464e-07, 
#'     0.000986460690910462, 0.00216437001654413, -0.00022441781766247, 
#'     0.000148247636388987, 0.000145519152283669, 3.97903932025656e-05, 
#'     0, 0.0156976511789253, 0.0678746800986119, -1.11413100967184e-05, 
#'     -1.36424205265939e-06, -4.54747350886464e-07, 0.000398813426727429, 
#'     0.00241834641201422, -0.000161662683240138, 0.000145519152283669, 
#'     0.000145519152283669, 3.9676706364844e-05, 0, 0.0138438736030366, 
#'     0.0678776359563926, -3.18323145620525e-06, -2.27373675443232e-07, 
#'     -2.27373675443232e-07, 0.000117665877041873, 0.000652562448522076, 
#'     0.000257614374277182, 3.97903932025656e-05, 3.9676706364844e-05, 
#'     1.77351466845721e-05, 0, 0.00317095327773131, 0.0316927071253303, 
#'     0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.0182945996130002, 
#'     0.000335830918629654, 0.000176896719494835, 0.727547330825473, 
#'     0.841939595375152, -0.578908839088399, 0.0156976511789253, 
#'     0.0138438736030366, 0.00317095327773131, 0, 8.85630879565724, 
#'     4.44044781033881, -0.00510601694259094, -0.000448608261649497, 
#'     -0.000224190443987027, 0.194675862985605, 1.0472297162778, 
#'     1.08917492980254, 0.0678746800986119, 0.0678776359563926, 
#'     0.0316927071253303, 0, 4.44044781033881, 88.8524673428037
#'     ), .Dim = c(12L, 12L), .Dimnames = list(c("D1_2", "D2D1_2", 
#'     "D3D2_2", "A_2", "B_2", "C_2", "D1_1", "D2D1_1", "D3D2_1", 
#'     "A_1", "B_1", "C_1"), c("D1_2", "D2D1_2", "D3D2_2", "A_2", 
#'     "B_2", "C_2", "D1_1", "D2D1_1", "D3D2_1", "A_1", "B_1", "C_1"
#'     )))), .Names = c("par", "value", "counts", "convergence", 
#' "message", "hessian"), class = "Tagloss")
#' par(mar=c(4, 4, 1, 1))
#' plot(o, t=1:3000, model="2", scale=1000, ylim=c(0, 3), 
#'             col="red")
#' plot(o, t=1500:3000, model="1", scale=1000, 
#'             add=TRUE)
#' legend("topright", legend=c("2 -> 1", "1 -> 0"), col=c("red", "black"), lty=1)
#' 
#' plot(o, t=1:300, model="2", scale=1000, ylim=c(0, 3), 
#'             col="red", hessian=o$hessian)
#' plot(o, t=1:300, model="1", scale=1000, 
#'             add=TRUE, hessian=o$hessian)
#' legend("topright", legend=c("2 -> 1", "1 -> 0"), col=c("red", "black"), lty=1)
#' 
#' ###### Example with fixed.parameters
#' 
#' data_f_21 <- Tagloss_format(outLR, model="21")
#' # Without the N20 the computing is much faster
#' data_f_21_fast <- subset(data_f_21, subset=(is.na(data_f_21$N20)))
#' par <- structure(c(49.5658922243074, 5.22150592456511, 8.00608716525864, 
#'                    50.612916258503, 6, 9), 
#'                 .Names = c("D1_2",  "A_2", "B_2", 
#'                            "D1_1",  "A_1", "B_1"))
#' pfixed <- c(D2D1_2=10000, D3D2_2=10000, C_2=0, D2D1_1=10000, D3D2_1=10000, C_1=0)
#' o <- Tagloss_fit(data=data_f_21_fast, fitted.parameters=par, fixed.parameters=pfixed)
#' # Here it is the result of the previous function
#' o <- structure(list(par = structure(c(55.2184044121564, 5.2630294044259, 
#' 8.13359029885985, 14269.9757684677, 21.8702023948044, 6.46586480967269
#' ), .Names = c("D1_2", "A_2", "B_2", "D1_1", "A_1", "B_1")), value = 5853.64634357369, 
#'     counts = structure(c(757L, NA), .Names = c("function", "gradient"
#'     )), convergence = 0L, message = NULL, hessian = structure(c(0.036636720324168, 
#'     -2.26385645873961, -1.2330608569755, -2.95585778076202e-06, 
#'     -2.27373675443232e-07, -0.0399197688238928, -2.26385645873961, 
#'     232.345637869003, 47.1904784262733, 0.000118689058581367, 
#'     7.50333128962666e-06, 1.69928603099834, -1.2330608569755, 
#'     47.1904784262733, 304.432723851278, 0.000196678229258396, 
#'     1.36424205265939e-06, 2.8553522497532, -2.95585778076202e-06, 
#'     0.000118689058581367, 0.000196678229258396, 4.54747350886464e-07, 
#'     0, 0.00741636085876962, -2.27373675443232e-07, 7.50333128962666e-06, 
#'     1.36424205265939e-06, 0, 4.00177668780088e-05, 8.79936123965308e-05, 
#'     -0.0399197688238928, 1.69928603099834, 2.8553522497532, 0.00741636085876962, 
#'     8.79936123965308e-05, 107.941018768543), .Dim = c(6L, 6L), .Dimnames = list(
#'         c("D1_2", "A_2", "B_2", "D1_1", "A_1", "B_1"), c("D1_2", 
#'         "A_2", "B_2", "D1_1", "A_1", "B_1")))), .Names = c("par", 
#' "value", "counts", "convergence", "message", "hessian"), class = "Tagloss")
#' par(mar=c(4, 4, 1, 1))
#' plot(o, t=1:3000, model="2", scale=1000, ylim=c(0, 3), 
#'             col="red")
#' plot(o, t=1500:3000, model="1", scale=1000, 
#'             add=TRUE)
#' legend("topright", legend=c("2 -> 1", "1 -> 0"), col=c("red", "black"), lty=1)
#' 
#' plot(o, t=1:300, model="2", scale=1000, ylim=c(0, 3), 
#'             col="red", hessian=o$hessian)
#' plot(o, t=1:300, model="1", scale=1000, 
#'             add=TRUE, hessian=o$hessian)
#' legend("topright", legend=c("2 -> 1", "1 -> 0"), col=c("red", "black"), lty=1)
#' 
#' ###### Example with delta
#' 
#' data_f_21 <- Tagloss_format(outLR, model="21")
#' # Without the N20 the computing is much faster
#' data_f_21_fast <- subset(data_f_21, subset=(is.na(data_f_21$N20)))
#' par <- structure(c(45.8764973711504, 5.22489974562498, 8.07602162728874, 
#' -0.865444694177429), .Names = c("D1_2", "A_2", "B_2", "delta"
#' ))
#' pfixed <- c(D2D1_2=10000, D3D2_2=10000, C_2=0)
#' o <- Tagloss_fit(data=data_f_21_fast, fitted.parameters=par, fixed.parameters=pfixed)
#' # Here it is the result of the previous function
#' o <- structure(list(par = structure(c(45.9035484983855, 5.22576211343279, 
#' 8.07585745169786, -0.865706100004634), .Names = c("D1_2", "A_2", 
#' "B_2", "delta")), value = 5913.716964613, counts = structure(c(91L, 
#' NA), .Names = c("function", "gradient")), convergence = 0L, message = NULL, 
#'     hessian = structure(c(0.0644593001197791, -2.88983483187621, 
#'     -1.49161280660337, -0.0875163550517755, -2.88983483187621, 
#'     221.02317802819, 45.3729608125286, 3.73816044429987, -1.49161280660337, 
#'     45.3729608125286, 440.129730122862, 30.4781699469459, -0.0875163550517755, 
#'     3.73816044429987, 30.4781699469459, 9.47964940678503), .Dim = c(4L, 
#'     4L), .Dimnames = list(c("D1_2", "A_2", "B_2", "delta"), c("D1_2", 
#'     "A_2", "B_2", "delta")))), .Names = c("par", "value", "counts", 
#' "convergence", "message", "hessian"), class = "Tagloss")
#' par(mar=c(4, 4, 1, 1))
#' plot(o, t=1:3000, model="2", scale=1000, ylim=c(0, 3), 
#'             col="red")
#' plot(o, t=1:3000, model="1", scale=1000, col="blue", 
#'             add=TRUE, hessian=o$hessian)
#' legend("topright", legend=c("2 -> 1", "1 -> 0"), col=c("red", "black"), lty=1)
#' 
#' ###### Example with model_after
#' data_f_LR <- Tagloss_format(outLR, model="LR")
#' par <- structure(c(72.0399239978454, 58.1034231071992, 645.068735669251, 
#'                    5.10791337470247, 3538.47220045768, 7.83358940767931), 
#'                 .Names = c("D1_L2", "D2D1_L2", "D3D2_L2", "A_L2", "B_L2", "C_L2"))
#' pfixed <- NULL
#' # A progress bar can be shown when one core is used
#' system.time(
#' print(Tagloss_L(individuals=data_f_LR, par=par, days.maximum=Tagloss_daymax(data_f_LR), 
#'           fixed.parameters=pfixed, cores=1, model_after="pR2=pL2;pR1=pL2;pL1=pL2", 
#'           progressbar = TRUE))
#' )
#' # When parallel computing is done, no progress bar can be shown
#' system.time(
#' print(Tagloss_L(individuals=data_f_LR, par=par, days.maximum=Tagloss_daymax(data_f_LR), 
#'           fixed.parameters=pfixed, model_after="pR2=pL2;pR1=pL2;pL1=pL2"))
#' )
#' # The NLR_00 are very long to calculate
#' data_f_LR_fast <- subset(data_f_LR, subset=(is.na(data_f_LR$NLR_00)))
#' system.time(
#' print(Tagloss_L(individuals=data_f_LR_fast, par=par, days.maximum=Tagloss_daymax(data_f_LR_fast), 
#'           fixed.parameters=pfixed, model_after="pR2=pL2;pR1=pL2;pL1=pL2"))
#' )
#' o <- Tagloss_fit(data=data_f_LR_fast, 
#'                  fitted.parameters=par, fixed.parameters=pfixed, 
#'                   model_after="pR2=pL2;pR1=pL2;pL1=pL2")
#' 
#' par(mar=c(4, 4, 1, 1))
#' plot(o, t=1:3000, model="2", scale=1000, ylim=c(0, 3), 
#'             col="red")
#' }
#' @export

Tagloss_L <- function(individuals, par, days.maximum=NULL, fixed.parameters=NULL, 
                      model_before=NULL, model_after=NULL, 
                      names.par=NULL, 
                      groups=NULL, 
                      cores=detectCores(all.tests = FALSE, logical = TRUE), 
                      progressbar=FALSE) {
  
  # individuals <- data_f_LR
  # par <- structure(c(72.0399239978454, 58.1034231071992, 645.068735669251, 
  #                    5.10791337470247, 3538.47220045768, 7.83358940767931), 
  #                 .Names = c("D1_L2", "D2D1_L2", "D3D2_L2", "A_L2", "B_L2", "C_L2"))
  # fixed.parameters <- NULL
  # days.maximum=Tagloss_daymax(data_f_LR)
  # model_before=NULL
  # model_after="pR2=pL2;pR1=pL2;pL1=pL2"
  # groups=NULL
  # cores=1
  
  # days.maximum=NULL; fixed.parameters=NULL; model_before=NULL; model_after=NULL; names.par=NULL; groups=NULL; cores=4; progressbar=FALSE
  # par <- structure(c(45.1329605991567, 450.853116508527, -0.00575380257262254, 4.50622019412649, 26.1727077014068, 7.98271155711244), .Names = c("D1_2", "D2D1_2", "D3D2_2", "A_2", "B_2", "C_2"))
  # pfixed <- NULL
  # load(file=file.path("/Users/marcgirondot/Dropbox/Stephanie Kalberer", "dataOut", "m1_1.Rdata"))
  # data_f_21_fast <- m1_1$data
  # par <- m1_1$par
  # individuals <- data_f_21_fast
  
  class(individuals) <- "data.frame"
  
  if (!is.null(names.par)) names(par) <- names.par
  
  if (cores > 1) progressbar=FALSE
  
  Tagloss_Lind <- getFromNamespace(".Tagloss_Lind", ns="phenology")
  
  p1 <- p2 <- pR1 <- pR2 <- pL1 <- pL2 <- Q1 <- Q2 <- Q3 <- Q4 <- Q5 <- Q6 <- Q7 <- Q8 <- LC_Q1 <- NA

  if ((class(individuals) != "matrix") & (class(individuals) != "data.frame")) individuals <- matrix(data = individuals, nrow=1, dimnames = list(NULL, names(individuals)))

  par <- c(par, fixed.parameters)
  
  # L <- NULL
  
  cn <- colnames(individuals)
  totnames <- c("NLR_LR", "NLR_L0", "NLR_0R", "NL0_L0", "N0R_0R", "NL0_00", "N0R_00", "NLR_00", "N22", "N21", "N11", "N10", "N20")
  newnames <- totnames[!(totnames %in% cn)]
  
  if (length(newnames) != 0) {
    individuals <- cbind(individuals, matrix(data=rep(NA, nrow(individuals)*length(newnames)), nrow = nrow(individuals), dimnames = list(NULL, newnames)))
  }
  
  if (is.null(days.maximum)) {
    days.maximum <- Tagloss_daymax(individuals)
  }
  t <- 1:days.maximum
  
  if (!is.null(model_before)) eval(parse(text=model_before), envir= environment())
  
  if (any(grepl("_R2", names(par)))) pR2 <- Tagloss_model(1:days.maximum, par, model_before = model_before, model_after = model_after, model="R2") else pR2 <- NA
  if (any(grepl("_R1", names(par)))) pR1 <- Tagloss_model(1:days.maximum, par, model_before = model_before, model_after = model_after, model="R1") else pR1 <- NA
  if (any(grepl("_L2", names(par)))) pL2 <- Tagloss_model(1:days.maximum, par, model_before = model_before, model_after = model_after, model="L2") else pL2 <- NA
  if (any(grepl("_L1", names(par)))) pL1 <- Tagloss_model(1:days.maximum, par, model_before = model_before, model_after = model_after, model="L1") else pL1 <- NA
  
  if (any(grepl("_2", names(par)))) p2 <- Tagloss_model(1:days.maximum, par, model_before = model_before, model_after = model_after, model="2") else p2 <- NA
  if (any(grepl("_1", names(par)))) p1 <- Tagloss_model(1:days.maximum, par, model_before = model_before, model_after = model_after, model="1") else p1 <- NA
  
  if (!any(grepl("_", names(par)))) p1 <- p2 <- Tagloss_model(1:days.maximum, par, model_before = model_before, model_after = model_after, model=NULL)
  
  if (!is.null(model_after)) eval(parse(text=model_after), envir= environment())
  
  
  if (is.na(p1[1]) & !is.na(p2[1])) p1 <- p2
  if (is.na(p2[1]) & !is.na(p1[1])) p2 <- p1
  
  if (is.na(pR1[1]) & !is.na(pR2[1])) pR1 <- pR2
  if (is.na(pR2[1]) & !is.na(pR1[1])) pR2 <- pR1
  
  if (is.na(pL1[1]) & !is.na(pL2[1])) pL1 <- pL2
  if (is.na(pL2[1]) & !is.na(pL1[1])) pL2 <- pL1
  
  if (!is.na(p1[1]) & !is.na(p2[1])) {
    # Là j'ai le modèle avec N22, N21, N20, N10
    Q1 <- (1 - p2)^2 # probabilité de garder les deux/2
    Q3 <- 1 - Q1 # probabilité de perdre 1/2
    Q2 <- Q3 * p1 # probabilité de perdre 2/2
    Q4 <- 1 - p1 # probabilité de garder 1/1
    Q5 <- p1 # probabilité de perdre 1/1
    if (any(Q1<=0) | any(Q1>=1)) {
      # print(dput(par))
      return(+Inf) 
    }
    LC_Q1 <- cumsum(log(Q1))
  }
  
  if (!is.na(pR2[1]) & !is.na(pR1[1]) & !is.na(pL2[1]) & !is.na(pL1[1])) {
    Q1 <- (1-pL2)*(1-pR2)
    Q2 <- pR2*(1-pL1)
    Q3 <- pL2*(1-pR1)
    Q4 <- pR2 * pL1 + pL2 * pR1
    Q5 <- 1 - pL1
    Q6 <- pL1
    Q7 <- 1 - pR1
    Q8 <- pR1
  }
  
  dfq <- data.frame(p1=p1, p2=p2, pR1=pR1, pR2=pR2, pL1=pL1, pL2=pL2, Q1=Q1, Q2=Q2, Q3=Q3, Q4=Q4, 
                    Q5=Q5, Q6=Q6, Q7=Q7, Q8=Q8, LC_Q1=LC_Q1)
  
  
  if (is.null(groups)) groups <- cores # nrow(individuals)
  
  if ((cores > 1) & (groups != 1)) {
  
  # library("parallel")
    if (groups > nrow(individuals)) groups <- nrow(individuals)
    
    # # print(individu)
    # Je sépare les individus en cores groupes
    nrind <- nrow(individuals)
    
    if (nrind <= groups) {
      group <- as.list(1:nrind)
    } else {
    # Ajouté le 24-09-2017. 6.0 pour les cas où les N20 sont tous groupés sur un coeur
    # Maintenant je les randomise avant de les envoyer
      nindsample <- sample(1:nrind, nrind)
      group <- list()
      nbpargroup <- floor(nrind/groups)
    for (core in 0:(groups-2)) {
      group <- c(group, 
                 list(nindsample[(core*nbpargroup+1):((core+1)*nbpargroup)]))
    }
    # Corrigé le 29-09-2017. 6.0.1
    group <- c(group, 
               list(nindsample[((groups-1)*nbpargroup+1):nrind]))
    }
    
    if (.Platform$OS.type == "unix") {
      L <- mclapply(1:groups,
                    FUN=function(individu) Tagloss_Lind(individuals[group[[individu]], , drop=FALSE],
                                                       dfq, progressbar=FALSE),
                    mc.cores = cores)
      
    } else {
      cl <- makeCluster(cores )
      # If you must use other package in the parallel function; use
      # invisible(clusterEvalQ(cl = cl , library(xxxxxx)))
      L <- parLapply(cl=cl, X=1:groups,
                     fun=function(individu) Tagloss_Lind(individuals[group[[individu]], , drop=FALSE], 
                                                        dfq, progressbar=FALSE)
      )
      stopCluster(cl)
    }
    
  } else {
    L <- Tagloss_Lind(individuals,
                     dfq, 
                     progressbar)
  }


    
    L <- (- sum(unlist(L)))
 
  
  return(L)
  # return(L)
}
