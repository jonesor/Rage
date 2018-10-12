#' Derive mean vital rates from a matrix population model
#' 
#' Derive mean vital rates corresponding to separate demographic processes from
#' a matrix population model. Specifically, this function decomposes vital rates
#' of survival, progression, retrogression, sexual reproduction and clonal
#' reproduction, with various options for weighting and grouping stages of the
#' life cycle.
#' 
#' @param matU The survival component of a matrix population model (i.e. a
#'   square projection matrix reflecting survival-related transitions; e.g.
#'   progression, stasis, and retrogression)
#' @param matF The sexual component of a matrix population model (i.e. a square
#'   projection matrix reflecting transitions due to sexual reproduction)
#' @param matC The clonal component of a matrix population model (i.e. a square
#'   projection matrix reflecting transitions due to clonal reproduction).
#'   Defaults to \code{NULL}, indicating no clonal reproduction (i.e.
#'   \code{matC} is a matrix of zeros).
#' @param weights Vector of stage-specific weights to apply while averaging
#'   vital rates. Default is \code{NULL} reflecting equal weighting for all
#'   stages. May also be \code{"SSD"} to weight vital rates by the stable
#'   distribution of \code{matA}.
#' @param splitStages What groups should vital rates be averaged over. Either:
#' 
#' \code{"all"}: all stages grouped
#' 
#' \code{"ontogeny"}: group juvenile stages (all stages prior to the first stage
#' with sexual reproduction) and adult stages
#' 
#' \code{"matrixStages"}: group according to a standardized set of stage classes
#' (propagule, active, and dormant). If \code{splitStages = "matrixStages"},
#' must also specify separate argument \code{matrixStages}.
#' 
#' @param matrixStages Vector of stage-specific standardized matrix classes
#'   ("prop" for propagule, "active", and/or "dorm" for dormant). Only used if
#'   \code{splitStages = "matrixClass"}.
#' @return A list of averaged vital rates.
#' 
#' @author Roberto Salguero-Gomez <rob.salguero@@zoo.ox.ac.uk>
#' @references Caswell, H. (2001) Matrix Population Models: Construction,
#'   Analysis, and Interpretation. Sinauer Associates; 2nd edition. ISBN:
#'   978-0878930968
#' @examples
#' matU <- rbind(c(0.1,   0,   0,   0),
#'               c(0.5, 0.2, 0.1,   0),
#'               c(  0, 0.3, 0.3, 0.1),
#'               c(  0,   0, 0.5, 0.6))
#' 
#' matF <- rbind(c(  0,   0, 1.1, 1.6),
#'               c(  0,   0, 0.8, 0.4),
#'               c(  0,   0,   0,   0),
#'               c(  0,   0,   0,   0))
#'               
#' matC <- rbind(c(  0,   0, 0.4, 0.5),
#'               c(  0,   0, 0.3, 0.1),
#'               c(  0,   0,   0,   0),
#'               c(  0,   0,   0,   0))
#' 
#' # Vital rate outputs without weights
#' vitalRates(matU, matF, matC, splitStages = 'all')
#' vitalRates(matU, matF, matC, splitStages = 'ontogeny')
#' 
#' # Group vital rates according to specified matrixStages
#' ms <- c('prop', 'active', 'active', 'active')
#' vitalRates(matU, matF, matC, splitStages = 'matrixStages', matrixStages = ms)
#' 
#' # Vital rate outputs weighted by the stable stage distribution of 'matA'
#' vitalRates(matU, matF, matC, splitStages = 'all', weights = 'SSD')
#' 
#' @importFrom popbio stable.stage
#' @export vitalRates
vitalRates <- function(matU, matF, matC = NULL, weights = NULL,
                       splitStages = "all", matrixStages = NULL) {
  
  # validate arguments
  checkValidMat(matU)
  checkValidMat(matF)
  if (!is.null(matC)) checkValidMat(matC, warn_all_zero = FALSE)
  if (!is.null(weights) && weights != "SSD" &&
        length(weights) != nrow(matU)) {
    stop("If weights are provided, length(weights) should be of the same ",
         "dimension as matU", call. = FALSE)
  }
  if (!splitStages %in% c("all", "ontogeny", "matrixStages")) {
    stop("Argument splitStages must be one of 'all', 'ontogeny', or ",
         "'matrixStages'", call. = FALSE)
  }
  if (splitStages == "matrixStages") {
    if (is.null(matrixStages)) {
      stop("If splitStages = 'matrixStages', argument matrixStages must be ",
           "provided", call. = FALSE)
    }
    if (length(matrixStages) != nrow(matU)) {
      stop("length(matrixStages) should be of the same dimension as matU",
           call. = FALSE)
    }
  }
  
  if (is.null(matC)) {
    matC <- matrix(0, nrow = nrow(matU), ncol = ncol(matU))
  }
  
  matDim <- dim(matU)[1]
  matA <- matU + matF + matC
  
  surv <- colSums(matU)
  fec <- colSums(matF)
  clo <- colSums(matC)
  
  matUIndep <- matrix(NA_real_, matDim, matDim)
  for (i in 1:matDim) {matUIndep[,i] <- matU[,i] / surv[i]}
  prog <- retr <- matUIndep
  prog[is.nan(prog)] <- 0
  retr[is.nan(retr)] <- 0
  
  prog[which(upper.tri(matUIndep, diag = TRUE))] <- 0
  retr[which(lower.tri(matUIndep, diag = TRUE))] <- 0
  
  prog <- colSums(prog)
  retr <- colSums(retr)
  
  if (is.null(weights)) {
    weights <- rep(1.0, matDim)
  } else if (weights[1] == "SSD") {
    weights <- popbio::stable.stage(matA)
  }
  
  weight <- weights / sum(weights)
  
  surv1 <- surv * weight
  fec1  <- fec  * weight
  clo1  <- clo  * weight
  prog1 <- prog * weight
  retr1 <- retr * weight
  
  out <- NULL
  
  if (splitStages == "all") {
    out$surv <- sum(surv1)
    out$retr <- sum(retr1)
    out$prog <- sum(prog1)
    out$fec  <- sum(fec1)
    out$clo  <- sum(clo1)
  }
  
  if (splitStages == "ontogeny") {
    #This adu classification does not account for non- and post-reproductive
    adu <- which(colSums(matF) > 0)
    juv <- which(colSums(matF) == 0)
    
    out$survJuv <- mean(surv1[juv], na.rm = TRUE)
    out$retrJuv <- mean(retr1[juv], na.rm = TRUE)
    out$progJuv <- mean(prog1[juv], na.rm = TRUE)
    out$cloJuv  <- mean(clo1[juv], na.rm = TRUE)
    
    out$survAdu <- mean(surv1[adu], na.rm = TRUE)
    out$retrAdu <- mean(retr1[adu], na.rm = TRUE)
    out$progAdu <- mean(prog1[adu], na.rm = TRUE)
    out$fecAdu  <- mean(fec1[adu], na.rm = TRUE)
    out$cloAdu  <- mean(clo1[adu], na.rm = TRUE)
  }
  
  if (splitStages == "matrixStages") {
    prop <- which(matrixStages == "prop")
    active <- which(matrixStages == "active")
    dorm <- which(matrixStages == "dorm")
    
    out$survProp <- mean(surv1[prop], na.rm = TRUE)
    out$progProp <- mean(prog1[prop], na.rm = TRUE)
    
    out$survActive <- mean(surv1[active], na.rm = TRUE)
    out$retrActive <- mean(retr1[active], na.rm = TRUE)
    out$progActive <- mean(prog1[active], na.rm = TRUE)
    out$fecActive  <- mean(fec1[active], na.rm = TRUE)
    out$cloActive  <- mean(clo1[active], na.rm = TRUE)
    
    out$survDorm <- mean(surv1[dorm], na.rm = TRUE)
    out$retrDorm <- mean(retr1[dorm], na.rm = TRUE)
    out$progDorm <- mean(prog1[dorm], na.rm = TRUE)
  }
  
  return(out)
}
