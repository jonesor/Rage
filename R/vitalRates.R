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
#'   progression, stasis, and retrogression).
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
#' 
#' @references Caswell, H. 2001. Matrix Population Models: Construction,
#'   Analysis, and Interpretation. Sinauer Associates; 2nd edition. ISBN:
#'   978-0878930968
#'   
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
  
  # possible transitions
  i_surv <- surv > 0
  i_prog <- prog > 0
  i_retr <- retr > 0
  i_fec <- fec > 0
  i_clo <- clo > 0
  
  if (is.null(weights)) {
    weights <- rep(1.0, matDim)
  } else if (weights[1] == "SSD") {
    weights <- popbio::stable.stage(matA)
  }
  
  weights <- weights / sum(weights)
  
  surv1 <- surv * weights
  prog1 <- prog * weights
  retr1 <- retr * weights
  fec1  <- fec  * weights
  clo1  <- clo  * weights
  
  out <- NULL
  
  if (splitStages == "all") {
    out$surv <- sum(surv1[i_surv]) / sum(weights[i_surv])
    out$retr <- sum(retr1[i_retr]) / sum(weights[i_retr])
    out$prog <- sum(prog1[i_prog]) / sum(weights[i_prog])
    out$fec  <- sum(fec1[i_fec])   / sum(weights[i_fec])
    out$clo  <- sum(clo1[i_clo])   / sum(weights[i_clo])
  }
  
  if (splitStages == "ontogeny") {
    #This adu classification does not account for non- and post-reproductive
    adu <- colSums(matF) > 0
    juv <- colSums(matF) == 0
    
    out$survJuv <- sum(surv1[juv & i_surv]) / sum(weights[juv & i_surv])
    out$retrJuv <- sum(retr1[juv & i_retr]) / sum(weights[juv & i_retr])
    out$progJuv <- sum(prog1[juv & i_prog]) / sum(weights[juv & i_prog])
    out$cloJuv  <- sum(clo1[juv & i_clo])   / sum(weights[juv & i_clo])
    
    out$survAdu <- sum(surv1[adu & i_surv]) / sum(weights[adu & i_surv])
    out$retrAdu <- sum(retr1[adu & i_retr]) / sum(weights[adu & i_retr])
    out$progAdu <- sum(prog1[adu & i_prog]) / sum(weights[adu & i_prog])
    out$fecAdu  <- sum(fec1[adu & i_fec])   / sum(weights[adu & i_fec])
    out$cloAdu  <- sum(clo1[adu & i_clo])   / sum(weights[adu & i_clo])
  }
  
  if (splitStages == "matrixStages") {
    prop <- matrixStages == "prop"
    acti <- matrixStages == "active"
    dorm <- matrixStages == "dorm"
    
    out$survProp <- sum(surv1[prop & i_clo]) / sum(weights[prop & i_clo])
    out$progProp <- sum(prog1[prop & i_clo]) / sum(weights[prop & i_clo])
    
    out$survActive <- sum(surv1[acti & i_surv]) / sum(weights[acti & i_surv])
    out$retrActive <- sum(retr1[acti & i_retr]) / sum(weights[acti & i_retr])
    out$progActive <- sum(prog1[acti & i_prog]) / sum(weights[acti & i_prog])
    out$fecActive  <- sum(fec1[acti & i_fec]) / sum(weights[acti & i_fec])
    out$cloActive  <- sum(clo1[acti & i_clo]) / sum(weights[acti & i_clo])
    
    out$survDorm <- sum(surv1[dorm & i_surv]) / sum(weights[dorm & i_surv])
    out$retrDorm <- sum(retr1[dorm & i_retr]) / sum(weights[dorm & i_retr])
    out$progDorm <- sum(prog1[dorm & i_prog]) / sum(weights[dorm & i_prog])
  }
  
  # convert NaN to 0 (non-existent vital rate within group will yield NaN)
  out <- lapply(out, function(x) ifelse(is.nan(x), 0.0, x))
  
  return(out)
}
