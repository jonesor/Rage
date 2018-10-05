#' Derive mean vital rates from a matrix population model
#' 
#' Derive mean vital rates from a matrix population model corresponding to
#' separate demographic processes. Specifically, this function decomposes vital
#' rates of survival, progression, retrogression, sexual reproduction and clonal
#' reproduction according to various ways of weighted means and organization of
#' stages along the life cycle represented in the matrix population model.
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
#' @param splitStages Splits vital rates according to some pre-determined
#'   criteria (below).
#' @param weighted Allows to weight mean vital rates according to various
#'   criteria (below).
#' @return - 'Weighted': This argument allows to weight mean values of vital
#' rates (survival 'surv', progression 'prog', retrogression 'retr', sexual
#' reproduction 'fec' and clonal reproduction 'clo') with an equal contribution
#' for all stages (default), by the stable st/age distribution ('SSD'), or by a
#' given population vector chosen by the user, so long as it is congruent with
#' the dimensions of the chosen 'matU', 'matF', and 'matC'.
#' 
#' - 'splitStages': This argument allows to split the values of vital rates
#' according to recognizable stages in the matrix. When 'all', all vital rates
#' are averaged all existing stages, if 'ontogeny', they are averaged as
#' juveniles ('Juv') and adults ('Adu'), and if by 'MatrixClassOrganized', it
#' takes a vector with the pre-established stages of
#' 'compadre$matrixClass[[i]]$MatrixClassOrganized' or
#' 'compadre$matrixClass[[i]]$MatrixClassOrganized', where 'i' is the index of
#' the chosen study in 'COMPADRE' or 'COMADRE'.
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
#' #Vital rate outputs without weights:
#' vitalRates(matU, matF, matC, splitStages = 'all', weighted = FALSE)
#' vitalRates(matU, matF, matC, splitStages = 'ontogeny', weighted = FALSE)
#' vitalRates(matU, matF, matC,
#'            splitStages = c('prop', 'active', 'active', 'active'),
#'            weighted = FALSE)
#' 
#' 
#' #Vital rate outputs weighted by the stable stage distribution of 'matA':
#' vitalRates(matU, matF, matC, splitStages = 'all', weighted = 'SSD')
#' vitalRates(matU, matF, matC, splitStages = 'ontogeny', weighted = 'SSD')
#' vitalRates(matU, matF, matC,
#'            splitStages = c('prop', 'active', 'active', 'active'),
#'            weighted = 'SSD')
#' 
#' @importFrom popbio stable.stage
#' @export vitalRates
vitalRates <- function(matU, matF, matC = NULL, splitStages = FALSE,
                       weighted = FALSE) {
  #Function to quantify vital rates values
  
  if (missing(matU)) {
    stop('matU missing')
  }
  if (missing(matF) & missing(matC)) {
    warning('matF or matC missing. These have been coerced to matrices of zero')
  }
  # if (sum(weighted)>0 & length(weighted) != dim(matU)[i]) {
  #   stop('Population vector does not agree with matrix dimension')
  # }
  if (is.null(matC)) {
    matC <- matrix(0, nrow = dim(matU)[1], ncol = dim(matU)[1])
  }
  
  matDim <- dim(matU)[1]
  matA <- matU + matF + matC
  
  surv <- colSums(matU)
  fec <- colSums(matF)
  clo <- colSums(matC)
  
  matUIndep <- matrix(NA, matDim, matDim)
  for (i in 1:matDim) {matUIndep[,i] <- matU[,i] / surv[i]}
  prog <- retr <- matUIndep
  prog[is.nan(prog)] <- 0
  retr[is.nan(retr)] <- 0
  
  prog[which(upper.tri(matUIndep, diag = TRUE))] <- 0
  retr[which(lower.tri(matUIndep, diag = TRUE))] <- 0
  
  prog <- colSums(prog)
  retr <- colSums(retr)
  
  if (weighted[1] == FALSE) {
    weight <- rep(1,matDim)
  }
  
  if (weighted[1] == 'SSD') {
    weight <- popbio::stable.stage(matA)
  }
  
  weight <- weight / sum(weight)
  
  surv1 <- surv * weight
  fec1  <- fec  * weight
  clo1  <- clo  * weight
  prog1 <- prog * weight
  retr1 <- retr * weight
  
  out <- NULL
  
  if (splitStages[1] == 'all') {
    out$surv <- sum(surv1)
    out$retr <- sum(retr1)
    out$prog <- sum(prog1)
    out$fec  <- sum(fec1)
    out$clo  <- sum(clo1)
  }
  
  if (splitStages[1] == 'ontogeny') {
    #This adu classification does not account for non- and post-reproductive
    adu <- which(colSums(matF) > 0)
    juv <- which(colSums(matF) == 0)
    
    out$survJuv <- mean(surv1[juv], na.rm=TRUE)
    out$retrJuv <- mean(retr1[juv], na.rm=TRUE)
    out$progJuv <- mean(prog1[juv], na.rm=TRUE)
    out$cloJuv  <- mean(clo1[juv], na.rm=TRUE)
    
    out$survAdu <- mean(surv1[adu], na.rm=TRUE)
    out$retrAdu <- mean(retr1[adu], na.rm=TRUE)
    out$progAdu <- mean(prog1[adu], na.rm=TRUE)
    out$fecAdu  <- mean(fec1[adu], na.rm=TRUE)
    out$cloAdu  <- mean(clo1[adu], na.rm=TRUE)
    }
  
  if (splitStages[1] %in% c('prop','active','dorm')) {
    prop <- which(splitStages == "prop")
    active <- which(splitStages == "active")
    dorm <- which(splitStages == "dorm")
    
    out$survProp <- mean(surv1[prop], na.rm=TRUE)
    out$progProp <- mean(prog1[prop], na.rm=TRUE)
    
    out$survActive <- mean(surv1[active], na.rm=TRUE)
    out$retrActive <- mean(retr1[active], na.rm=TRUE)
    out$progActive <- mean(prog1[active], na.rm=TRUE)
    out$fecActive  <- mean(fec1[active], na.rm=TRUE)
    out$cloActive  <- mean(clo1[active], na.rm=TRUE)
    
    out$survDorm <- mean(surv1[dorm], na.rm=TRUE)
    out$retrDorm <- mean(retr1[dorm], na.rm=TRUE)
    out$progDorm <- mean(prog1[dorm], na.rm=TRUE)
  }
  
	return(out)
 }
