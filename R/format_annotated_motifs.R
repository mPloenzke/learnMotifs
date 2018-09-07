#' Convert PFM to IGM.
#'
#' Convert position frequency matrix into an information content matrix to be consistent with the de novo filters.
#'
#' @param motif.raw List of annotated motifs as PFMs.
#' @param motif.pos Vector of indices to extract motifs. If \code{NULL} all motifs will be selected.
#' @param mmax Maximum motif length, found via \code{\link{find_max_length}}.
#' @param bg Vector of background nucleotide probabilities (A, C, G, T)
#'
#' @return 4-D array containing IGMs of the extracted motifs.
#'
#' @author Matthew Ploenzke, \email{ploenzke@@g.harvard.edu}
#' @seealso \code{\link{find_max_length}}
#' @keywords pfm IGM information content motif
#'
#' @export
pfm_to_igm <- function(motif.raw, motif.pos=NULL, mmax, bg=c(.25,.25,.25,.25)) {
  if (is.null(motif.pos)) {motif.pos <- 1:length(motif.raw)}
  annotated.motifs <- array(NA,dim=c(4,mmax,1,length(motif.pos)))
  for (j in 1:length(motif.pos)) {
    pfm <- as.matrix(motif.raw[[motif.pos[j]]])
    ppm <- pfm/matrix(colSums(pfm),nrow=4,ncol=ncol(pfm),byrow=T)
    if (ncol(pfm)<mmax) {
      new <- matrix(0,nrow=4,ncol=mmax-ncol(pfm))
      pfm <- cbind(pfm,new)
      ppm <- cbind(ppm,new)
    }
    R_i = pmax(apply(ppm, 2, function(col) sum(col * log2(col/bg),na.rm=T)),0)
    annotated.motifs[,,,j] <- ppm*matrix(R_i,nrow=4,ncol=mmax,byrow=T)
  }
  return(annotated.motifs)
}
#' Convert PFM to PWM.
#'
#' Convert position frequency matrix into a position weight matrix.
#'
#' @param motif.raw List of annotated motifs as PFMs.
#' @param motif.pos Vector of indices to extract motifs. If \code{NULL} all motifs will be selected.
#' @param mmax Maximum motif length, found via \code{\link{find_max_length}}.
#' @param bg Vector of background nucleotide probabilities (A, C, G, T)
#' @param pseudoprob Value to add to zero-valued position probability matrix entries prior to converting to PWM.
#'
#' @return 4-D array containing PWMs of the extracted motifs.
#'
#' @author Matthew Ploenzke, \email{ploenzke@@g.harvard.edu}
#' @seealso \code{\link{find_max_length}}
#' @keywords pfm pwm position weight motif
#'
#' @export
pfm_to_pwm <- function(motif.raw, motif.pos=NULL, mmax, bg=c(.25,.25,.25,.25), pseudoprob=1e-3) {
  if (is.null(motif.pos)) {motif.pos <- 1:length(motif.raw)}
  annotated.motifs <- array(NA,dim=c(4,mmax,1,length(motif.pos)))
  for (j in 1:length(motif.pos)) {
    pfm <- as.matrix(motif.raw[[motif.pos[j]]])
    ppm <- pfm/matrix(colSums(pfm),nrow=4,ncol=ncol(pfm),byrow=T)
    if (ncol(pfm)<mmax) {
      new <- matrix(0,nrow=4,ncol=mmax-ncol(pfm))
      pfm <- cbind(pfm,new)
      ppm <- cbind(ppm,new)
    }
    ppm[ppm==0] <- pseudoprob
    annotated.motifs[,,,j] <- log2(ppm/matrix(bg,nrow=4,ncol=mmax,byrow=F))
  }
  return(annotated.motifs)
}
#' Find maximum motif length.
#'
#' Loop through all motifs and return the length of the longest motif.
#'
#' @param motif.raw List of annotated motifs as PFMs.
#' @param motif.pos Vector of indices to extract motifs. If \code{NULL} all motifs will be selected.
#'
#' @return Length of the longest motif.
#'
#' @author Matthew Ploenzke, \email{ploenzke@@g.harvard.edu}
#' @seealso \code{\link{pfm_to_igm}}
#' @keywords pfm igm information gain motif
#'
#' @export
find_max_length <- function(motif.raw, motif.pos) {
  mmax <- 1
  if (is.null(motif.pos)) {motif.pos <- 1:length(motif.raw)}
  for(j in motif.pos) {
    mmax <- ifelse(ncol(motif.raw[[j]])>mmax,ncol(motif.raw[[j]]),mmax)
  }
  return(mmax)
}
