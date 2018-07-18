#' Convert PFM to ICM.
#'
#' Convert position frequency matrix into an information content matrix to be consistent with the de novo filters.
#'
#' @param motif.raw List of annotated motifs as PFMs.
#' @param motif.pos Vector of indices to extract motifs. If \code{NULL} all motifs will be selected.
#' @param mmax Maximum motif length, found via \code{\link{find_max_length}}.
#'
#' @return 4-D array containing ICMs of the extracted motifs.
#'
#' @author Matthew Ploenzke, \email{ploenzke@@g.harvard.edu}
#' @seealso \code{\link{find_max_length}}
#' @keywords pfm icm information content motif
#'
#' @examples
#'
#'
#' @export
pfm_to_icm <- function(motif.raw, motif.pos, mmax) {
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
    R_i = pmax(2 + apply(ppm, 2, function(col) sum(col * log2(col),na.rm=T)),0)
    annotated.motifs[,,,j] <- ppm*matrix(R_i,nrow=4,ncol=mmax,byrow=T)
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
#' @seealso \code{\link{pfm_to_icm}}
#' @keywords pfm icm information content motif
#'
#' @examples
#'
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
