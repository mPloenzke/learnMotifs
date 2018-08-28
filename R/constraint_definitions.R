#' Constrains filter weights within information bounds.
#'
#' Constrains the individual filters weights to sum column-wise less than or equal to 2.
#'
#' @param w Filter weights.
#'
#' @return Constraint.
#'
#' @author Matthew Ploenzke, \email{ploenzke@@g.harvard.edu}
#' @keywords information constraint
#'
#' @export
info_constraint <- function(w) {
  w <- w * k_cast(k_greater_equal(w, 0), k_floatx()) # force nonnegative
  w * (2/k_maximum(k_sum(w,axis=1,keepdims=T),2)) # force sum less-than equal to two (bits)
}
#' Constrains weights to be negative.
#'
#' Constrains weights to be negative.
#'
#' @param w Weights.
#'
#' @return Constraint.
#'
#' @author Matthew Ploenzke, \email{ploenzke@@g.harvard.edu}
#' @keywords negative constraint
#'
#' @export
negative_constraint <- function(w) {
  w <- w * k_cast(k_less_equal(w, 0), k_floatx()) # force negative
}
#' Constrains weights to be non-negative.
#'
#' Constrains weights to be non-negative.
#'
#' @param w Weights.
#'
#' @return Constraint.
#'
#' @author Matthew Ploenzke, \email{ploenzke@@g.harvard.edu}
#' @keywords negative constraint
#'
#' @export
nonnegative_constraint <- function(w) {
  w <- w * k_cast(k_greater_equal(w, 0), k_floatx()) # force nonnegative
}
#' Combination constraint to set non-negative and fix non-de novo weights.
#'
#' Constrains weights to be non-negative and fixes the weights for the fixed motifs.
#'
#' @param w Weights.
#'
#' @return Constraint.
#'
#' @author Matthew Ploenzke, \email{ploenzke@@g.harvard.edu}
#' @keywords combo constraint
#'
#' @importFrom keras k_cast k_greater_equal k_floatx k_expand_dims k_concatenate
#' @export
combo_constraint <- function(w) {
  w <- w * k_cast(k_greater_equal(w, 0), k_floatx()) # force nonnegative
  w <- k_expand_dims(k_concatenate(list(w[1:opt$n_filters,1],beta.initializer[-(1:opt$n_filters)]),axis=1),axis=-1)
}
#' Constrains weights to fixed motifs.
#'
#' Constrains weights to to the fixed motifs. These may be the annotated database motifs or perhaps de novo motifs.
#'
#' @param w Weights.
#'
#' @return Constraint.
#'
#' @author Matthew Ploenzke, \email{ploenzke@@g.harvard.edu}
#' @keywords fixed motif constraint
#'
#' @export
fixedMotif_constraint <- function(w) {
  w <- fixed.motifs
}
#' Constrains weights to be valid information measures.
#'
#' Constrains filter weights to be valid ICMs by rescaling the weights.
#'
#' @param wts Convolutional filter weights.
#' @param offset Filter offset weights.
#' @param bg Vector of background nucleotide probabilities (A, C, G, T)
#'
#' @return Rescaled weights.
#'
#' @author Matthew Ploenzke, \email{ploenzke@@g.harvard.edu}
#' @keywords ICM constraint information
#'
#' @export
restrict_to_ICM <- function(wts,offset,bg=c(.25,.25,.25,.25)) {
  curr.info <- apply(wts,4,sum)
  newwts <- lapply(1:dim(wts)[4], function(j) {
    i <- wts[,,,j,drop=FALSE]
    i[,colSums(i)<1,1,1] <- i[,colSums(i)<1,1,1] + .05
    # convert to ppm
    ppm <- (i)/array(matrix(colSums(i),nrow=4,ncol=ncol(i),byrow=T),dim=dim(i))
    ppm[is.nan(ppm)] <- .25
    # convert back to valid icm
    #R_i = pmax(2 + apply(ppm, 2, function(col) sum(col * log2(col),na.rm=T)),0)
    R_i = pmax(apply(ppm, 2, function(col) sum(col * log2(col/bg),na.rm=T)),0)
    icm <- ppm*array(matrix(R_i,nrow=4,ncol=ncol(ppm),byrow=T),dim=dim(ppm))
  })
  newwts <- do.call(abind,newwts)
  new.info <- apply(newwts,4,sum)
  newoffset <- offset*(new.info/ifelse(curr.info==0,0.1,curr.info))
  return(list(newwts,newoffset))
}
