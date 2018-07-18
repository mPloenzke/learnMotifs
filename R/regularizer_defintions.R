#' Sparse group lasso penalty applied to convolutional filters.
#'
#' L1 regularization penalty to constrain individual filter weights plus the L1/L2 regularization to collectively regularize the entire convolutional filter as a group.
#'
#' @param weight_matrix Filter weights.
#'
#' @return Sparse group lasso penalty.
#'
#' @author Matthew Ploenzke, \email{ploenzke@@g.harvard.edu}
#' @seealso \code{\link{position_regularizer}} \code{\link{total_regularizer}} \code{\link{Z_regularizer}}
#' @keywords filter regularizer sparse group
#'
#' @importFrom keras k_sum k_l2_normalize k_sum
#' @export
filter_regularizer <- function(weight_matrix) {
  filter_penalty <- lambda_filter*k_sum(k_l2_normalize(k_sum(weight_matrix,axis=0L),axis=1L))
  weight_penalty <- lambda_l1*k_sum(k_abs(weight_matrix))
  return(filter_penalty+weight_penalty)
}

#' Position-wise lasso penalty applied to convolutional filters weights.
#'
#' L1 regularization penalty to constrain individual filter weights with decaying strength nearing the center of tte filter.
#'
#' @param weight_matrix Filter weights.
#'
#' @return Lasso penalty.
#'
#' @author Matthew Ploenzke, \email{ploenzke@@g.harvard.edu}
#' @seealso \code{\link{filter_regularizer}} \code{\link{total_regularizer}} \code{\link{Z_regularizer}}
#' @keywords regularizer sparse position
#'
#' @importFrom keras k_cast k_concatenate k_arange k_sum
#' @export
position_regularizer <- function(weight_matrix) {
  location_lambda <- k_cast(k_concatenate(list(k_arange(opt$filter_len/2,stop=0,step=-1),
                                               k_arange(start=1,stop=(opt$filter_len/2+1))),
                                          axis=1),'float32')*(lambda_pos/(filter_len/2))
  location_penalty <- k_sum(location_lambda*k_sum(k_abs(weight_matrix),axis=c(1L,3L,4L)))
  return(location_penalty)
}

#' Total regularization function wrapper.
#'
#' Sum of the filter regularization and the position regularization.
#'
#' @param weight_matrix Filter weights.
#'
#' @return Total regularization penalty.
#'
#' @author Matthew Ploenzke, \email{ploenzke@@g.harvard.edu}
#' @seealso \code{\link{filter_regularizer}} \code{\link{position_regularizer}} \code{\link{Z_regularizer}}
#' @keywords regularizer sparse total position filter
#'
#' @export
total_regularizer <- function(weight_matrix) {
  return(filter_regularizer(weight_matrix) +
           position_regularizer(weight_matrix)
  )
}

#' Regularization of the activation values around 0.5.
#'
#' Regularizes the max-pooled activation value (post-sigmoid) around 0.5.
#'
#' @param weight_matrix Filter weights.
#'
#' @return Regularization penalty.
#'
#' @author Matthew Ploenzke, \email{ploenzke@@g.harvard.edu}
#' @seealso \code{\link{filter_regularizer}} \code{\link{position_regularizer}} \code{\link{total_regularizer}}
#' @keywords regularizer sparse activation
#'
#' @importFrom keras k_sum k_max k_abs
#' @export
Z_regularizer <- function(weight_matrix) {
  return(lambda_Z_offset*k_sum(k_abs(k_max(weight_matrix,axis=c(3L))-.5),axis=c(1L,2L,3L)))
}

