#' Top max-pooled hits
#'
#' Calculate the top N motif hits max-pooled over the entire sequence. This function uses the \code{keras} \code{k_max} function which internally
#' converts the input to its pre-activation scale. Therefore we recommend using the \code{linear} activation for the layer feeding into
#' this one and then using the \code{layer_activation} function after the operation.
#'
#' @param object Keras model.
#' @param n_max Number of top hits to return.
#' @param name Layer name (optional).
#'
#' @return Top max-pooled hits layer.
#'
#' @author Matthew Ploenzke, \email{ploenzke@@g.harvard.edu}
#' @seealso \code{\link{topN_max_pooling}}
#' @keywords motif hits max pooling
#'
#' @export
layer_topN_max_pooling <- function(object, n_max, name = NULL) {
  keras::create_layer(topN_max_pooling, object, list(
    n_max = as.integer(n_max),
    name = name
  ))
}
#' Top max-pooled hits R6 class wrapper
#'
#' R6 class container for the top max-pooled motif hits
#'
#' @author Matthew Ploenzke, \email{ploenzke@@g.harvard.edu}
#' @keywords callback model
#'
#' @importFrom keras k_max k_cast k_not_equal k_concatenate
#' @export
topN_max_pooling <- R6::R6Class("topN_max_pooling",
                                inherit = KerasLayer,
                                public = list(
                                  n_max = NULL,
                                  initialize = function(n_max) {
                                    self$n_max <- n_max
                                  },
                                  call = function(x, mask = NULL) {
                                    a <- k_max(x, axis=3, keepdims = TRUE)
                                    x <- x*k_cast(k_not_equal(x,a),dtype = 'float')
                                    if (self$n_max > 1) {
                                      for (maxnum in 2:self$n_max) {
                                        b <- k_max(x, axis=3, keepdims = TRUE)
                                        x <- x*k_cast(k_not_equal(x,b),dtype = 'float')
                                        a <- k_concatenate(c(a, b), axis = 3)
                                      }
                                    }
                                    return(a)
                                  },
                                  compute_output_shape = function(input_shape) {
                                    input_shape[[3]] <- self$n_max
                                    input_shape
                                  }
                                )
)
