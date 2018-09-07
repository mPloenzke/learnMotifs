#' Calculate offset
#'
#' Calculate the offset (bias) as the mean activation over a sampling of the data.
#'
#' @param filter.motifs 4-D array of IGMs or PWMs.
#' @param motif.pos Vector of indices to extract motifs. If \code{NULL} all motifs will be selected.
#' @param motif.maxlen Maximum motif length, found via \code{\link{find_max_length}}.
#' @param onehot_array One-hot encoded training data to sample a batch from.
#' @param batchsize Size of batch to sample.
#'
#' @return Vector of offsets of length equal to the fourth dimension of the input array of IGMs/PWMs.
#'
#' @author Matthew Ploenzke, \email{ploenzke@@g.harvard.edu}
#' @seealso \code{\link{find_max_length}}
#' @keywords offset bias
#'
#' @importFrom keras keras_model_sequential layer_conv_2d layer_max_pooling_2d layer_flatten compile
#' @export
calculate_offsets <- function(filter.motifs, motif.pos=NULL, motif.maxlen, onehot_array, batchsize=256) {
  if (is.null(motif.pos)) {motif.pos <- 1:dim(filter.motifs)[4]}
  model_offset <- keras_model_sequential()
  model_offset %>%
    layer_conv_2d(filters = length(motif.pos),
                  kernel_size = c(4,motif.maxlen),
                  kernel_initializer = initializer_constant(filter.motifs),
                  kernel_constraint = fixedMotif_constraint,
                  use_bias = FALSE,
                  input_shape = c(4, opt$max_len+2*motif.maxlen-2, 1),
                  padding = "valid",
                  activation='linear') %>%
    layer_max_pooling_2d(pool_size = c(1,(opt$max_len+motif.maxlen-1))) %>%
    layer_flatten() %>%
    compile(optimizer = 'sgd',loss = 'binary_crossentropy',metrics = 'accuracy')
  offsets <- colMeans(predict(model_offset,array(onehot_array[sample(1:nrow(onehot_array),batchsize),,,],
                                                 dim=c(batchsize,dim(onehot_array)[2:4]))))
  return(offsets)
}
