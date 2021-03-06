#' Custom 2-D convolution wrapper with to learn IGMs or PWMs.
#'
#' Custom 2-D convolution wrapper with regularization and constraints to learn de novo motifs in the form of IGMs or PWMs.
#'
#' @param object Keras model object.
#' @param filters Number of convolutional filters.
#' @param filter_len Length of the motifs.
#' @param lambda_pos Position-wise L1 sparsity penalty.
#' @param lambda_filter Sparse group lasso penalty.
#' @param lambda_l1 Sparse L1 penalty applied to filter weights.
#' @param lambda_offset L1 penalty applied the offsets (bias) terms.
#' @param strides An integer or list of 2 integers, specifying the strides of the convolution along the width and height. Can be a single integer to specify the same value for all spatial dimensions. Specifying any stride value != 1 is incompatible with specifying any \code{dilation_rate} value != 1.
#' @param padding one of "valid" (default) or "same" (case-insensitive).
#' @param activation Activation function to use. Default is the sigmoid activation.
#' @param use_bias Boolean, whether the layer uses a bias vector. Default is to include a bias.
#' @param kernel_initializer Initializer for the \code{kernel} weights matrix. Default is random uniform on \code{(0,0.5)} (i.e. \code{initializer_random_uniform(0,.5)}).
#' @param bias_initializer Initializer for the bias vector. Default is \code{'zeros'}.
#' @param kernel_regularizer Regularizer function applied to the \code{kernel} weights matrix. Options include \code{NULL}, \code{total_regularizer} (default), \code{position_regularizer}, or \code{filter_regularizer}.
#' @param kernel_constraint Constraint function applied to the kernel matrix. Options include \code{NULL}, \code{info_constraint} (default), \code{pwm_constraint}, \code{negative_constraint}, or \code{nonnegative_constraint}.
#' @param bias_constraint Constraint function applied to the bias vector.Options include \code{NULL}, \code{negative_constraint} (default), or \code{nonnegative_constraint}.
#' @param input_shape Dimensionality of the input (integer) not including the samples axis. This argument is required when using this layer as the first layer in a model.
#' @param name An optional name string for the layer. Should be unique in a model (do not reuse the same name twice). It will be autogenerated if it isn't provided. Default is \code{'deNovo_conv'}.
#' @param trainable Whether the layer weights will be updated during training. Must be specified.
#'
#' @return A Keras model with the added layer.
#'
#' @author Matthew Ploenzke, \email{ploenzke@@g.harvard.edu}
#' @keywords conv convolution layer denovo
#'
#' @export
layer_deNovo <- function (object,
                          filters,
                          filter_len,
                          lambda_pos=0,
                          lambda_filter=0,
                          lambda_l1=0,
                          lambda_offset=0,
                          strides = c(1L, 1L),
                          padding = "valid",
                          activation = 'sigmoid',
                          use_bias = TRUE,
                          kernel_initializer = initializer_random_uniform(0,.5),
                          bias_initializer = "zeros",
                          kernel_regularizer = total_regularizer,
                          kernel_constraint = info_constraint,
                          bias_constraint = negative_constraint,
                          input_shape = NULL,
                          name = 'deNovo_conv', trainable = NULL) {
  lambda_filter <<- lambda_filter
  lambda_l1 <<- lambda_l1
  lambda_pos <<- lambda_pos
  filter_len <<- filter_len
  keras::create_layer(keras:::keras$layers$Conv2D, object, list(filters = as.integer(filters),
                                                 kernel_size = keras:::as_integer_tuple(c(4,filter_len)), strides = keras:::as_integer_tuple(strides),
                                                 padding = padding, data_format = NULL, dilation_rate = c(1L, 1L),
                                                 activation = activation, use_bias = use_bias, kernel_initializer = kernel_initializer,
                                                 bias_initializer = bias_initializer, kernel_regularizer = kernel_regularizer,
                                                 bias_regularizer = keras::regularizer_l1(l=lambda_offset), activity_regularizer = NULL,
                                                 kernel_constraint = kernel_constraint, bias_constraint = bias_constraint,
                                                 input_shape = keras:::normalize_shape(input_shape), batch_input_shape = keras:::normalize_shape(NULL),
                                                 batch_size = keras:::as_nullable_integer(NULL), dtype = NULL,
                                                 name = name, trainable = trainable, weights = NULL))
}

