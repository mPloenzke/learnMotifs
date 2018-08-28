#' One-hot encoding fit generator
#'
#' Streaming generator to one-hot encode samples during training or predicting.
#'
#' @param X_data Character strings composed of nucleotides A, C, G, T.
#' @param Y_data Vector of labels.
#' @param batch_size Integer value denoting batch size.
#' @param zeros_len Integer value denoting number of zeros to pad vector with. Commonly the convolutional filter length - 1.
#' @param max_len Maximum sequence length. Calculated per batch if NULL.
#'
#' @return Zero-padded array of one-hot encodings.
#'
#' @author Matthew Ploenzke, \email{ploenzke@@g.harvard.edu}
#' @seealso \code{\link{one_hot_vector}}
#' @keywords one-hot generator
#'
#' @export
one_hot_generator <- function(X_data, Y_data, batch_size, zeros_len, max_len) {
  function() {
    rows <- sample(1:length(X_data), batch_size, replace = TRUE)
    X_array <- one_hot(X_data[rows],zeros_len,max_len)
    list(X_array, Y_data[rows])
  }
}

#' One-hot encoding of vector of DNA nucleotide strings
#'
#' Inputs vector of strings composed of nucleotides and outputs a zero-padded array of encodings.
#'
#' @param data Character strings composed of nucleotides A, C, G, T.
#' @param zeros_len Integer value denoting number of zeros to pad vector with. Commonly the convolutional filter length - 1.
#' @param max_len Maximum sequence length. If unspecified defaults to the maximum number of characters present in the data character strings.
#'
#' @return Zero-padded array of one-hot encodings.
#'
#' @author Matthew Ploenzke, \email{ploenzke@@g.harvard.edu}
#' @seealso \code{\link{one_hot_vector}}
#' @keywords one-hot
#'
#' @examples
#' one_hot(rep(paste(sample(c('A','C','G','T'),100,replace=TRUE),collapse=''),50),10)
#'
#' @export
one_hot <- function(data,zeros_len,max_len=NULL) {
  if (is.null(max_len)) {max_len <- max(nchar(data))}
  tensor <- aperm(array(c(sapply(data,one_hot_vector,zeros_len=zeros_len,max_len=max_len)),dim=c((max_len+2*zeros_len-2),4,length(data),1)),perm=c(3,2,1,4))
  return(tensor)
}

#' One-hot encoding of DNA nucleotide string
#'
#' Inputs a string composed of nucleotides and outputs a zero-padded vector of encodings.
#'
#' @param dnastr Character string composed of nucleotides A, C, G, T.
#' @param zeros_len Integer value denoting number of zeros to pad vector with. Commonly the convolutional filter length - 1.
#' @param max_len Integer value denoting the maximum sequence length. Used in the padding calculation to account for shorter sequences.
#'
#' @return Zero-padded vector of one-hot encodings. If the input is of length \eqn{L} and with padding of length Z the output is of length \eqn{4(L+2Z)}.
#'
#' @author Matthew Ploenzke, \email{ploenzke@@g.harvard.edu}
#' @seealso \code{\link{one_hot}}
#' @keywords one-hot
#'
#' @examples
#' one_hot(paste(sample(c('A','C','G','T'),100,replace=TRUE),collapse=''),10)
#'
#' @export
one_hot_vector <- function(dnastr,zeros_len, max_len) {
  a <- strsplit(as.character(dnastr),"")[[1]]
  buffer.left <- rep(0, ceiling((max_len+2*zeros_len-2-nchar(dnastr))/2))
  buffer.right <- rep(0, floor((max_len+2*zeros_len-2-nchar(dnastr))/2))
  oH <- c(buffer.left, as.numeric(a=="A"), buffer.right,
          buffer.left, as.numeric(a=="C"), buffer.right,
          buffer.left, as.numeric(a=="G"), buffer.right,
          buffer.left, as.numeric(a=="T"), buffer.right)
  return(oH)
}
