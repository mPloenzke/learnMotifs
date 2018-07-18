#' One-hot encoding of vector of DNA nucleotide strings
#'
#' Inputs vector of strings composed of nucleotides and outputs a zero-padded array of encodings.
#'
#' @param data Character strings composed of nucleotides A, C, G, T.
#' @param zeros_len Integer value denoting number of zeros to pad vector with. Commonly the convolutional filter length - 1.
#'
#' @return Zero-padded array of one-hot encodings.
#'
#' @author Matthew Ploenzke, \email{ploenzke@@g.harvard.edu}
#' @seealso \code{\link{one_hot_vector}}
#' @keywords one-hot
#'
#' @examples
#' one_hot(rep(paste(sample(c('A','C','G','T'),100,replace=T),collapse=''),50),10)
#'
#' @export
one_hot <- function(data,zeros_len) {
  max_len <- max(nchar(data))
  tensor <- aperm(array(c(sapply(data,one_hot_vector,zeros_len=zeros_len)),dim=c((max_len+2*zeros_len-2),4,length(data),1)),perm=c(3,2,1,4))
  return(tensor)
}

#' One-hot encoding of DNA nucleotide string
#'
#' Inputs a string composed of nucleotides and outputs a zero-padded vector of encodings.
#'
#' @param dnastr Character string composed of nucleotides A, C, G, T.
#' @param zeros_len Integer value denoting number of zeros to pad vector with. Commonly the convolutional filter length - 1.
#'
#' @return Zero-padded vector of one-hot encodings. If the input is of length \eqn{L} and with padding of length Z the output is of length \eqn{4(L+2Z)}.
#'
#' @author Matthew Ploenzke, \email{ploenzke@@g.harvard.edu}
#' @seealso \code{\link{one_hot}}
#' @keywords one-hot
#'
#' @examples
#' one_hot(paste(sample(c('A','C','G','T'),100,replace=T),collapse=''),10)
#'
#' @export
one_hot_vector <- function(dnastr,zeros_len) {
  a <- strsplit(as.character(dnastr),"")[[1]]
  buffer.left <- rep(0, ceiling((opt$max_len+2*zeros_len-2-nchar(dnastr))/2))
  buffer.right <- rep(0, floor((opt$max_len+2*zeros_len-2-nchar(dnastr))/2))
  oH <- c(buffer.left, as.numeric(a=="A"), buffer.right,
          buffer.left, as.numeric(a=="C"), buffer.right,
          buffer.left, as.numeric(a=="G"), buffer.right,
          buffer.left, as.numeric(a=="T"), buffer.right)
  return(oH)
}
