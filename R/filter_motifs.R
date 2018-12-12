#' Format filters for motif plotting.
#'
#' Convert the information content motifs to a list with named rows for plotting with \code{\link{plot_motifs}}.
#'
#' @param conv_filters 4-dimensional array containing filter weights.
#' @param type String declaring motif representation type (either `IGM` or `PWM`).
#' @param bg Vector of background nucleotide probabilities (A, C, G, T)
#'
#' @return List with named motifs.
#'
#' @author Matthew Ploenzke, \email{ploenzke@@g.harvard.edu}
#' @keywords motif filter format
#'
#' @export
format_filter_motifs <- function(conv_filters, type='IGM', bg=c(.25,.25,.25,.25)) {
  W_conv_list <- list()
  if (type=='PWM') {
    for (filt in 1:dim(conv_filters)[4]) {
      i <- conv_filters[,,,filt,drop=FALSE]
      bg.array <- array(matrix(bg,nrow=4,ncol=ncol(i),byrow=F),dim=dim(i))
      ppm <- (2^i)*bg.array
      ppm <- ppm/array(matrix(colSums(ppm),nrow=4,ncol=ncol(i),byrow=T),dim=dim(i))
      R_i = pmax(apply(ppm, 2, function(col) sum(col * log2(col/bg),na.rm=T)),0)
      conv_filters[,,,filt] <- ppm*array(matrix(R_i,nrow=4,ncol=ncol(ppm),byrow=T),dim=dim(ppm))
    }
  }
  for (filt in 1:dim(conv_filters)[4]) {
    if (all(conv_filters[,,,filt]==0)) {conv_filters[,,,filt] <- conv_filters[,,,filt]+1e-10}
    W_conv_list[[paste('Filter',filt,sep=' ')]] <- matrix(conv_filters[,,,filt],nrow=4,dimnames=list(c('A', 'C', 'G', 'T')))
  }
  return(W_conv_list)
}
#' Find maximally-activating sequences for motif plotting.
#'
#' Find the sequeneces which maximally activate (maximially cross-correlated with) the filters and format to a list with named rows for plotting with \code{\link{plot_motifs}}.
#'
#' @param activations 4-D array of activation values per filter.
#' @param test_seqs Data to extract nucleotide sequences from. Must be a character array.
#' @param filter_len Length of the sequence motif to extract; the length of the convolutional filters.
#' @param method Either \code{'alipinahi'} (extract the sequences with activation greater than \code{threshold*max_activation}),
#' \code{'max-pooled'} (extract the sequences which gave rise to the max-pooled value), or
#' \code{'quantile'} (extract the sequences with activation greater than \code{quantile(activations,probs=threshold)}.
#' @param threshold Threshold above which to extract sequences.
#'
#' @return List with named motifs.
#'
#' @author Matthew Ploenzke, \email{ploenzke@@g.harvard.edu}
#' @keywords motif filter format
#'
#' @importFrom stats quantile
#' @importFrom stringi stri_dup
#' @export
format_activation_motifs <- function(activations, test_seqs, filter_len, method, threshold) {
  all_dna_strings <- list()
  for (filt in 1:dim(activations)[4]) {
    act <- activations[,,,filt]
    if (method=='alipinahi') {
      act[act < threshold*max(act)] <- 0
    } else if (method=='max-pooled') {
      maxidxs <- apply(act,1,which.max)
      pos.values <- cbind(1:nrow(act),maxidxs)
      act <- do.call(rbind,lapply(1:nrow(act),function(i) {
        act[i,-pos.values[i,2]] <- 0
        act[i,]
      }))
    } else if (method=='quantile') {
      act[act < quantile(act[act!=0],probs=.99,na.rm=TRUE)] <- 0
    }
    pos.values <- which(act>0,arr.ind=TRUE)
    pfm <- matrix(0,nrow=4,ncol=filter_len)
    if (nrow(pos.values)>0) {
      seqs <- lapply(1:nrow(pos.values), function(seqq) {
        temp_str <- test_seqs[pos.values[seqq,1]]
        temp_idx <- pos.values[seqq,2]
        temp_idx <- ifelse(temp_idx-filter_len<1,temp_idx-(temp_idx-filter_len),temp_idx)
        strseq <- substr(temp_str, temp_idx-(filter_len-1), temp_idx)
        strseq <- ifelse(nchar(strseq)<filter_len,paste(strseq,stri_dup("N",filter_len-nchar(strseq)),sep=''),strseq)
        return(strseq)
      })
      dna_strings <- do.call(c,seqs)
      all_dna_strings[[paste('Motif',filt,sep=' ')]] <- dna_strings
      for (seqq in 1:length(dna_strings)) {
        pfm <- pfm + IUPAC2Matrix(dna_strings[seqq])
      }
      colnames(pfm) <- paste("POS",1:filter_len,sep="")
    }
  }
  return(all_dna_strings)
}
