#' Run Tomtom to match motifs.
#'
#' Checks estimated filter motifs against JASPAR database motifs and returns significant motif hits. Must have the Memesuit installed \url{http://meme-suite.org/}.
#'
#' @param W_conv_list List of motifs to run through Tomtom.
#' @param pfm Logical indicating whether position frequency matrices are provided (\code{TRUE}) or information content matrices (\code{FALSE}).
#'
#' @return \code{tomtom_out} directory containing Tomtom output.
#'
#' @author Matthew Ploenzke, \email{ploenzke@@g.harvard.edu}
#' @keywords tomtom motif jaspar
#'
#' @export
run_tomtom <- function(W_conv_list, pfm=F) {
  if (!file.exists(opt$JASPAR)) {
    download.file('http://jaspar.genereg.net/download/CORE/JASPAR2018_CORE_vertebrates_non-redundant_pfms_meme.txt',opt$JASPAR)
  }
  if (pfm) {
    system(paste('sites2meme ',epoch_dir,' > ',paste(epoch_dir,'PFM_meme_format.txt',sep="/"),sep=''))
    fl <- 'PFM_meme_format.txt'
  } else {
    W_conv_list <- lapply(W_conv_list, function(i) {
      ppm <- i/matrix(colSums(i),nrow=4,ncol=ncol(i),byrow=T)
      ppm <- ppm[,!is.nan(colSums(ppm)),drop=F]
      if (ncol(ppm)<1) {return(NULL)}
      return(ppm)
    })
    invisible(lapply(seq_along(W_conv_list),function(i) {
      if (!is.null(W_conv_list[[i]])) {
        write(gsub(' ','_',names(W_conv_list)[[i]]),file=paste(epoch_dir,"PPM_filters.txt",sep="/"),append=TRUE)
        row.names(W_conv_list[[i]]) <- paste(row.names(W_conv_list[[i]]),':',sep='')
        write.table(W_conv_list[[i]],file=paste(epoch_dir,"PPM_filters.txt",sep="/"),
                    append=TRUE,row.names=T,col.names=F,sep='\t',quote=F)
        write('\n',file=paste(epoch_dir,"PPM_filters.txt",sep="/"),append=T)
      }
    }))
    system(paste('uniprobe2meme ',paste(epoch_dir,"PPM_filters.txt",sep="/"),' > ',paste(epoch_dir,'PPM_meme_format.txt',sep='/'),sep=''))
    fl <- 'PPM_meme_format.txt'
  }
  system(paste('tomtom -oc ',paste(epoch_dir,'/tomtom_out ',sep=''),'-thresh ',opt$motif_qval_threshold,' ',paste(epoch_dir,fl,sep='/'),opt$JASPAR,sep=''))
}
