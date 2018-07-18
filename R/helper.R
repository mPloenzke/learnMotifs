#' Initialize output directory.
#'
#' Creates the initial output directory, optionally backing up any files which may be overwritten.
#'
#' @param opt Options list.
#'
#' @return NULL. Output directory will be initialized.
#'
#' @author Matthew Ploenzke, \email{ploenzke@@g.harvard.edu}
#' @seealso \code{\link{make_epoch_dir}}
#' @keywords log directory
#'
#' @examples
#' setup_log_dir(list(log_dir='output_directory',cache_Old=F))
#'
#' @export
setup_log_dir <- function(opt) {
  folds <- list.dirs(path=opt$log_dir)[!grepl(pattern='past_runs',list.dirs(path=opt$log_dir)) &
                                         grepl(pattern='epoch',list.dirs(path=opt$log_dir))]
  if (file.exists(opt$log_dir)) {
    if (length(folds) > 0 & opt$cache_Old) {
      if (!file.exists(file.path(opt$log_dir,'past_runs'))) {
        dir.create(file.path(opt$log_dir,'past_runs'))
      }
      tempdir <- format(Sys.time(),format="%Y-%m-%d_%H:%M:%S")
      dir.create(file.path(opt$log_dir,'past_runs',tempdir))
      for (fi in folds) {
        file.copy(from=fi,to=file.path(opt$log_dir,'past_runs',tempdir),recursive = T)
        unlink(fi,recursive=T)
      }
      file.copy(from=list.files(opt$log_dir,pattern='.csv',full.names=T),to=file.path(opt$log_dir,'past_runs',tempdir))
      file.remove(list.files(opt$log_dir,pattern='.csv',full.names=T))
    } else if (length(folds)>0 & !opt$cache_Old) {for (fi in folds) {unlink(fi,recursive = T)}}
  } else {
    dir.create(opt$log_dir)
  }
  capture.output(opt, file = file.path(opt$log_dir,'model_Opts.csv'))
}

#' Create output directory for current training epoch.
#'
#' Creates the epoch output directory.
#'
#' @param log_dir Initial log directory.
#' @param epoch Current training epoch.
#' @param wts Current weights to output.
#'
#' @return NULL. Epoch directory will be initialized.
#'
#' @author Matthew Ploenzke, \email{ploenzke@@g.harvard.edu}
#' @seealso \code{\link{setup_log_dir}}
#' @keywords epoch log directory
#'
#' @export
make_epoch_dir <- function(log_dir, epoch, wts) {
  epoch_dir <- paste(log_dir,'/epoch_',epoch,sep='')
  dir.create(epoch_dir)
  capture.output(wts, file = file.path(epoch_dir,'current_weights.csv'))
  return(epoch_dir)
}
#' Convert dna strings to PFM
#'
#' Converts nucleotide strings to position frequency matrices. Function available in the TFBSTools package (\code{https://github.com/ge11232002/TFBSTools}).
#'
#' @param x nucleotide strings.
#'
#' @return A matrix.
#'
#' @author Matthew Ploenzke, \email{ploenzke@@g.harvard.edu}
#' @keywords iupac biostrings tfbstools
#'
#' @export
IUPAC2Matrix <- function(x){
  x <- as.character(x)
  x <- strsplit(x, "")[[1]]
  if(!all(x %in% names(IUPAC_CODE_MAP))){
    stop("All characters must be in IUPAC_CODE_MAP!")
  }
  ans <- matrix(0L, nrow=4, ncol=length(x),
                dimnames=list(c("A", "C", "G", "T")))
  for(i in 1:length(x)){
    dnaCharacters <- strsplit(IUPAC_CODE_MAP[x[i]], "")[[1]]
    ans[dnaCharacters, i] <- 1L
  }
  return(ans)
}
