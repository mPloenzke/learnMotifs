#' Read HOCOMOCO file-format line.
#'
#' Reads in a single motif in HOCOMOCO motif in JASPAR format. Slight modification of the original function available in the TFBSTools package (\code{https://github.com/ge11232002/TFBSTools}).
#'
#' @param text JASPAR format line of text
#'
#' @return A PFMatrix object.
#'
#' @author Matthew Ploenzke, \email{ploenzke@@g.harvard.edu}
#' @keywords read motif annotated hocomoco
#'
#' @export
processHOCOMOCOText <- function(text){
  DNA_BASES <- c('A','C','G','T')
  ID <- sub("^>", "", strsplit(text[1], "\t")[[1]][1])
  name <- ID
  profileMatrix <- do.call(rbind,strsplit(text[2:5], "\t"))
  mode(profileMatrix) <- "integer"
  rownames(profileMatrix) <- DNA_BASES
  ans <- PFMatrix(ID=ID, name=name, profileMatrix=profileMatrix)
}

#' Read HOCOMOCO file in JASPAR-format.
#'
#' Reads in HOCOMOCO motifs. Slight modification of the original function available in the TFBSTools package (\code{https://github.com/ge11232002/TFBSTools}).
#'
#' @param fn Filename containing JASPAR-formatted HOCOMOCO PFMs.
#' @param type Does the file contain \code{'individual'} PFMs or a list (\code{'all'})?
#'
#' @return List of PFM motifs.
#'
#' @author Matthew Ploenzke, \email{ploenzke@@g.harvard.edu}
#' @keywords read hocomoco motif annotated
#'
#' @export

my_readHOCOMOCOMatrix <- function(fn, type=c("individual", "all")){
  type <- match.arg(type)
  text <- readLines(fn)
  if(type == "individual"){
    if(length(text) != 5L){
      stop("The `individual` format is supposed to have 5 lines!")
    }
    ans <- processHOCOMOCOText(text)
  }else{
    if(length(text) %% 5 != 0L){
      stop("The `all` format is supposed to have a number of lines",
           "mutipled by 5!")
    }
    text2 <- split(text, rep(1:(length(text)/5), rep(5, length(text)/5)))
    ans <- lapply(text2, processHOCOMOCOText)
    ans.mot <- lapply(ans, function(i) {slot(i,'profileMatrix')})
    names(ans.mot) <- lapply(ans, function(i) {slot(i,'name')})
  }
  return(ans)
}
