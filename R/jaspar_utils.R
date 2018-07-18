#' Read JASPAR file-format line.
#'
#' Reads in a single motif in JASPAR format. Slight modification of the original function available in the TFBSTools package (\code{https://github.com/ge11232002/TFBSTools}).
#'
#' @param text JASPAR format line of text
#'
#' @return A PFMatrix object.
#'
#' @author Matthew Ploenzke, \email{ploenzke@@g.harvard.edu}
#' @keywords read jaspar motif annotated
#'
#' @export
processJASPARText <- function(text){
  DNA_BASES <- c('A','C','G','T')
  ID <- sub("^>", "", strsplit(text[1], "\t")[[1]][1])
  name <- strsplit(text[1], "\t")[[1]][2]
  if(!identical(substr(text[2:5], 1, 1), DNA_BASES)){
    stop("The second to fifth lines of the file must start with",
         "`A`, `C`, `G`, `T`, respectively.")
  }
  profileMatrix <- do.call(rbind, strsplit(sub(" *]$", "",
                                               sub("^(A|C|G|T)  \\[ *", "",
                                                   text[2:5])), " +"))
  mode(profileMatrix) <- "integer"
  rownames(profileMatrix) <- DNA_BASES
  ans <- PFMatrix(ID=ID, name=name, profileMatrix=profileMatrix)
}

#' Read JASPAR file-format.
#'
#' Reads in JASPAR format motifs. Slight modification of the original function available in the TFBSTools package (\code{https://github.com/ge11232002/TFBSTools}).
#'
#' @param fn Filename containing JASPAR-formatted PFMs
#'
#' @return List of PFM motifs.
#'
#' @author Matthew Ploenzke, \email{ploenzke@@g.harvard.edu}
#' @keywords read jaspar motif annotated
#'
#' @export
my_readJASPARMatrix <- function(fn, type=c("individual", "all")){
  type <- match.arg(type)
  text <- readLines(fn)
  if(type == "individual"){
    if(length(text) != 5L){
      stop("The `individual` format is supposed to have 5 lines!")
    }
    ans <- processJASPARText(text)
  }else{
    if(length(text) %% 5 != 0L){
      stop("The `all` format is supposed to have a number of lines",
           "mutipled by 5!")
    }
    text2 <- split(text, rep(1:(length(text)/5), rep(5, length(text)/5)))
    ans <- lapply(text2, processJASPARText)
    ans.mot <- lapply(ans, function(i) {slot(i,'profileMatrix')})
    names(ans.mot) <- lapply(ans, function(i) {slot(i,'name')})
  }
  return(ans.mot)
}

#' Create PFMatrix object.
#'
#' Reads in a JASPAR-formatted motif and outputs as a PFMatrix class. Slight modification of the original function available in the TFBSTools package (\code{https://github.com/ge11232002/TFBSTools}).
#'
#' @param ID Identifier
#' @param name Motif name
#' @param matrixClass Type
#' @param strand Strand +/-
#' @param bg Background probabilities
#' @param tags Tags
#' @param profileMatrix Current profileMatrix to update
#'
#' @return ProfileMatrix PFMatrix.
#'
#' @author Matthew Ploenzke, \email{ploenzke@@g.harvard.edu}
#' @keywords pfmatrix matrix
#'
#' @importFrom methods new
#' @export
PFMatrix = function(ID="Unknown", name="Unknown", matrixClass="Unknown",
                    strand="+", bg=c(A=0.25, C=0.25, G=0.25, T=0.25),
                    tags=list(), profileMatrix=matrix()){
    mode(profileMatrix) <- "integer"
    new("PFMatrix", ID=ID, name=name, matrixClass=matrixClass,
        strand=strand, bg=bg,
        tags=tags,
        profileMatrix=profileMatrix)
}

setClass("PFMatrix",
         slots=c(ID="character",
                 name="character",
                 matrixClass="character",
                 strand="character",
                 bg="numeric",
                 tags="list",
                 profileMatrix="matrix")
         )

