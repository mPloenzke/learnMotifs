#' Nucleotide sequences sampled from known motifs.
#'
#' A dataset containing nucleotide sequences of length 200 basepairs composed of nucleotides A, C, G, T.
#' Each sequence contains between 0-3 ocurrences of the MYC, CTCF, or IRF motifs.
#'
#' @format A data frame with 5000 rows and 3 variables:
#' \describe{
#'   \item{seqname}{Sequence identifier}
#'   \item{sequence}{Nucleotide sequences}
#'   \item{embeddings}{Locations and descriptions of embedded motifs.}
#' }
#' @source \url{https://github.com/kundajelab/simdna}
"MYC_CTCF_IRF"
#' Nucleotide sequences sampled from random background.
#'
#' A dataset containing nucleotide sequences of length 200 basepairs composed of nucleotides A, C, G, T.
#' Each sequence is sampled from random genome background with 60% GC-content.
#'
#' @format A data frame with 5000 rows and 3 variables:
#' \describe{
#'   \item{seqname}{Sequence identifier}
#'   \item{sequence}{Nucleotide sequences}
#'   \item{embeddings}{Locations and descriptions of embedded motifs.}
#' }
#' @source \url{https://github.com/kundajelab/simdna}
"EmptyBackground"
