% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/jaspar_utils.R
\name{processJASPARText}
\alias{processJASPARText}
\title{Read JASPAR file-format line.}
\usage{
processJASPARText(text)
}
\arguments{
\item{text}{JASPAR format line of text}
}
\value{
A PFMatrix object.
}
\description{
Reads in a single motif in JASPAR format. Slight modification of the original function available in the TFBSTools package (\code{https://github.com/ge11232002/TFBSTools}).
}
\author{
Matthew Ploenzke, \email{ploenzke@g.harvard.edu}
}
\keyword{annotated}
\keyword{jaspar}
\keyword{motif}
\keyword{read}
