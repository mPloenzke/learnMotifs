% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/format_annotated_motifs.R
\name{pfm_to_pwm}
\alias{pfm_to_pwm}
\title{Convert PFM to PWM.}
\usage{
pfm_to_pwm(motif.raw, motif.pos = NULL, mmax, bg = c(0.25, 0.25, 0.25,
  0.25), pseudoprob = 0.001)
}
\arguments{
\item{motif.raw}{List of annotated motifs as PFMs.}

\item{motif.pos}{Vector of indices to extract motifs. If \code{NULL} all motifs will be selected.}

\item{mmax}{Maximum motif length, found via \code{\link{find_max_length}}.}

\item{bg}{Vector of background nucleotide probabilities (A, C, G, T)}

\item{pseudoprob}{Value to add to zero-valued position probability matrix entries prior to converting to PWM.}
}
\value{
4-D array containing PWMs of the extracted motifs.
}
\description{
Convert position frequency matrix into a position weight matrix.
}
\seealso{
\code{\link{find_max_length}}
}
\author{
Matthew Ploenzke, \email{ploenzke@g.harvard.edu}
}
\keyword{motif}
\keyword{pfm}
\keyword{position}
\keyword{pwm}
\keyword{weight}
