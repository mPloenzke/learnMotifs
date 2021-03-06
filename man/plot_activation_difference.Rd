% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot_definitions.R
\name{plot_activation_difference}
\alias{plot_activation_difference}
\title{Plot effect size versus activation difference.}
\usage{
plot_activation_difference(betas, combo = F, fl = NULL)
}
\arguments{
\item{betas}{Tibble containing effect size information, filter name, mean activation difference, and information.}

\item{combo}{Whether both de novo and annotated motifs are included.}

\item{fl}{Filename to save plots to. If not provided, a plot object will be returned.}
}
\value{
A ggplot object.
}
\description{
Plots the estimated motif effect size (beta) against the mean between-class activation difference.
}
\examples{
a <- runif(10,0,1)
b <- rnorm(10)
c <- paste(1:10)
d <- rnorm(10,1)
mat <- dplyr::tibble(Difference=a,Effect=b,Filter=c,Information=d)
plot_activation_difference(mat)

}
\author{
Matthew Ploenzke, \email{ploenzke@g.harvard.edu}
}
\keyword{beta}
\keyword{betas}
\keyword{plot}
