% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/constraint_definitions.R
\name{pwm_constraint}
\alias{pwm_constraint}
\title{Constrains filter weights within position weight bounds.}
\usage{
pwm_constraint(w)
}
\arguments{
\item{w}{Filter weights.}
}
\value{
Constraint.
}
\description{
Constrains the individual filters weights less than or equal to 2.
}
\author{
Matthew Ploenzke, \email{ploenzke@g.harvard.edu}
}
\keyword{constraint}
\keyword{position}
\keyword{pwm}
\keyword{weight}
