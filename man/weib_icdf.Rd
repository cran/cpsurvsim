% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/weibullfns.R
\name{weib_icdf}
\alias{weib_icdf}
\title{Inverse CDF value generation for the Weibull distribution}
\usage{
weib_icdf(n, gamma, theta)
}
\arguments{
\item{n}{Number of output Weibull values}

\item{gamma}{Shape parameter \eqn{\gamma}}

\item{theta}{Scale parameter \eqn{\theta}}
}
\value{
Output is a value or vector of values
from the Weibull distribution.
}
\description{
\code{weib_icdf} returns a value from the Weibull distribution by
using the inverse CDF.
}
\details{
This function uses the Weibull density of the form
\deqn{f(t)=\theta t^(\gamma - 1)exp(-\theta/\gamma t^(\gamma))}
to get the inverse CDF
\deqn{F^(-1)(u)=(-\gamma/\theta log(1-u))^(1/\gamma)} where \eqn{u}
is a uniform random variable. It can be implemented directly and is
also called by the function \code{\link{weib_memsim}}.
}
\examples{
simdta <- weib_icdf(n = 10, theta = 0.05, gamma = 2)

}
