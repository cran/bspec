\name{ppsample}
\alias{ppsample}
\alias{ppsample.bspec}
\title{Posterior predictive sampling}
\description{
  Draws a sample from the posterior predictive distribution specified by
  the supplied \code{bspec} object.
}
\usage{
ppsample(x, ...)
\method{ppsample}{bspec}(x, start=x$start, ...)
}
\arguments{
  \item{x}{a \code{bspec} object specifying the posterior
    distribution from which to sample.}
  \item{start}{the start time of the resulting time series.}
  \item{...}{currently unused.}
}
\value{
  A time series (\code{ts}) object of the same kind (with respect to
  sampling rate and sample size) as the data the posterior distribution
  is based on.
}
\author{Christian Roever, \email{christian.roever@med.uni-goettingen.de}}
\seealso{\code{\link[=bspec.default]{bspec}}, \code{\link[=sample.bspec]{sample.bspec}}}
\examples{

par(mfrow=c(2,1))
plot(lh, main="'lh' data")
plot(ppsample(bspec(lh)), main="posterior predictive sample")
par(mfrow=c(1,1))
}
\keyword{ts}
