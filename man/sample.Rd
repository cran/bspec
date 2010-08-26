\name{sample.bspec}
\alias{sample.bspec}
\alias{sample.default}
\alias{sample}
\title{Posterior sampling}
\description{
  Function to generate samples from the spectrum's posterior
  distribution specified through the supplied \code{bspec}
  object argument.
}
\usage{
  \method{sample}{bspec}(x, size = 1, two.sided = x$two.sided, ...)
}
\arguments{
  \item{x}{a \code{bspec} object.}
  \item{size}{the sample size.}
  \item{two.sided}{a \code{logical} flag indicating whether the drawn
    samples are supposed to correspond to the \emph{one-sided} or
    \emph{two-sided} spectrum.}
  \item{...}{currently unused.}
}
\details{
  The posterior distribution is a product of independent \emph{scaled inverse
    \eqn{\chi^2}{chi-squared} distributions}.
}
\value{
  A (\code{numerical}) vector of samples from the posterior distribution
  of the spectral parameters, of the same length as and corresponding to
  the \code{$freq} element of the supplied \code{bspec} object.
}
\author{Christian Roever, \email{bspec@web.de}}
\seealso{\code{\link[=bspec.default]{bspec}},
  \code{\link[=ppsample.bspec]{ppsample}}}
\examples{
# determine spectrum's posterior distribution:
lhspec <- bspec(lh)

# plot 95 percent central intervals and medians:
plot(lhspec)

# draw and plot two samples from posterior distribution:
lines(lhspec$freq, sample(lhspec), type="b", pch=20, col="red")
lines(lhspec$freq, sample(lhspec), type="b", pch=20, col="green")
}
\keyword{ts}
