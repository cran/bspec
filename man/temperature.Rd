\name{temperature}
\alias{temperature}
\alias{temperature.bspec}
\title{Querying the tempering parameter}
\description{
  Retrieving the \dQuote{temperature} parameter of (\sQuote{tempered})
  \code{bspec} objects
}
\usage{
  temperature(x, ...)
  \method{temperature}{bspec}(x, ...)
}
\arguments{
  \item{x}{a \code{bspec} object.}
  \item{...}{currently unused.}
}
\value{
  The (\code{numeric}) value of the \code{temperature} element of the
  supplied \code{\link[=bspec.default]{bspec}} object, if present, and 1
  otherwise.
}
\references{Roever, C. (2007):
  \href{http://hdl.handle.net/2292/2356}{Bayesian inference on
    astrophysical binary inspirals based on gravitational-wave
    measurements}.
  PhD thesis, Department of Statistics, The University of Auckland,
  New Zealand.}
\author{Christian Roever, \email{bspec@web.de}}
\seealso{\code{\link[=temper.bspec]{temper}},
  \code{\link[=bspec.default]{bspec}}}
\examples{
lhspec1 <- bspec(lh)
lhspec2 <- temper(lhspec1, 1.23)

print(lhspec2$temperature)
print(lhspec1$temperature)

print(temperature(lhspec2))
print(temperature(lhspec1))
}
\keyword{ts}
