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
\references{Roever, C., Meyer, R., Christensen, N. (2008):
  Modelling coloured noise.
  Arxiv preprint 0804.3853 [stat.ME],
  \url{http://arxiv.org/abs/0804.3853},
  (\emph{submitted for publication}).}
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
