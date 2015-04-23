\name{one.sided}
\alias{one.sided}
\alias{two.sided}
\alias{one.sided.bspec}
\alias{two.sided.bspec}
\title{Conversion between one- and two-sided spectra}
\description{
  Functions to convert between one- and two-sided
  \code{bspec} objects.
}
\usage{
one.sided(x, ...)
two.sided(x, ...)
\method{one.sided}{bspec}(x, ...)
\method{two.sided}{bspec}(x, ...)
}
\arguments{
  \item{x}{a \code{bspec} object.}
  \item{...}{currently unused.}
}
\details{
  The conversion only means that the \code{$two.sided} element of the
  returned \code{bspec} object is set correspondingly, as internally always
  the same (one-sided) spectrum is used.
}
\value{
  A \code{\link[=bspec.default]{bspec}} object
  (see the help for the \code{\link[=bspec.default]{bspec}} function).
}
\author{Christian Roever, \email{christian.roever@med.uni-goettingen.de}}
\seealso{\code{\link[=bspec.default]{bspec}}}
\examples{
lhspec <- bspec(lh)

# compare distributions visually:
par(mfrow=c(2,1))
  plot(lhspec)
  plot(two.sided(lhspec))
par(mfrow=c(1,1))

# ...and numerically:
print(cbind("frequency"=lhspec$freq,
            "median-1sided"=quantile(lhspec,0.5),
            "median-2sided"=quantile(two.sided(lhspec),0.5)))
}
\keyword{ts}
