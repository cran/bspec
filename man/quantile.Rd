\name{quantile.bspec}
\alias{quantile.bspec}
\title{Quantiles of the posterior spectrum}
\description{
  Function to compute quantiles of the spectrum's posterior
  distribution specified through the supplied \code{bspec}
  object argument.
}
\usage{
\method{quantile}{bspec}(x, probs = c(0.025, 0.5, 0.975),
  two.sided = x$two.sided, ...)
}
\arguments{
  \item{x}{a \code{bspec} object.}
  \item{probs}{a \code{numerical} vector of probabilities.}
  \item{two.sided}{a \code{logical} flag indicating whether quantiles
    are supposed to correspond to the \emph{one-sided} or
    \emph{two-sided} spectrum.}
  \item{...}{currently unused.}
}
\details{
  The posterior distribution is a product of independent \emph{scaled inverse
    \eqn{\chi^2}{chi-squared} distributions}.
}
\value{
  A matrix with columns corresponding to elements of \code{probs}, and
  rows corresponding to the Fourier frequencies \code{x$freq}.
  If \code{probs} is of length 1, a vector is returned instead.
}
\author{Christian Roever, \email{bspec@web.de}}
\seealso{\code{\link[=bspec.default]{bspec}}, \code{\link[stats]{quantile}}}
\examples{
lhspec <- bspec(lh)

# posterior medians:
print(cbind("frequency"=lhspec$freq,
            "median"=quantile(lhspec, 0.5)))
}
\keyword{ts}
