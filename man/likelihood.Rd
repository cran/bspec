\name{likelihood}
\alias{dprior}
\alias{likelihood}
\alias{dposterior}
\alias{dprior.bspec}
\alias{likelihood.bspec}
\alias{dposterior.bspec}
\title{Prior, likelihood and posterior}
\description{
  Prior density, likelihood, and posterior density functions for the
  posterior distributions specified through a \code{bspec} object.
}
\usage{
dprior(x, ...)
likelihood(x, ...)
dposterior(x, ...)
\method{dprior}{bspec}(x, theta, two.sided=x$two.sided, log=FALSE, ...)
\method{likelihood}{bspec}(x, theta, two.sided=x$two.sided, log=FALSE, ...)
\method{dposterior}{bspec}(x, theta, two.sided=x$two.sided, log=FALSE, ...)
}
\arguments{
  \item{x}{a \code{bspec} object.}
  \item{theta}{a \code{numeric} vector of parameter values,
    corresponding to the Fourier frequencies in the \code{x$freq} element.}
  \item{two.sided}{a \code{logical} flag indicating whether the
    parameters \code{theta} correspond to the \emph{one-sided} or
    \emph{two-sided} spectrum.}
  \item{log}{a \code{logical} flag indicating whether to return
    logarithmic density (or likelihood) values.}
  \item{...}{}
}
\details{
  Prior and posterior are both \emph{scaled inverse
  \eqn{\chi^2}{chi-squared} distributions},
  and the likelihood is Normal.
}
\value{
  A \code{numeric} function value.
}
\references{Roever, C., Meyer, R., Christensen, N. (2008):
  Modelling coloured noise.
  (\emph{in preparation}).}
\author{Christian Roever, \email{bspec@web.de}}
\seealso{\code{\link[=bspec.default]{bspec}},
  \code{\link[=quantile.bspec]{quantile.bspec}},
  \code{\link[=expectation.bspec]{expectation}}}
\examples{
lhspec <- bspec(lh, priordf=1, priorscale=0.6)

# draw sample from posterior:
posteriorsample <- sample(lhspec)

# plot the sample:
plot(lhspec)
lines(lhspec$freq, posteriorsample, type="b", col="red")

# compute prior, likelihood, posterior:
print(c("prior"     = dprior(lhspec, posteriorsample),
        "likelihood"= likelihood(lhspec, posteriorsample),
        "posterior" = dposterior(lhspec, posteriorsample)))
}
\keyword{ts}
