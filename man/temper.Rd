\name{temper}
\alias{temper}
\alias{temper.bspec}
\title{Tempering of (posterior) distributions}
\description{
  Setting the tempering parameter of (\sQuote{tempered})
  \code{bspec} objects.
}
\usage{
  temper(x, ...)
  \method{temper}{bspec}(x, temperature = 2, likelihood.only = TRUE, ...)
}
\arguments{
  \item{x}{a \code{bspec} object.}
  \item{temperature}{a (positive) \sQuote{temperature} value.}
  \item{likelihood.only}{a \code{logical} flag indicating whether to
    apply the tempering to the \sQuote{complete} posterior density, or
    to the likelihood only (default).}
  \item{...}{}
}
\details{
  In the context of Markov chain Monte Carlo (MCMC) applications it is
  often desirable to apply \emph{tempering} to the distribution of
  interest, as it is supposed to make the distribution more easily
  tractable. Examples where tempering is utilised are \emph{simulated
  annealing}, \emph{parallel tempering} or \emph{evolutionary MCMC}
  algorithms. In the context of Bayesian inference, tempering may be
  done by specifying a \sQuote{temperature} \eqn{T}{T} and then
  manipulating the original posterior distribution
  \eqn{p(\theta|y)}{p(theta|y)} by applying an exponent
  \eqn{\frac{1}{T}}{1/T} either to the complete posterior distribution:
  \deqn{p_T(\theta) \propto p(\theta|y)^\frac{1}{T}%
        = (p(y|\theta)p(\theta))^\frac{1}{T}}{%
        p_T(theta)  ~  p(theta|y)^(1/T)  =  (p(y|theta)*p(theta))^(1/T)}
  or to the likelihood part only:
  \deqn{p_T(\theta) \propto p(y|\theta)^\frac{1}{T}p(\theta).}{%
        p_T(theta)  ~  p(theta|y)^(1/T) * p(theta).}
  In this context, where the posterior distribution is a product of
  \emph{scaled inverse \eqn{\chi^2}{chi-squared} distributions}, the
  tempered distributions in both cases turn out to be again of the same
  family, just with different parameters. For more details see also the
  references.
}
\value{
  An object of class \code{\link[=bspec.default]{bspec}}
  (see the help for the \code{\link[=bspec.default]{bspec}} function),
  but with an additional \code{temperature} element.
}
\references{Roever, C., Meyer, R., Christensen, N. (2008):
  Modelling coloured noise.
  (\emph{in preparation}).}
\author{Christian Roever, \email{bspec@web.de}}
\note{Tempering with the \code{likelihood.only} flag set to \code{FALSE}
  only works as long as the \code{temperature} is less than
  \code{min((x$df+2)/2)}.
}
\seealso{\code{\link[=temperature.bspec]{temperature}},
  \code{\link[=bspec.default]{bspec}}}
\examples{
lhspec <- bspec(lh, priorscale=0.6, priordf=1)

# details of the regular posterior distribution:
str(lhspec)

# details of the tempered distribution
# (note the differing scale and degrees-of-freedom):
str(temper(lhspec, 1.23))
}
\keyword{ts}
