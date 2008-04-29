\name{acf.bspec}
\alias{acf.bspec}
\alias{acf.default}
\alias{acf}
\alias{plot.bspecACF}
\alias{print.bspecACF}
\alias{is.bspecACF}
\title{Posterior autocovariances}
\description{Deriving (posterior) autocovariances or autocorrelations
  from the spectrum's posterior distribution.
}
\usage{
\method{acf}{bspec}(x, spec = NULL,
   type = c("covariance", "correlation"),
   two.sided = x$two.sided, ...)
}
\arguments{
  \item{x}{a \code{bspec} object.}
  \item{spec}{(optional) a \code{numeric} vector giving \emph{fixed}
    values of the spectral parameters (e.g. derived by the
    \code{\link[=sample.bspec]{sample}} function) for which the
    autocovariances then are deterministic.}
  \item{type}{a \code{character} string specifying the desired type of
    output.}
  \item{two.sided}{a \code{logical} flag indicating whether the
    \code{spec} values are to be interpreted as \emph{one-sided} or
    \emph{two-sided}.}
  \item{...}{}
}
\details{
  If \code{spec} is supplied, the autocovariance (or autocorrelation)
  function corresponding to that specific spectrum will be returned.
  As this is a completely deterministic relationship, the
  \dQuote{\code{stderr}} slot of the result will be zero in this case.

  If \code{spec} is \emph{not} supplied, the \emph{(posterior) expected
  autocovariance} is returned in the \dQuote{\code{acf}} element, and its
  \emph{(posterior) standard deviation} is returned in the
  \dQuote{\code{stderr}} element.
  The posterior expectation of the autocovariance is only finite if
  \emph{all (!)} posterior degrees-of-freedom parameters in the
  \code{bspec} object are \eqn{>2}{>2}. The posterior
  variance (and with that the \code{stderr} element) is only finite if all
  these are \eqn{>4}{>4}.

  Autocorrelations are only returned if \code{spec} is supplied.
}
\value{
  A list of class \code{bspecACF} containing the following components:
  \item{lag}{a \code{numeric} vector giving the lags
    corresponding to the (discrete) autocovariance / autocorrelation
    values.}
  \item{acf}{a \code{numeric} vector giving the values of the
    autocovariance / autocorrelation function correponding to the
    above lags.}
  \item{stderr}{a \code{numeric} vector giving the standard
    errors (posterior standard deviations) of the above autocovariance
    values.}
  \item{type}{a \code{character} string giving the
    nature of the above \code{acf} element: either \code{"covariance"}
    or \code{"correlation"}.}
  \item{N}{an \code{integer} giving the sample
    size of the original time series.}
  \item{bspec}{a \code{character} string giving the
    name of the \code{bspec} object the \code{bspecACF} object was
    generated from.}
}
\references{Roever, C., Meyer, R., Christensen, N. (2008):
  Modelling coloured noise.
  Arxiv preprint 0804.3853 [stat.ME],
  \url{http://arxiv.org/abs/0804.3853},
  (\emph{submitted for publication}).}
\author{Christian Roever, \email{bspec@web.de}}
\note{
  (Posterior) expectation and standard deviation of the spectrum may in
  many cases not be finite (see above).
  Autocorrelations are only returned if \code{spec} is supplied.
}
\seealso{\code{\link[=bspec.default]{bspec}},
  \code{\link[=expectation.bspecACF]{expectation}},
  \code{\link[=sample.bspec]{sample.bspec}},
  \code{\link{acf}}
}
\examples{
lhspec1 <- bspec(lh)

# without any prior specifications,
# autocovariances are not finite:
print(acf(lhspec1))
str(acf(lhspec1))

# for given values of the spectral parameters,
# the autocovariances are fixed:
str(acf(lhspec1, spec=sample(lhspec1)))

# for all the prior degrees-of-freedom greater than one,
# the expected autocovariance is finite, its variance isn't:
lhspec2 <- bspec(lh, priordf=2, priorscale=0.6, intercept=FALSE)
print(acf(lhspec2))
str(acf(lhspec2))
plot(acf(lhspec2))
}
\keyword{ts}
