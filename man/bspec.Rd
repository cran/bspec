\name{bspec}
\alias{bspec}
\alias{bspec.default}
\alias{print.bspec}
\alias{plot.bspec}
\alias{is.bspec}
\title{Computing the spectrum's posterior distribution}
\description{
  Derives the posterior distribution of the spectrum
  of one or several time series,
  based on data and prior specifications.
}
\usage{
  bspec(x, ...)
  \method{bspec}{default}(x, priorscale=1, priordf=0, intercept=TRUE,
    two.sided=FALSE, ...)
}
\arguments{
  \item{x}{a time series object of the data to be analysed.
    May be a univariate (\code{\link{ts}} object) or multivariate
    (\code{\link[=ts]{mts}} object) time series.}
  \item{priorscale}{\emph{either} a \code{numeric} vector giving the
    scale parameters of the spectrum's prior distribution; recycled if
    of length 1.

    \emph{Or} a \code{function} of frequency.}
  \item{priordf}{\emph{either} a \code{numeric} vector giving the
    degrees-of-freedom parameters of the spectrum's prior distribution;
    recycled if of length 1.

    \emph{Or} a \code{function} of frequency.}
  \item{intercept}{a \code{logical} flag indicating whether to include
    the \sQuote{intercept} (zero frequency) term.}
  \item{two.sided}{a \code{logical} flag indicating whether to refer to
    a one-sided or a two-sided spectrum. In particular affects the
    interpretation of the prior scale parameters, and sets the default
    for some methods applied to the resulting \code{bspec} object via
    its \code{two.sided} element.}
  \item{...}{}
}
\details{
  Based on the assumptions of a zero mean and a finite spectrum,
  the posterior distribution of the (discrete) spectrum is derived.
  The data are modeled using the \emph{Maximum Entropy} (Normal)
  distribution for the above constraints, and based on the
  prior information about the spectrum specified in terms of the
  (conjugate) \emph{scaled inverse \eqn{\chi^2}{chi-squared} distribution}.

  For more details, see the references.
}
\value{
  A list of class \code{bspec} containing the following elements:
  \item{freq}{a \code{numeric} vector giving the
    (Fourier-) frequencies that the spectral parameters
    correspond to.}
  \item{scale}{a \code{numeric} vector giving the scale
    parameters of the posterior distributions of the spectral
    parameters corresponding to the above frequencies.
    These -internally- always correspond to the \emph{one-sided}
    spectrum, regardless of the \code{two.sided} flag (see below).}
  \item{df}{a \code{numeric} vector giving the
    degrees-of-freedom parameters of the posterior distributions of
    the spectral parameters corresponding to the above frequencies.}
  \item{priorscale}{a \code{numeric} vector giving the
    prior scale parameters.}
  \item{priordf}{a \code{numeric} vector giving the
    prior degrees-of-freedom parameters.}
  \item{datassq}{a \code{numeric} vector giving the
    sum-of-squares contributed by the data.}
  \item{datadf}{a \code{numeric} vector giving the
    degrees-of-freedom contributed by the data.}
  \item{N}{the sample size of the original time series.}
  \item{deltat}{the sampling interval of the original time series.}
  \item{deltaf}{the frequency interval of the Fourier-transformed time series.}
  \item{start}{the time of the first observation in the original time series.}
  \item{call}{an object of class \code{call} giving the
    function call that generated the \code{bspec} object.}
  \item{two.sided}{a \code{logical} flag indicating
    whether the spectrum is to be interpreted as one-sided or two-sided.}
}
\references{Roever, C., Meyer, R., Christensen, N. (2008):
  Modelling coloured noise.
  (\emph{in preparation}).}
\author{Christian Roever, \email{bspec@web.de}}
\seealso{\code{\link[=expectation.bspec]{expectation}},
  \code{\link[=quantile.bspec]{quantile.bspec}},
  \code{\link[=sample.bspec]{sample.bspec}},
  \code{\link[=ppsample.bspec]{ppsample}},
  \code{\link[=acf.bspec]{acf.bspec}},
  \code{\link{spectrum}}
}
\examples{
# determine spectrum's posterior distribution
# (for noninformative prior):
lhspec <- bspec(lh)
print(lhspec)

# show some more details:
str(lhspec)

# plot 95 percent central intervals and medians:
plot(lhspec)

# draw and plot a sample from posterior distribution:
lines(lhspec$freq, sample(lhspec), type="b", pch=20)

########

# compare the default outputs of "bspec()" and "spectrum()":
bspec1    <- bspec(lh)
spectrum1 <- spectrum(lh, plot=FALSE)
plot(bspec1) 
lines(spectrum1$freq, spectrum1$spec, col="blue")
# (note -among others- the factor 2 difference)

# match the outputs:
# Need to suppress  tapering, padding and de-trending
# (see help for "spec.pgram()"):
spectrum2 <- spectrum(lh, taper=0, fast=FALSE, detrend=FALSE, plot=FALSE)
# Need to drop intercept (zero frequency) term:
bspec2    <- bspec(lh, intercept=FALSE)
# plot the "spectrum()" output:
plot(spectrum2)
# draw the "bspec()" scale parameters, adjusted
# by the corresponding degrees-of-freedom,
# so they correspond to one-sided spectrum:
lines(bspec2$freq, bspec2$scale/bspec2$datadf,
      type="b", col="green", lty="dashed")



}
\keyword{ts}
