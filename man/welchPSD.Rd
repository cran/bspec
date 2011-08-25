\name{welchPSD}
\alias{welchPSD}
\title{Power spectral density estimation using Welch's method.}
\description{
  Estimates a time series' power spectral density using \emph{Welch's
    method}, i.e., by subdividing the data into segments, computing
  spectra for each, and averaging over the results.
}
\usage{
welchPSD(x, seglength, two.sided = FALSE, windowfun = tukeywindow,
         method = c("mean", "median"), windowingPsdCorrection = TRUE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{a time series (\code{\link{ts}} object).}
  \item{seglength}{the length of the subsegments to be used (in units of
    time relative to \code{x}).}
  \item{two.sided}{a \code{logical} flag indicating whether the result is
    supposed to be one-sided (default) or two-sided.}
  \item{windowfun}{The windowing function to be used.}
  \item{method}{The "averaging" method to be used -- either
    \code{"mean"} or \code{"median"}.}
  \item{windowingPsdCorrection}{a \code{logical} flag indicating whether
  an overall correction for the windowing is supposed to be applied to
  the result -- this essentially specifies whether the result is
  supposed to reflect the PSD of the \emph{windowed} or of the
  "\emph{un-windowed"} time series.}
  \item{\dots}{other parameters passed to the windowing function.}
}
\details{
  The time series will be divided into overlapping sub-segments, each
  segment is windowed and its "empirical" spectrum is computed. The
  average of these spectra is the resulting PSD estimate. For
  robustness, the median may also be used instead of the mean.
}
\value{
  A list containing the following elements:
  \item{frequency}{the Fourier frequencies.}
  \item{power}{the estimated spectral power.}
  \item{kappa}{the number of (by definition) non-zero imaginary
    components of the Fourier series.}
  \item{two.sided}{a \code{logical} flag indicating one- or
    two-sidedness.}
  \item{segments}{the number of (overlapping) segments used.}
}
\references{Welch, P. D. (1967):
  \href{http://dx.doi.org/10.1109/TAU.1967.1161901}{The use of Fast
  Fourier Transform for the estimation of Power Spectra:
  A method based on time averaging over short, modified
  periodograms}. \emph{IEEE Transactions on Audio and Electroacoustics},
  AU-15(2):70--73.

  Press, W. H., Teukolsky, S. A., Vetterling, W. T., Flannery,
  B. P. (1992): \emph{Numerical recipes in C}. Cambridge University
  Press.
}
\author{Christian Roever, \email{bspec@web.de}}
\seealso{\code{\link{empiricalSpectrum}},
  \code{\link{tukeywindow}},
  \code{\link{spectrum}}
}
\examples{
# load example data:
data(sunspots)
# compute and plot the "plain" spectrum:
spec1 <- empiricalSpectrum(sunspots)
plot(spec1$frequency, spec1$power, log="y", type="l")

# plot Welch spectrum using segments of length 10 years:
spec2 <- welchPSD(sunspots, seglength=10)
lines(spec2$frequency, spec2$power, col="red")

# use 20-year segments and a flatter Tukey window:
spec3 <- welchPSD(sunspots, seglength=20, r=0.2)
lines(spec3$frequency, spec3$power, col="blue")
}
\keyword{ts}
