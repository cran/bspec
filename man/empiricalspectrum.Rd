\name{empiricalSpectrum}
\alias{empiricalSpectrum}
\title{Compute the "empirical" spectrum of a time series.}
\description{
  Computes the "empirical power" of a time series via a discrete Fourier
  transform.
}
\usage{
empiricalSpectrum(x, two.sided=FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{a time series (\code{\link[stats:ts]{ts}} object).}
  \item{two.sided}{a \code{logical} flag indicating whether the output
    is supposed to correspond to the \emph{one-sided} (default) or
    \emph{two-sided} spectrum.}
}
\details{
  Performs a Fourier transform, and then derives (based on the
  additional information on sampling rate etc. provided via the time
  series' attributes) the spectral power as a function of frequency.
  The result is simpler (in a way) than the \code{\link[stats:spectrum]{spectrum}()}
  function's output, see also the example below. What is returned is the
  real-valued frequency series
  \deqn{\kappa_j\frac{\Delta_t}{N}\bigl|\tilde{x}(f_j)\bigr|^2}{kappa[j]
    * (deltat/N) * abs(DFT(x)(f[j]))^2}
  where \eqn{j=0,...,N/2+1}{j=0,...,N/2+1},
  and \eqn{f_j=\frac{j}{N \Delta_t}}{f[j]=j / (N*deltat)} are the
  Fourier frequencies. \eqn{\Delta_t}{deltat} is the time series'
  sampling interval and \eqn{N}{N} is its
  length. \eqn{\kappa_j}{kappa[j]} is =1 for zero and Nyquist
  frequencies, and =2 otherwise, and denotes the number of (by
  definition) non-zero Fourier coefficients. In case
  \code{two.sided=TRUE}, the \eqn{\kappa_j}{kappa[j]} prefactor is
  omitted.

  For actual spectral estimation purposes, the use of a windowing
  function (see e.g. the \code{\link{tukeywindow}()} function) is highly
  recommended.
}
\value{
  A list containing the following elements:
  \item{frequency}{the Fourier frequencies.}
  \item{power}{the spectral power.}
  \item{kappa}{the number of (by definition) non-zero imaginary
    components of the Fourier series.}
  \item{two.sided}{a \code{logical} flag indicating one- or
    two-sidedness.}
}
\references{Roever, C., Meyer, R., Christensen, N. (2011):
  \href{http://dx.doi.org/10.1088/0264-9381/28/1/015010}{Modelling
    coloured residual noise in gravitational-wave signal processing}.
  \emph{Classical and Quantum Gravity}, 28(1):015010.
  See also \href{http://arxiv.org/abs/0804.3853}{arXiv preprint 0804.3853}.}
\author{Christian Roever, \email{bspec@web.de}}
\seealso{\code{\link[stats:fft]{fft}},
  \code{\link[stats:spectrum]{spectrum}},
  \code{\link{tukeywindow}},
  \code{\link{welchPSD}}
}
\examples{
# load example data:
data(lh)

# compute spectrum:
spec1 <- empiricalSpectrum(lh)
plot(spec1$frequency, spec1$power, log="y", type="b")

# plot "spectrum()" function's result in comparison:
spec2 <- spectrum(lh, plot=FALSE)
lines(spec2$freq, spec2$spec, col="red")

# make both spectra match:
spec3 <- empiricalSpectrum(lh, two.sided=TRUE)
spec4 <- spectrum(lh, plot=FALSE, taper=0, fast=FALSE, detrend=FALSE)
plot(spec3$frequency, spec3$power, log="y", type="b")
lines(spec4$freq, spec4$spec, col="green")
}
\keyword{ts}
