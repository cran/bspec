\name{snr}
\alias{snr}
\title{Compute the signal-to-noise ratio (SNR) of a signal}
\description{
  Compute the SNR for a given signal and noise power spectral density.
}
\usage{
snr(x, psd, two.sided = FALSE)
}
\arguments{
  \item{x}{the signal waveform, a time series (\code{\link{ts}}) object.}
  \item{psd}{the noise power spectral density. May be a vector of
    appropriate length (\code{length(x)/2+1}) or a function of frequency.}
  \item{two.sided}{a \code{logical} flag indicating whether the
      \code{psd} argument is to be interpreted as a one-sided or a
      two-sided spectrum.}
}
\details{For a signal \eqn{s(t)}{s(t)}, the complex-valued discrete
  Fourier transform \eqn{\tilde{s}(f)}{s(f)} is computed along the Fourier
  frequencies \eqn{f_j=\frac{j}{N \Delta_t} |
    j=0,\ldots,N/2+1}{f[j]=j/(N*delta_t) | j=0,...,N/2+1}, where
  \eqn{N}{N} is the sample size, and \eqn{\Delta_t}{delta_t} is the
  sampling interval.
  The SNR, as a measure of "signal strength" relative to the noise, then
  is given by
  \deqn{\varrho=\sqrt{\sum_{j=0}^{N/2+1}\frac{\bigl|\tilde{s(f_j)}\bigr|^2}{\frac{N}{4\Delta_t} S_1(f_j)}},}{rho=sqrt(sum(abs(s(f))^2 / ((N/(4*delta_t)) * S1(f))),}
  where \eqn{S_1(f)}{S1(f)} is the noise's one-sided power spectral
  density. For more on its interpretation, see e.g. Sec. II.C.4 in the
  reference below.
}
\value{
  The SNR \eqn{\varrho}{rho}.
}
\references{Roever, C. (2011): A Student-t based filter for robust
  signal detection. Submitted for publication;
  \href{http://arxiv.org/abs/1109.0442}{arXiv preprint 1109.0442}.}
\author{Christian Roever, \email{bspec@web.de}}
\seealso{\code{\link{matchedfilter}}, \code{\link{studenttfilter}}}
\examples{
# sample size and sampling resolution:
N       <- 1000
deltaT  <- 0.001

# For the coloured noise, use some AR(1) process;
# AR noise process parameters:
sigmaAR <- 1.0
phiAR   <- 0.9

# generate non-white noise
# (autoregressive AR(1) low-frequency noise):
noiseSample <- rnorm(10*N)
for (i in 2:length(noiseSample))
  noiseSample[i] <- phiAR*noiseSample[i-1] + noiseSample[i]
noiseSample <- ts(noiseSample, deltat=deltaT)

# estimate the noise spectrum:
PSDestimate <- welchPSD(noiseSample, seglength=1,
                        windowingPsdCorrection=FALSE)

# generate a (sine-Gaussian) signal:
t0    <- 0.6
phase <- 1.0
t <- ts((0:(N-1))*deltaT, deltat=deltaT, start=0)
signal <- exp(-(t-t0)^2/(2*0.01^2)) * sin(2*pi*150*(t-t0)+phase)
plot(signal)

# compute the signal's SNR:
snr(signal, psd=PSDestimate$power)
}
\keyword{ts}
