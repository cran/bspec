\name{matchedfilter}
\alias{matchedfilter}
\alias{studenttfilter}
\title{Filter a noisy time series for a signal of given shape}
\description{
  Computes the maximized likelihood ratio (as a test- or
  detection-statistic) of "\emph{signal}" vs. "\emph{noise only}"
  hypotheses. The signal is modelled as a linear combination of
  orthogonal basis vectors of unknown amplitude and arrival time. The
  noise is modelled either as Gaussian or Student-t-distributed, and
  coloured.
}
\usage{
matchedfilter(data, signal, noisePSD, timerange = NA,
    reconstruct = TRUE, two.sided = FALSE)

studenttfilter(data, signal, noisePSD, df = 10, timerange = NA,
    deltamax = 1e-06, itermax = 100,
    reconstruct = TRUE, two.sided = FALSE) 
}
\arguments{
  \item{data}{the data to be filtered, a time series (\code{\link[stats:ts]{ts}}) object.}
  \item{signal}{the signal waveform to be filtered for. May be a vector,
  a matrix, a time series (\code{\link[stats:ts]{ts}}) or a multivariate time
  series (\code{\link[stats:ts]{mts}}) object.}
  \item{noisePSD}{the noise power spectral density. May be a vector of
    appropriate length (\code{length(data)/2+1}) or a function of frequency.}
  \item{df}{the number of degrees-of-freedom (\eqn{\nu_j}{nu[j]}) for
    each frequency bin. May be a vector of appropriate length
    (\code{length(data)/2+1}) or a function of frequency.}
  \item{timerange}{the range of times (with respect to the \code{data}
    argument's time scale) to maximize the likelihood ratio over.}
  \item{deltamax}{the minimal difference in logarithmic likelihoods to
    be aimed for in the EM-iterations (termination criterion).}
  \item{itermax}{the maximum number of EM-iterations to be performed.}
  \item{reconstruct}{a \code{logical} flag indicating whether the
    output is supposed to include the best-fitting signal waveform.}
  \item{two.sided}{a \code{logical} flag indicating whether the
      \code{noisePSD} argument is to be interpreted as a one-sided or a
      two-sided spectrum.}
}
\details{The time series data \eqn{d(t)}{d(t)} is modelled as a
  superposition of signal \eqn{s(\beta,t_0,t)}{s(beta,t0,t)} and noise
  \eqn{n(t)}{n(t)}:
  \deqn{d(t)=s(\beta,t-t_0)+n(t),}{d(t)=s(beta,t-t0)+n(t),}
  where the signal is a linear combination of orthogonal (!) basis
  vectors \eqn{s_i(t)}{s[i](t)}, and whose time-of-arrival is given by
  the parameter \eqn{t_0}{t0}:
  \deqn{s(\beta,t-t_0)=\sum_{i=1}^k \beta_i s_i(t-t_0).}{s(beta,t-t0) =
    beta[1] s[1](t-t_0) + ... + beta[k] s[k](t-t_0).}
  The noise is modelled as either Gaussian (\code{matchedfilter()}) or
  Student-t distributed (\code{studenttfilter()}) with given power
  spectral density and, for the latter model only, degrees-of-freedom
  parameters.

  The filtering functions perform a likelihood maximization over the
  time-of-arrival \eqn{t_0}{t0} and coefficients \eqn{\beta}{beta}. In
  the Gaussian model, the conditional likelihood, conditional on
  \eqn{t_0}{t0}, can be maximized analytically, while the maximization
  over \eqn{t_0}{t0} is done numerically via a brute-force search. In
  the Student-t model, likelihood maximization is implemented using an
  EM-algorithm. The maximization over \eqn{t_0}{t0} is restricted to the
  range specified via the \code{timerange} argument.

  What is returned is the maximized (logarithmic)
  likelihood ratio of "signal" versus "noise-only" hypotheses (the
  result's \code{$maxLLR} component), and the corresponding
  ML-estimates \eqn{\hat{t}_0}{\code{tHat}} and
  \eqn{\hat{\beta}}{\code{betaHat}}, as well as the ML-fitted signal
  ("\code{$reconstruction}").
}
\value{
  A list containing the following elements:
  \item{maxLLR}{the maximized likelihood ratio of \emph{signal}
    vs. \emph{noise only} hypotheses.}
  \item{timerange}{the range of times to maximize the likelihood ratio
    over (see the \code{timerange} input argument).}
  \item{betaHat}{the ML-estimated vector of coefficients.}
  \item{tHat}{the ML-estimated signal arrival time.}
  \item{reconstruction}{the ML-fitted signal (a time series
    (\code{\link[stats:ts]{ts}}) object).}
  \item{call}{an object of class \code{call} giving the
    function call that generated the result.}
  elements only for the \code{matchedfilter()} function:
  \item{maxLLRseries}{the time series of (conditionally) maximized
    likelihood ratio for each given time point (the \emph{profile
    likelihood}).}
  elements only for the \code{studenttfilter()} function:
  \item{EMprogress}{a \code{matrix} indicating the progress of the EM-fitting.}
}
\references{Roever, C. (2011):
  \href{http://dx.doi.org/10.1103/PhysRevD.84.122004}{A Student-t based
    filter for robust signal detection}.
  \emph{Physical Review D}, 84(12):122004.
  See also \href{http://arxiv.org/abs/1109.0442}{arXiv preprint 1109.0442}.

  Roever, C., Meyer, R., Christensen, N. (2011):
  \href{http://dx.doi.org/10.1088/0264-9381/28/1/015010}{Modelling
    coloured residual noise in gravitational-wave signal processing}.
  \emph{Classical and Quantum Gravity}, 28(1):015010.
  See also \href{http://arxiv.org/abs/0804.3853}{arXiv preprint
    0804.3853}.}
\author{Christian Roever, \email{bspec@web.de}}
\examples{
# sample size and sampling resolution:
deltaT  <- 0.001
N       <- 1000

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

# show noise and noise PSD:
par(mfrow=c(2,1))
  plot(noiseSample, main="noise sample")
  plot(PSDestimate$freq, PSDestimate$pow, log="y", type="l",
       main="noise PSD", xlab="frequency", ylab="power")
par(mfrow=c(1,1))

# generate actual data:
noise <- rnorm(N)
for (i in 2:length(noise))
  noise[i] <- phiAR*noise[i-1] + noise[i]
noise <- ts(noise, start=0, deltat=deltaT)

# the "sine-Gaussian" signal to be injected into the noise:
t0    <- 0.6
phase <- 1.0
signal <- exp(-(time(noise)-t0)^2/(2*0.01^2)) * sin(2*pi*150*(time(noise)-t0)+phase)
plot(signal)

t <- seq(-0.1, 0.1, by=deltaT)
# the signal's orthogonal (sine- and cosine-) basis waveforms:
signalSine   <- exp(-t^2/(2*0.01^2)) * sin(2*pi*150*t)
signalCosine <- exp(-t^2/(2*0.01^2)) * sin(2*pi*150*t+pi/2)
signalBasis <- ts(cbind("sine"=signalSine, "cosine"=signalCosine),
                  start=-0.1, deltat=deltaT)
plot(signalBasis[,1], col="red", main="the signal basis")
lines(signalBasis[,2], col="green")
# (the sine- and cosine- components allow to span
#  signals of arbitrary phase)
# Note that "signalBasis" may be shorter than "data",
# but must be of the same time resolution.


# compute the signal's signal-to-noise ration (SNR):
signalSnr <- snr(signal, psd=PSDestimate$pow)

# scale signal to SNR = 6:
rho <- 6
data <- noise + signal * (rho/signalSnr)
data <- data * tukeywindow(length(data))
# Note that the data has (and should have!)
# the same resolution, size, and windowing applied
# as in the PSD estimation step.

# compute filters:
f1 <- matchedfilter(data, signalBasis, PSDestimate$power)
f2 <- studenttfilter(data, signalBasis, PSDestimate$power)

# illustrate the results:
par(mfrow=c(3,1))
  plot(data, ylab="", main="data")
  lines(signal* (rho/signalSnr), col="green")
  legend(0,max(data),c("noise + signal","signal only"),
         lty="solid", col=c("black","green"), bg="white")

  plot(signal * (rho/signalSnr), xlim=c(0.55, 0.65), ylab="",
       main="original & recovered signals")
  lines(f1$reconstruction, col="red")
  lines(f2$reconstruction, col="blue")
  abline(v=c(f1$tHat,f2$tHat), col=c("red", "blue"), lty="dashed")
  legend(0.55, max(signal*(rho/signalSnr)),
         c("injected signal","best-fitting signal (Gaussian model)",
           "best-fitting signal (Student-t model)"),
         lty="solid", col=c("black","red","blue"), bg="white")

  plot(f1$maxLLRseries, type="n", ylim=c(0, f1$maxLLR),
       main="profile likelihood (Gaussian model)",
       ylab="maximized (log-) likelihood ratio")
  lines(f1$maxLLRseries, col="grey")
  lines(window(f1$maxLLRseries, start=f1$timerange[1], end=f1$timerange[2]))
  abline(v=f1$timerange, lty="dotted")
  lines(c(f1$tHat,f1$tHat,-1), c(0,f1$maxLLR,f1$maxLLR), col="red", lty="dashed")
par(mfrow=c(1,1))
}
\keyword{ts}
\keyword{robust}
