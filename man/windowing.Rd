\name{tukeywindow}
\alias{tukeywindow}
\alias{squarewindow}
\alias{hannwindow}
\alias{welchwindow}
\alias{trianglewindow}
\alias{hammingwindow}
\alias{cosinewindow}
\alias{kaiserwindow}
\title{Compute windowing functions for spectral time series analysis.}
\description{
Several windowing functions for spectral or Fourier analysis of time
series data are provided.
}
\usage{
tukeywindow(N, r = 0.1)
squarewindow(N)
hannwindow(N)
welchwindow(N)
trianglewindow(N)
hammingwindow(N, alpha=0.543478261)
cosinewindow(N, alpha=1)
kaiserwindow(N, alpha=3)
}
\arguments{
  \item{N}{the length of the time series to be windowed}
  \item{r}{the Tukey window's parameter, denoting the its "non-flat"
    fraction.}
  \item{alpha}{additional parameter for Hamming-, cosine-, and Kaiser-windows.}
}
\details{
  \emph{Windowing} of time series data, i.e., multiplication with a tapering
  function, is often useful in spectral or Fourier analysis in order to
  reduce "leakage" effects due to the discrete and finite
  sampling. These functions provide windowing coefficients for a given
  sample size \code{N}.
}
\value{
  A vector (of length \code{N}) of windowing coefficients.
}
\references{Harris, F. J. (1978):
  \href{http://dx.doi.org/10.1109/PROC.1978.10837}{On the use of windows
  for harmonic analysis with the discrete Fourier transform}.
  \emph{Proceedings of the IEEE}, 66(1):51--83.

  Press, W. H., Teukolsky, S. A., Vetterling, W. T., Flannery,
  B. P. (1992): \emph{Numerical recipes in C}. Cambridge University
  Press.
}
\author{Christian Roever, \email{bspec@web.de}}
\seealso{\code{\link{welchPSD}}, \code{\link{empiricalSpectrum}}}
\examples{
# illustrate the different windows' shapes:
N <- 100
matplot(1:N,
        cbind(cosinewindow(N),
              hammingwindow(N),
              hannwindow(N),
              kaiserwindow(N),
              squarewindow(N),
              trianglewindow(N),
              tukeywindow(N,r=0.5),
              welchwindow(N)),
        type="l", lty="solid", col=1:8)
legend(N, 0.99, legend=c("cosine","hamming","hann","kaiser",
                         "square","triangle","tukey","welch"),
       col=1:8, lty="solid", xjust=1, yjust=1, bg="white")

# show their effect on PSD estimation:
data(sunspots)



spec1 <- welchPSD(sunspots, seglength=10, windowfun=squarewindow)
plot(spec1$frequency, spec1$power, log="y", type="l")

spec2 <- welchPSD(sunspots, seglength=10, windowfun=tukeywindow, r=0.25)
lines(spec2$frequency, spec2$power, log="y", type="l", col="red")

spec3 <- welchPSD(sunspots, seglength=10, windowfun=trianglewindow)
lines(spec3$frequency, spec3$power, log="y", type="l", col="green")
}
\keyword{ts}
