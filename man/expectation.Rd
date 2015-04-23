\name{expectation}
\alias{expectation}
\alias{variance}
\alias{expectation.bspec}
\alias{expectation.bspecACF}
\alias{variance.bspec}
\alias{variance.bspecACF}
\title{Expectations and variances of distributions}
\description{
  Functions to compute (posterior) expectations or variances of the
  distributions specified as arguments.
}
\usage{
  expectation(x, ...)
  variance(x, ...)
  \method{expectation}{bspec}(x, two.sided=x$two.sided, ...)
  \method{variance}{bspec}(x, two.sided=x$two.sided, ...)
  \method{expectation}{bspecACF}(x, ...)
  \method{variance}{bspecACF}(x, ...)
}
\arguments{
  \item{x}{A \code{\link[=bspec.default]{bspec}} or
    \code{\link[=acf.bspec]{bspecACF}} object.}
  \item{two.sided}{A \code{logical} flag to indicate whether to compute
    expectation / variance of one- or two-sided spectrum, \emph{if}
    the argument \code{x} is a \code{bspec} object.}
  \item{...}{currently unused.}
}
\value{
  A \code{numeric} vector giving the expectations/variances
  corresponding to the frequencies or lags of the argument.
}
\references{Roever, C., Meyer, R., Christensen, N.
  \href{http://dx.doi.org/10.1088/0264-9381/28/1/015010}{Modelling
    coloured residual noise in gravitational-wave signal processing}.
  Classical and Quantum Gravity, 28(1):015010, 2011.
  See also \href{http://arxiv.org/abs/0804.3853}{arXiv preprint 0804.3853}.}
\author{Christian Roever, \email{christian.roever@med.uni-goettingen.de}}
\seealso{\code{\link[=bspec.default]{bspec}}, \code{\link[=acf.bspec]{acf.bspec}}}
\examples{
# note the changing expectation
# with increasing prior/posterior degrees-of-freedom:
expectation(bspec(lh))
expectation(bspec(lh, priordf=1, priorscale=0.6))
expectation(bspec(lh, priordf=2, priorscale=0.6))

# similar for variance:
variance(bspec(lh, priordf=2, priorscale=0.6))
variance(bspec(lh, priordf=3, priorscale=0.6))

# and again similar for autocovariances:
expectation(acf(bspec(lh)))
expectation(acf(bspec(lh, priordf=2, priorscale=0.6)))
variance(acf(bspec(lh)))
variance(acf(bspec(lh, priordf=4, priorscale=0.6)))
}
\keyword{ts}
