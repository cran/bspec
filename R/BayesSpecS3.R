#
# bspec: Bayesian spectral inference.
#
# S3 Version.
#
# http://cran.r-project.org/package=bspec
#
#
# see also:
#
#   Roever, C., Meyer, R. and Christensen, N. (2010):
#   Modelling coloured noise.
#   Arxiv preprint 0804.3853
#   URL http://arxiv.org/abs/0804.3853
#

print.bspec <- function(x, ...)
{
  if (x$two.sided) cat(" 'bspec' posterior spectrum (two-sided).\n")
  else cat(" 'bspec' posterior spectrum (one-sided).\n")
  if (is.element("temperature", names(x))){
    cat(paste(" (tempered: T = ",x$temperature,")\n",sep=""))
  }
  cat(paste(" frequency range     : ",
            as.character(signif(x$freq[1],4)),
            "--",
            as.character(signif(x$freq[length(x$freq)],4)),
            "\n",sep=""))
  cat(paste(" number of parameters: ",length(x$scale),"\n",sep=""))
  expect <- expectation(x)
  cat(" finite expectations : ")
  if (all(is.finite(expect))) cat("all\n")
  else if (any(is.finite(expect))) cat("some\n")
  else cat("none\n")
  vari <- variance(x)
  cat(" finite variances    : ")
  if (all(is.finite(vari))) cat("all\n")
  else if (any(is.finite(vari))) cat("some\n")
  else cat("none\n")
  cat(" call: "); print(x$call)
  invisible(x)
}



print.bspecACF <- function(x, ...)
{
  if (x$type=="correlation")
    cat(paste("Autocorrelations of bspec object '",x$bspec,"', by lag\n", sep=""))
  else
    cat(paste("Autocovariances of bspec object '",x$bspec,"', by lag\n", sep=""))
  printvec <- x$acf
  names(printvec) <- x$lag
  print(printvec)
  invisible(x)
}



bspec.default <- function(x, priorscale=1, priordf=0,
                          intercept=TRUE, two.sided=FALSE, ...)
# x                   : time series object (uni- or multivariate)
# priorscale, priordf : a vector or a function of frequency
# intercept           : flag indicating whether to include zero frequency
# two.sided           : flag indicating whether to refer to one- or two-sided spectrum.
#                       only effect within _this_ function is interpretation of prior scale.
#                       Note that the 'two.sided' flag is "inherited" by other functions
#                       via their default arguments. See e.g. 'expectation.bspec()'.
{
  # some initial checks:
  if (is.vector(x) | (is.ts(x) & !is.mts(x)))
    N <- length(x)
  else if (is.mts(x) | is.matrix(x) | is.data.frame(x)) 
    N <- nrow(x)
  else warning("incompatible argument 'x'")
  FTlength <- (N %/% 2) + 1   # size of (nonredundant) FT output
  Neven <- ((N %% 2) == 0)    # indicator for even N
  stopifnot(is.function(priorscale) ||
            ((length(priorscale)==1) |
            (intercept & (length(priorscale)==FTlength)) |
            ((!intercept) & (length(priorscale)==(FTlength-1)))))
  stopifnot(is.function(priordf) ||
            ((length(priordf)==1) |
            (intercept & (length(priordf)==FTlength)) |
            ((!intercept) & (length(priordf)==(FTlength-1)))))
  stopifnot(is.function(priorscale) || (all(is.finite(priorscale)) & all(priorscale>0)),
            is.function(priordf) || (all(is.finite(priordf)) & all(priordf>=0)))
  if (!is.ts(x)) {
    x <- as.ts(x)
    warning("argument 'x' is not a time-series object, default conversion 'as.ts(x)' applied.")
  }
  deltat <- 1 / tsp(x)[3]
  deltaf <- 1 / (N*deltat)
  t0 <- tsp(x)[1]           # time stamp corresponding to 1st observation
  kappa <- c(0, rep(1,FTlength-2), ifelse(Neven,0,1))
  # 1-D case:
  if (!is.mts(x)) {
    # Fourier transform:
    y <- fft(as.vector(x))
    # (`fft()' yields unnormalised FT)
    nonredundant <- 1:FTlength
    # vector of (N/2 + 1) cosine coefficients:
    a <- (1+kappa) * sqrt(deltat/N) * Re(y[nonredundant])
    # vector of (N/2 + 1) sine coefficients:
    b <- -(1+kappa) * sqrt(deltat/N) * Im(y[nonredundant])
    datassq <- a^2 + b^2
    datadf <- c(1, rep(2,FTlength-2), ifelse(Neven,1,2))
  }
  # multidimensional case:
  else {
    # Fourier transform:
    #y <- apply(x, 2, fft)
    y <- mvfft(x)
    # (`fft()' yields unnormalised FT)
    nonredundant <- 1:FTlength
    # vector of (N/2 + 1) cosine coefficients:
    a <- sqrt(deltat/N) * Re(y[nonredundant,])
    for (i in 1:ncol(a)) a[,i] <- (1+kappa)*a[,i]
    # vector of (N/2 + 1) sine coefficients:
    b <- -sqrt(deltat/N) * Im(y[nonredundant,])
    for (i in 1:ncol(b)) b[,i] <- (1+kappa)*b[,i]
    datassq <- apply(a^2 + b^2, 1, sum)
    datadf <- c(1, rep(2,FTlength-2), ifelse(Neven,1,2)) * ncol(x)
  }
  # vector of corresponding frequencies:
  freq <- (0:(FTlength-1)) * deltaf
  
  # set up a-priori scale:
  if (is.function(priorscale))
    priorscalevec <- priorscale(freq)
  else if (length(priorscale)==1)
    priorscalevec <- rep(priorscale, FTlength)
  else if (!intercept)
    priorscalevec <- c(0, priorscale)
  else priorscalevec <- priorscale
  if (two.sided){
    #  /!\  different interpretation of prior scale  /!\
    priorscalevec <- priorscalevec * (1+kappa)
  }
  
  # set up a-priori degrees-of-freedom:
  if (is.function(priordf))
    priordfvec <- priordf(freq)
  if (length(priordf)==1)
    priordfvec <- rep(priordf, FTlength)
  else if (!intercept)
    priordfvec <- c(0, priordf)
  else priordfvec <- priordf
  
  #arg <- seq(from=tsp(x)[1], to=tsp(x)[2], le=500)
  #trigo <- arg*0
  #for (i in 1:length(a))
  # trigo<-trigo+(a[i]*cos(2*pi*freq[i]*(arg-t0)) + b[i]*sin(2*pi*freq[i]*(arg-t0)))
  #trigo <- trigo / sqrt(N*deltat)
  #plot(x,type="b")
  #lines(arg, trigo, col="red")

  if (! intercept) {
    freq <- freq[-1]
    priorscalevec <- priorscalevec[-1]
    priordfvec <- priordfvec[-1]
    datassq <- datassq[-1]
    datadf <- datadf[-1]
    kappa <- kappa[-1]
  }  
  
  # determine posterior distribution's parameters:
  # (posterior distn. of 1-sided spectrum S_1(f_j) ==  sigmasquared_j!)
  #  -->  S_2(f_j)  =  sigmasquared_j / (1+kappa(j))  =  S_1(f_j) / (1+kappa(j))
  scale <- (priordfvec*priorscalevec + datassq) / (priordfvec + datadf)
  df <- datadf + priordfvec
  
  result <- list(freq = freq,
                 scale = scale,
                 df = df,
                 priorscale = priorscalevec,
                 priordf = priordfvec,
                 datassq = datassq,
                 datadf = datadf,
                 N = N, deltat = deltat, deltaf = deltaf,
                 start = t0,
                 call = match.call(expand.dots=FALSE),
                 two.sided = two.sided)
  class(result) <- "bspec"
  return(result)
}



expectation.bspec <- function(x, two.sided=x$two.sided, ...)
{
  expect <- rep(0, length(x$scale))
  finite <- x$df > 2
  expect[finite] <- x$scale[finite] * (x$df[finite]/(x$df[finite]-2))
  if (two.sided){
    Neven <- ((x$N %% 2) == 0)
    index <- 2:ifelse(Neven, x$N %/% 2, (x$N %/% 2) + 1)
    if (x$freq[1] != 0) index <- index-1
    expect[index] <- expect[index] / 2
  }
  expect[!finite] <- Inf
  return(expect)
}



variance.bspec <- function(x, two.sided=x$two.sided, ...)
{
  vari <- rep(0, length(x$scale))
  finite <- (x$df > 4)
  vari[finite] <- (2*x$df[finite]^2*x$scale[finite]^2) / ((x$df[finite]-2)^2 * (x$df[finite]-4))
  if (two.sided){
    Neven <- ((x$N %% 2) == 0)
    index <- 2:ifelse(Neven, x$N %/% 2, (x$N %/% 2) + 1)
    if (x$freq[1] != 0) index <- index-1
    vari[index] <- vari[index] / 4
  }
  vari[!finite] <- Inf
  return(vari)
}



sample.bspec <- function(x, size=1, two.sided=x$two.sided, ...)
{
  rinvchisq <- function(n=1, df=1, scale=1)
  # sample n times from an Inverse-Chi-Squared distribution
  # with scale "scale" and degrees-of-freedom "df".
  {
    return(scale * (df/rchisq(n=n,df=df)))
  }
  sigmasq <- matrix(rinvchisq(n = size*length(x$freq),
                              df = x$df, scale = x$scale),
                    ncol=size)
  if (two.sided){
    Neven <- ((x$N %% 2) == 0)
    index <- 2:ifelse(Neven, x$N %/% 2, x$N %/% 2 + 1)
    if (x$freq[1] != 0) index <- index-1
    sigmasq[index,] <- sigmasq[index,] / 2
  }
  if (size==1) sigmasq <- as.vector(sigmasq)
  return(sigmasq)
}



quantile.bspec <- function(x, probs = c(0.025,0.5,0.975),
                           two.sided=x$two.sided, ...)
{
  qinvchisq <- function(p=0.5, df=1, scale=1)
  {
    return((df*scale)/qchisq(1-p,df=df))
  }
  sigmasq <- matrix(qinvchisq(p = rep(probs,each=length(x$freq)),
                              df = x$df, scale = x$scale),
                    ncol=length(probs))
  if (two.sided){
    Neven <- ((x$N %% 2) == 0)
    index <- 2:ifelse(Neven, x$N %/% 2, x$N %/% 2 + 1)
    if (x$freq[1] != 0) index <- index-1
    sigmasq[index,] <- sigmasq[index,] / 2
  }
  if (length(probs)==1) sigmasq <- as.vector(sigmasq)
  else colnames(sigmasq) <- as.character(signif(probs))
  return(sigmasq)
}



plot.bspec <- function(x, two.sided=x$two.sided, ...)
{
  oldpar <- par(no.readonly=TRUE)
  quant <- quantile(x, probs=c(0.025,0.5,0.975),
                    two.sided=two.sided)
  plot(range(x$freq), range(quant),type="n", axes=FALSE,
       log="y", xlab="", ylab="", ...)
  matlines(rbind(x$freq,x$freq),
           t(quant[,c(1,3)]), col="darkgrey", lty="solid")
  points(x$freq, quant[,2], pch=5)
  expect <- expectation(x, two.sided=two.sided)
  if (any(is.finite(expect)))
    points(x$freq, expect, pch=20)
  axis(1)
  mtext(expression("frequency "*italic(f)),
        side=1, line=par("mgp")[1], cex=par("cex.lab"))
  axis(2)
  if (two.sided)
    mtext(expression("(two-sided) posterior spectrum "*italic(S(f))),
          side=2, line=par("mgp")[1], cex=par("cex.lab"))
  else
    mtext(expression("(one-sided) posterior spectrum "*italic(S(f))),
          side=2, line=par("mgp")[1], cex=par("cex.lab"))
  sqrtrange <- par("usr")[3:4]/2
  ticks <- axTicks(side=4, usr=sqrtrange, log=par("ylog"))
  axlabel <- NULL
  for (i in 1:length(ticks))
    axlabel <- c(axlabel, eval(bquote(expression((.(ticks[i]))^2))))     
  axis(4, at=ticks^2, label=axlabel)
  box()
  result <- list("freq"=x$freq,
                 "spectrum"=cbind("0.5 %"=quant[,1],
                                  "median"=quant[,2],
                                  "99.5 %"=quant[,3],
                                  "mean"=expect))
  invisible(result)
}



ppsample.bspec <- function(x, start=x$start, ...)
# sample from posterior predictive distribution
{
  kappa <- c(0, rep(1,(x$N %/% 2) -1), ifelse(((x$N %% 2) == 0),0,1))
  # sample spectrum:
  spectrumsample <- sample(x, size=1, two.sided=TRUE)
  # sample data (in Fourier domain):
  Acoefsample <- rnorm(n=length(spectrumsample), mean=0, sd=sqrt(spectrumsample/(1+kappa)))
  Bcoefsample <- rnorm(n=length(spectrumsample), mean=0, sd=sqrt(spectrumsample/(1+kappa)))
  real <- sqrt(x$N/x$deltat) * Acoefsample
  imag <- -sqrt(x$N/x$deltat) * Bcoefsample
  real <- c(real, rev(real[kappa>0]))
  imag <- c(imag, -rev(imag[kappa>0]))
  noiseFT <- complex(x$N, real=real, imaginary=imag)
  # transform to time domain:
  noiseTS <- Re(fft(noiseFT, inverse=TRUE)) / x$N
  noiseTS <- ts(noiseTS, start=0, deltat=x$deltat)
  return(noiseTS)
}



acf.bspec <- function(x, spec=NULL,
                      type=c("covariance", "correlation"),
                      two.sided=x$two.sided, ...)
{
  if (is.null(spec)){
    spec <- expectation(x, two.sided=FALSE)
    vars <- variance(x, two.sided=FALSE)
    if (all(is.finite(vars)))
      estimate.errors <- TRUE
    else {
      stderr <- rep(Inf, x$N%/%2 + 1)
      estimate.errors <- FALSE
    }
  }
  else {
    stderr <- rep(0, x$N%/%2 + 1)
    estimate.errors <- FALSE
    if (two.sided){
      if ((x$N%%2)==0)
        multi <- c(1, rep(2,(x$N %/% 2)-1), 1)
      else
        multi <- c(1, rep(2,(x$N %/% 2)-1))
      spec <- multi * spec
    }        
  }
  type <- match.arg(type)
  lags <- (0:(x$N%/%2)) * x$deltat
  kappa <- c(0, rep(1,(x$N %/% 2) -1), ifelse(((x$N %% 2) == 0),0,1))
  if (x$freq[1] != 0) kappa <- kappa[-1]
  if (all(is.finite(spec))){
    autocov1 <- function(deltat)
    {
      expect <- sum(spec*(1+kappa)/2*cos(2*pi*x$freq * deltat)) / (x$N*x$deltat)
      return(c(expect, NA))
    }
    autocov2 <- function(deltat)
    {
      coef <- (1+kappa) * 0.5 * cos(2*pi*x$freq * deltat) / (x$N*x$deltat)
      expect <- sum(spec*coef) 
      stdev <- sqrt(sum(vars*coef^2))
      return(c(expect,stdev))
    }
    if (estimate.errors){
      dummy <- apply(matrix(lags,ncol=1),1,autocov2)
      acf <- dummy[1,]
      stderr <- dummy[2,]
    }
    else {
      dummy <- apply(matrix(lags,ncol=1),1,autocov1)
      acf <- dummy[1,]
    }
  }
  else acf <- rep(Inf, length(lags))
  result <- list(lag = lags,
                 acf = acf,
                 stderr = stderr,
                 type = type,
                 N = x$N,
                 bspec = deparse(substitute(x)))
  if (type=="correlation")
    result$acf <- result$acf/result$acf[1]
  class(result) <- "bspecACF"
  return(result)
}



expectation.bspecACF <- function(x, ...)
{
  if (x$type == "covariance")
    result <- x$acf
  else {
    result <- rep(NA, length(x$acf))
    warning("'expectation()' only defined for autocovariances, not for autocorrelations")
  }
  return(result)
}



variance.bspecACF <- function(x, ...)
{
  if (x$type == "covariance")
    result <- x$stderr^2
  else {
    result <- rep(NA, length(x$acf))
    warning("'variance()' only defined for autocovariances, not for autocorrelations")
  }
  return(result)
}



plot.bspecACF <- function(x, ci=0.95,
                          type = "h", xlab = NULL, ylab = NULL,
                          main = NULL, ci.col = "blue", ...)
{
  if (is.null(xlab)) xlab <- "lag"
  if (is.null(ylab)) ylab <- ifelse(x$type=="correlation", "autocorrelation", "autocovariance")
  plot(x$lag, x$acf, type=type,
       xlab=xlab, ylab=ylab, ...)
  abline(h=0)
  if ((x$type=="covariance") && all(is.finite(x$stderr)) && (any(x$stderr>0))) {
    z <- qnorm(1-(1-ci)/2)
    cilines <- cbind(x$acf-z*x$stderr, x$acf+z*x$stderr)
    matlines(x$lag, cilines, col=ci.col, lty="dashed")
    
    z <- 1/sqrt(1-ci)
    cilines <- cbind(x$acf-z*x$stderr, x$acf+z*x$stderr)
    matlines(x$lag, cilines, col="red", lty="dashed")
  }
  invisible(x)
}



temper.bspec <- function(x, temperature=2.0, likelihood.only=TRUE, ...)
{
  tempered <- list(freq = x$freq,
                   scale = NA,
                   df = NA,
                   priorscale = x$priorscale,
                   priordf = x$priordf,
                   datassq = x$datassq,
                   datadf = x$datadf,
                   N = x$N,
                   deltat = x$deltat,
                   deltaf = x$deltaf,
                   start = x$start,
                   call = x$call,
                   two.sided = x$two.sided)
  if (temperature != 1) {
    tempered <- c(tempered, "temperature"=temperature)
    if (likelihood.only){
      tempered$scale <- (x$priordf*x$priorscale + (x$datadf/temperature)*x$datassq) / (x$priordf + x$datadf/temperature)
      tempered$df    <- x$priordf + (x$datadf)/temperature
    }
    else {
      stopifnot(temperature < (min(x$df)+2)/2)
      # ...for T < (nu+2)/2...
      tempered$scale <- (x$df*x$scale) / (2 + x$df - 2*temperature)
      tempered$df    <- (x$df+2)/temperature - 2
    }
  }
  else{
    tempered$scale <- (x$priordf*x$priorscale + x$datassq) / (x$priordf + x$datadf)
    tempered$df    <- x$datadf + x$priordf
  }
  class(tempered) <- "bspec"
  return(tempered)
}



temperature.bspec <- function(x, ...)
{
  if (is.element("temperature", names(x)))
    result <- x$temperature
  else
    result <- 1.0
  return(result)
}



dprior.bspec <- function(x, theta,
                         two.sided=x$two.sided, log=FALSE, ...)
# (log-) prior density
# `theta' is the 2-sided PSD (S_2(f_j))
{
  stopifnot(length(theta) == length(x$freq))
  if (any(theta<=0))
    prior <- -Inf
  else{
    if (two.sided) {
      Neven <- ((x$N %% 2) == 0)
      index <- 2:ifelse(Neven, x$N %/% 2, (x$N %/% 2) + 1)
      if (x$freq[1] != 0) index <- index-1
      theta[index] <- theta[index] * 2
    } 
    halfdf <- x$priordf/2
    prior <- sum(halfdf*log(halfdf) - lgamma(halfdf) + halfdf*log(x$priorscale) - (halfdf+1)*log(theta) - (x$priordf*x$priorscale)/(2*theta))
  }
  if (!log)
    prior <- exp(prior)
  return(prior)
}



likelihood.bspec <- function(x, theta,
                             two.sided=x$two.sided, log=FALSE, ...)
# (log-) likelihood
# `theta' is the 2-sided PSD (S_2(f_j))
{
  stopifnot(length(theta) == length(x$freq))
  if (any(theta<=0))
    likeli <- -Inf
  else {
    if (two.sided) {
      Neven <- ((x$N %% 2) == 0)
      index <- 2:ifelse(Neven, x$N %/% 2, (x$N %/% 2) + 1)
      if (x$freq[1] != 0) index <- index-1
      theta[index] <- theta[index] * 2
    } 
    likeli <- -0.5*sum(x$datadf)*log(2*pi) + sum(-(x$datadf/2)*log(theta) - x$datassq/(2*theta))
  }
  if (! log)
    likeli <- exp(likeli)
  return(likeli)
}



dposterior.bspec <- function(x, theta,
                             two.sided=x$two.sided, log=FALSE, ...)
# (log-) posterior density
# `theta' is the 2-sided PSD (S_2(f_j))
{
  stopifnot(length(theta) == length(x$freq))
  if (any(theta<=0))
    posterior <- -Inf
  else {
    if (two.sided) {
      Neven <- ((x$N %% 2) == 0)
      index <- 2:ifelse(Neven, x$N %/% 2, (x$N %/% 2) + 1)
      if (x$freq[1] != 0) index <- index-1
      theta[index] <- theta[index] * 2
    } 
    halfdf <- x$df/2
    posterior <- sum(halfdf*log(halfdf) - lgamma(halfdf) + halfdf*log(x$scale) - (halfdf+1)*log(theta) - (x$df*x$scale)/(2*theta))
  }
  if (!log)
    posterior <- exp(posterior)
  return(posterior)
}



is.bspec <- function(object)
{
  is.element("bspec", class(object))
}



is.bspecACF <- function(object)
{
  is.element("bspecACF", class(object))
}



one.sided.bspec <- function(x, ...)
{
  x$two.sided <- FALSE
  return(x)
}



two.sided.bspec <- function(x, ...)
{
  x$two.sided <- TRUE
  return(x)
}


##################################

bspec <- function(x, ...)
{
  UseMethod("bspec")
}


expectation <- function(x, ...)
{
  UseMethod("expectation")
}


variance <- function(x, ...)
{
  UseMethod("variance")
}


one.sided <- function(x, ...)
{
  UseMethod("one.sided")
}


two.sided <- function(x, ...)
{
  UseMethod("two.sided")
}


sample <- function(x, ...)
{
  UseMethod("sample")
}


sample.default <- base::sample
formals(sample.default) <- c(formals(sample.default), alist(... = ))


ppsample <- function(x, ...)
{
  UseMethod("ppsample")
}


acf <- function(x, ...)
{
  UseMethod("acf")
}

acf.default <- stats::acf


temper <- function(x, ...)
{
  UseMethod("temper")
}


temperature <- function(x, ...)
{
  UseMethod("temperature")
}


dprior <- function(x, ...)
{
  UseMethod("dprior")
}


dposterior <- function(x, ...)
{
  UseMethod("dposterior")
}


likelihood <- function(x, ...)
{
  UseMethod("likelihood")
}
