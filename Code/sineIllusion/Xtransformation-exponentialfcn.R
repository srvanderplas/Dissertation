library(ggplot2)

createSine <- function(n=200, len=1, f=f, fprime=fprime, getquadapprox=FALSE, f2prime=getquadapprox) {
  if(getquadapprox & !is.function(f2prime)) f2prime <- function(x) -1*f(x) # for backwards compatibility
  x <- seq(0, 2*pi, length=n+2)[(2:(n+1))]
  l <- rep(len, length=length(x))
  fx <- f(x)
  ystart <- fx - .5*l
  yend <- fx + .5*l
  ell <- yend-ystart
  # now correct for line illusion in vertical direction
  dy <- diff(range(fx))
  dx <- diff(range(x))
  # fprime works in framework of dx and dy, but we represent it in framework of dx and dy+len
  # needs to be fixed by factor a:  
  a <- dy/(dy + len) 
  # ellx is based on the "trig" correction
  ellx <- ell / cos(atan(abs(a*fprime(x))))
  # ellx2 is based on linear approximation of f  
  ellx2 <- ell * sqrt(1 + a*fprime(x)^2)
  
  # make this a data frame - ggplot2 doesn't do well with floating vectors
  dframe <- data.frame(x=x, xstart=x, xend=x, y=fx, ystart=ystart, yend=yend, ell=ell, ellx = ellx, ellx2=ellx2)
  
  # third adjustment is based on quadratic approximation of f.
  # this needs two parts: correction above and below f(x)  
  if(getquadapprox & is.function(f2prime)){
    secseg <- do.call("rbind", lapply(dframe$x, function(i) getSecantSegment(i, dframe, f, fprime, f2prime)))
    dframe$ellx3.u <- secseg$sec.ell1
    dframe$ellx3.l <- secseg$sec.ell2
  }
  dframe
}

correctx <- function(z, fprime, a=0, b=2*pi, shrink=0) {
  const <- integrate(function(x) abs(fprime(x)), a, b)$value
  trans <- sapply(z, function(i) integrate(function(x) abs(fprime(x)), a, i)$value*(b-a)/const+a)
  rowMeans(matrix(c(rep(z, times=shrink), trans), nrow=length(z), byrow=FALSE))
}

f <- exp
fprime <- exp
dframe <- data.frame(x=seq(0, 2, by=.1), y=f(seq(0, 2, by=.1)))
dframe$yend <- dframe$y-.5
dframe$ystart <- dframe$y+.5
dframe$xstart <- dframe$x
dframe$xend <- dframe$x
minor.axis.correction <- correctx(seq(0, 2, .1), fprime)

dframe$xtrans <- correctx(dframe$x, fprime=fprime)

dots <- data.frame(x = rep(minor.axis.correction, times=1), y=rep(c(0), each=length(minor.axis.correction)))

qplot(x=xtrans, xend=xtrans, y = ystart, yend=yend, colour=I("blue"), geom="segment", data=dframe) + theme_bw() + xlab("x") + ylab("y")+ geom_point(data=dots, aes(x=x, y=y), inherit.aes=FALSE) 



f <- function(x) (x-.5)*(x-2)^2
fprime <- function(x) (x-2)^2 + 2*(x)*(x-2) - x + 2
dframe <- data.frame(x=seq(0, 3, by=.1), y=f(seq(0, 3, by=.1)))
dframe$yend <- dframe$y-.5
dframe$ystart <- dframe$y+.5
dframe$xstart <- dframe$x
dframe$xend <- dframe$x
minor.axis.correction <- correctx(seq(0, 3, .1), fprime, a=0, b=3)
major.breaks <- correctx(seq(0, 3, .5), fprime, a=0, b=3)

dframe$xtrans <- correctx(dframe$x, fprime=fprime, a=0, b=3, shrink=0)
dframe$xint.f <- dframe$x*(f(dframe$x)==0)
dframe$xint.fprime <- dframe$x*(fprime(dframe$x)==0)

dots <- data.frame(x = rep(minor.axis.correction, times=1), y=rep(c(0), each=length(minor.axis.correction)))

qplot(x=x, xend=x, y=ystart, yend=yend, colour=I("blue"), geom="segment", data=dframe) + theme_bw() + xlab("x") + ylab("y") 
qplot(x=xtrans, xend=xtrans, y = ystart, yend=yend, colour=I("blue"), geom="segment", data=dframe) + theme_bw() + xlab("x") + ylab("y")+ geom_point(data=dots, aes(x=x, y=y), inherit.aes=FALSE) + scale_x_continuous(minor_breaks=minor.axis.correction, breaks=major.breaks, labels=seq(0, 3, .5))
