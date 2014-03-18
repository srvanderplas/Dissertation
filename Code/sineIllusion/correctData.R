# Functions to correct actual data in x and y
library(ggplot2)
correct.x <- function(x, y, data, model, f, fprime, f2prime=NULL, weight=NULL, ratio=1){ 
  if(length(weight)<1) weight=.36
  # since .361 is the midpoint of the left and right estimates of optimal values.

  a <- min(data[,x])
  b <- max(data[,x])
  
  const <- integrate(function(z) abs(fprime(z)), a, b, subdivisions=1000*length(unique(data[,x])), abs.tol=0.001, rel.tol=0.001)$value
  trans <- sapply(data[,x], function(i) integrate(function(z) abs(fprime(z)), a, i, subdivisions=1000*length(unique(data[,x])), abs.tol=0.001, rel.tol=0.001)$value*(b-a)/const + a)
  #   const <- sum(abs(fprime(data[,x])))
  #   trans <- cumsum(abs(fprime(data[,x])))*(b-a)/const+a                  
  
  data[,paste(x, "correctx", sep=".")] <- trans*weight + data[,x]*(1-weight)
  data[,paste(y, "fit", sep=".")] <- f(data[,x])
  
  data
}

correct.y <- function(x, y, data, model, f, fprime, f2prime=NULL, weight=NULL, ratio=1){
  if(length(weight)<1) weight = .40
  # since .407 is the midpoint of the left and right estimates of optimal values.
  
  
  df <- data.frame(x=data[,x], y=data[,y], 
                   pred=f(data[,x]), 
                   deriv=fprime(data[,x]), 
                   deriv2 = f2prime(data[,x]))
  df$resid <- df$y - df$pred
  
  dy <- diff(range(df$pred))
  dy1 <- diff(range(df$y))
  dx <- diff(range(df$x))
#   a <- dy/dy1  # aspect ratio : line "length" correction
  a <- ratio
  df$ell <-abs(df$resid)
  
  # linear correction
  cor <- (weight*sqrt(1 + a^2*fprime(df$x)^2) + (1-weight))
  data[,paste(y, "correcty", sep=".")] <- df$pred + df$resid*cor
#   

#   # quadratic correction
#   fp <- a*fprime(data[,x])
#   f2p <- a*f2prime(data[,x])
#   v <- 1 + fp^2
#   
#   a2 <- .5*fp^2*f2p
#   b2 <- v
#   c2 <- df$ell/2
#   lambdapinv <- 0.5*(sqrt(v^2-f2p*fp^2*df$ell) + v)    
#   lambdaminv <- 0.5*(sqrt(v^2+f2p*fp^2*df$ell) + v)
#   
#   lower <- abs(lambdapinv)/sqrt(v)
#   upper <- abs(lambdaminv)/sqrt(v)
#   data[,paste(y, "correcty", sep=".")] <- df$pred + (
#     (df$resid>=0)*upper*df$resid + (df$resid<0)*lower*df$resid
#     )
  
  data[,paste(y, "fit", sep=".")] <- f(data[,x])
  data

}

correct <- function(x, y, data, type, model, f, fprime, f2prime=NULL, weight=NULL, ratio=1){
  if(type=="x") correct.x(x, y, data, model, f, fprime, f2prime=NULL, weight, ratio)
  else correct.y(x, y, data, model, f, fprime, f2prime=NULL, weight, ratio)
}

# datasub <- read.csv("data/Ozone/Ozone-subset.csv")
# temp <- correct("Tmax", "Ozone", datasub, type="x")
# qplot(data=temp, x=Tmax.correctx, y=Ozone)
# temp2 <- correct("Tmax", "Ozone", datasub, type="y")
# qplot(data=temp2, x=Tmax, y=Ozone.correcty)
