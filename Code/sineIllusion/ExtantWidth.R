f <- function(x) sin(x)
fprime <- function (x) cos(x)

createSine <- function(n=200, len=1, f=f, fprime=fprime) {
  x <- seq(0, 2*pi, length=n)
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
  ellx <- ell / cos(atan(abs(a*fprime(x))))
  # make this a data frame - ggplot2 doesn't do well with floating vectors
  dframe <- data.frame(x=x, xstart=x, xend=x, y=fx, ystart=ystart, yend=yend, ell=ell, ellx = ellx)
  
  dframe
}

adjY <- function(df){
  dx <- diff(df$x)
  dy <- diff(df$y)
  #aspratio <- abs(diff(range(c(df$ystart, df$yend)))/diff(range(df$x)))
  aspratio <- 1
  adj <- atan(dy/dx)
  adj <- c(adj, rev(adj)[1])
  adj <- 1/cos(atan(abs(aspratio*adj)))
  adj <- adj-mean(adj)+1
  
  df$ell <- df$ell*adj
  df$ystart <- df$y - .5*df$ell
  df$yend <- df$y + .5*df$ell
  
  df
}


getSecantSegment <- function(x1, df, f, fprime){
  secSlope <- -1/fprime(x1)
  elltemp <- df$ell[which.min(abs(df$x-x1))]
  temp <- seq(min(df$x), max(df$x), .0001)
  leftend <- temp[which.min(abs(f(temp) + elltemp/2 - secSlope*(temp-x1)))]
  rightend <- temp[which.min(abs(f(temp) - elltemp/2 - secSlope*(temp-x1)))]
  yintercept <- -secSlope*x1
  df2 <- data.frame(xstart=leftend, x=x1, xend=rightend, ystart=(yintercept + secSlope*(leftend)), y=secSlope*x1, yend = (yintercept + secSlope*(rightend)))
  df2$ell1 <- with(df2, sqrt((yend-y)^2+(xend-x)^2))
  df2$ell2 <- with(df2, sqrt((y-ystart)^2+(x-xstart)^2))  
  df2$type <- "Perceived Width"
  return(df2)
}

suppressMessages(library(ggplot2))
dframe <- createSine(n = 40, len = 1, f=f, fprime=fprime)
dframe$type <- "Segments"

secantlines <- do.call("rbind", lapply(seq(0, 2*pi, length=40), function(i) getSecantSegment(i, dframe, sin, cos)))
secantlines2 <- do.call("rbind", lapply(seq(0, 2*pi, length=40), function(i) getSecantSegment(i, adjY(dframe), sin, cos)))

library(plyr)
dframe2 <- rbind.fill(dframe, secantlines)

qplot(x=xstart, xend=xend, y = ystart, yend=yend, geom="segment", data=dframe2, colour=type) +
  scale_colour_manual(values=c("blue", "black"))+
  theme_bw() + coord_fixed(ratio=1)+ 
  geom_segment(aes(x=xstart, xend=xend, y=ystart, yend=yend), data=secantlines, colour="blue")

qplot(x=x, y=ell1, geom="line", data=subset(secantlines, type=="Perceived Width")) + 
  geom_line(aes(y=ell2), colour="blue")

qplot(x=x, xend=xend, y=ystart, yend=yend, geom="segment", data=adjY(dframe))+
  theme_bw() + coord_fixed(ratio=1)+ 
  geom_segment(aes(x=xstart, xend=xend, y=ystart, yend=yend), data=secantlines2, colour="blue") 