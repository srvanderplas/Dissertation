library(ggplot2)
library(reshape2)

correctx <- function(z, fprime, a=0, b=2*pi, shrink=0) {
  const <- integrate(function(x) abs(fprime(x)), a, b)$value
  trans <- sapply(z, function(i) integrate(function(x) abs(fprime(x)), a, i)$value*(b-a)/const + a)
  rowMeans(matrix(c(rep(z, times=shrink), trans), nrow=length(z), byrow=FALSE))
}

f <- function(x) 0.2*(x^3 - 2*x^2 - 5*x + 6)
fprime <- function(x) 0.2*(3*x^2 - 4*x - 5)
f2prime <- function(x) 0.2*(6*x - 4)

x <-  seq(-2.5, 4, .001)
data <- data.frame(x=x, y=f(x), yprime=fprime(x))
data$xtrans <- correctx(x, fprime, -2.5, 4, shrink=0)

dots <- data.frame(x=seq(-2.5, 4, .25), dots.y=-2)
dots$xtrans <- correctx(dots$x, fprime, -2.5, 4, shrink=0)

data <- merge(data, dots, all.x=TRUE)
names(data) <- c("Original", "Transformed", "y", "yprime", "dots.y")
data2 <- melt(data, id.vars=c("y", "yprime", "dots.y"), variable.name="type", value.name="x")

qplot(data=data2, x=x, y=y, geom="line") + geom_point(aes(x=x, y=dots.y)) + facet_grid(.~type)
