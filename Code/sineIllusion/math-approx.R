library(ggplot2)
x <- seq(0, 2*pi, length=40)
f <- function(x) sin(x)
fprime <- function(x) cos(x)
f2prime <- function(x) -sin(x)

ell <- 1
qplot(x, f(x)) + geom_segment(aes(x=x, xend=x, y=f(x)-ell/2, yend=f(x)+ell/2), size=0.5) + 
  geom_line(aes(x,f(x)+ell/2)) +
  geom_line(aes(x,f(x)-ell/2)) +
# linear approximation - results in small deviations from high curvature areas  
  geom_segment(aes(x=x-0.5*ell/(fprime(x)^2+1)*fprime(x), 
                   xend=x+0.5*ell/(fprime(x)^2+1)*fprime(x), 
                   yend=f(x) - 0.5*ell/(fprime(x)^2+1), 
                   y=f(x) + 0.5*ell/(fprime(x)^2+1)), size=1, colour="steelblue")

lambda1 <- (-(fprime(x)^2 + 1) + sqrt((fprime(x)^2 + 1)^2 - 2* fprime(x)^2*f2prime(x)*ell/2))/(fprime(x)^2*f2prime(x))
lambda2 <- (-(fprime(x)^2 + 1) + sqrt((fprime(x)^2 + 1)^2 + 2* fprime(x)^2*f2prime(x)*ell/2))/(fprime(x)^2*f2prime(x))

qplot(x, f(x)) + geom_segment(aes(x=x, xend=x, y=f(x)-ell/2, yend=f(x)+ell/2), size=0.5) + 
  geom_line(aes(x,f(x)+ell/2)) +
  geom_line(aes(x,f(x)-ell/2)) +
  # quadratic approximation - almost imperceptible difference  
  geom_segment(aes(x=x+lambda1*fprime(x), 
                   xend=x+lambda2*fprime(x), 
                   yend=f(x) -lambda2, 
                   y=f(x) -lambda1), size=1, colour="red") + 
  # linear approximation - results in small deviations from high curvature areas  
  geom_segment(aes(x=x-0.5*ell/(fprime(x)^2+1)*fprime(x), 
                   xend=x+0.5*ell/(fprime(x)^2+1)*fprime(x), 
                   yend=f(x) - 0.5*ell/(fprime(x)^2+1), 
                   y=f(x) + 0.5*ell/(fprime(x)^2+1)), size=1, colour="steelblue")  
  
# now correct for it
dframe <- data.frame(x=x, y=f(x), x1=x+lambda1*fprime(x), y1=f(x) -lambda1, x2=x+lambda2*fprime(x), y2=f(x) -lambda2)
dframe$ell1 <- with(dframe, sqrt((x-x1)^2+(y-y1)^2))
dframe$ell2 <- with(dframe, sqrt((x-x2)^2+(y-y2)^2))
dframe$ellx1 <- ell/2 * 1/dframe$ell1 * ell/2
dframe$ellx2 <- ell/2 * 1/dframe$ell2 * ell/2
qplot(x, f(x)) + geom_segment(aes(x, xend=x, y=f(x)-ellx1, yend=f(x)+ellx2), data=dframe, size=1, colour="steelblue") +
  geom_segment(aes(x, xend=x, y=f(x)-ell/2, yend=f(x)+ell/2), colour="white", data=dframe) +
  coord_equal(ratio=1)
qplot(x, 0) + geom_segment(aes(x, xend=x, y=-ellx1, yend=+ellx2), data=dframe) + coord_equal(ratio=1)


#######################
x <- seq(0, 2*pi, length=100)

f <- function(x) 5*sin(x)
fprime <- function(x) 5*cos(x)
f2prime <- function(x) -5*sin(x)

ell <- 1
qplot(x, f(x)) + geom_segment(aes(x=x, xend=x, y=f(x)-ell/2, yend=f(x)+ell/2), size=0.5) + 
  geom_line(aes(x,f(x)+ell/2)) +
  geom_line(aes(x,f(x)-ell/2)) +
  # linear approximation - results in small deviations from high curvature areas  
  geom_segment(aes(x=x-0.5*ell/(fprime(x)^2+1)*fprime(x), 
                   xend=x+0.5*ell/(fprime(x)^2+1)*fprime(x), 
                   yend=f(x) - 0.5*ell/(fprime(x)^2+1), 
                   y=f(x) + 0.5*ell/(fprime(x)^2+1)), size=1, colour="steelblue") + coord_equal(ratio=1)

# need to include catches for values with fprime = 0 and f2prime = 0
lambda1 <- (-(fprime(x)^2 + 1) + sqrt((fprime(x)^2 + 1)^2 - 2* fprime(x)^2*f2prime(x)*ell/2))/(fprime(x)^2*f2prime(x))
lambda2 <- (-(fprime(x)^2 + 1) + sqrt((fprime(x)^2 + 1)^2 + 2* fprime(x)^2*f2prime(x)*ell/2))/(fprime(x)^2*f2prime(x))

qplot(x, f(x)) + #geom_segment(aes(x=x, xend=x, y=f(x)-ell/2, yend=f(x)+ell/2), size=0.5) + 
  geom_line(aes(x,f(x)+ell/2)) +
  geom_line(aes(x,f(x)-ell/2)) + coord_equal(ratio=1) +
  # quadratic approximation - almost imperceptible difference  
  geom_segment(aes(x=x+lambda1*fprime(x), 
                   xend=x+lambda2*fprime(x), 
                   yend=f(x) -lambda2, 
                   y=f(x) -lambda1), size=1, colour="red") + 
  # linear approximation - results in small deviations from high curvature areas  
  geom_segment(aes(x=x-0.5*ell/(fprime(x)^2+1)*fprime(x), 
                   xend=x+0.5*ell/(fprime(x)^2+1)*fprime(x), 
                   yend=f(x) - 0.5*ell/(fprime(x)^2+1), 
                   y=f(x) + 0.5*ell/(fprime(x)^2+1)), size=1, colour="steelblue") 
  

dframe <- data.frame(x=x, y=f(x), x1=x+lambda1*fprime(x), y1=f(x) -lambda1, x2=x+lambda2*fprime(x), y2=f(x) -lambda2)
dframe$ell1 <- with(dframe, sqrt((x-x1)^2+(y-y1)^2))
dframe$ell2 <- with(dframe, sqrt((x-x2)^2+(y-y2)^2))
dframe$ellx1 <- ell/2 * 1/dframe$ell1 * ell/2
dframe$ellx2 <- ell/2 * 1/dframe$ell2 * ell/2
qplot(x, f(x)) + geom_segment(aes(x, xend=x, y=f(x)-ellx1, yend=f(x)+ellx2), data=dframe, size=1, colour="steelblue") +
  geom_segment(aes(x, xend=x, y=f(x)-ell/2, yend=f(x)+ell/2), colour="white", data=dframe) +
  coord_equal(ratio=1)
qplot(x, 0) + geom_segment(aes(x, xend=x, y=-ellx1, yend=+ellx2), data=dframe) + coord_equal(ratio=1)
