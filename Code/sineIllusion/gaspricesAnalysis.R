library(lubridate)
library(plyr)
library(ggplot2)
library(scales)
library(gridExtra)
library(quantreg)

gasprices <- read.csv("./Submission2/data/GasPrices.csv", stringsAsFactors=F)
gasprices$date <- ymd(gasprices$date)
gasprices$year <- year(gasprices$date)
gasprices <- subset(gasprices, year>=1995)
gasprices$month <- floor_date(gasprices$date, "month")+days(14)
gasprices$season <- c(rep("Winter", 2), rep("Spring", 3), rep("Summer", 3), rep("Fall", 3), "Winter")[month(gasprices$date)]
gasprices$round.year <- ifelse(month(gasprices$date)==12, gasprices$year+1, gasprices$year)
gasprices <- ddply(gasprices, .(round.year, season), transform, season.date=median(date))
gasprices <- ddply(gasprices, .(date), transform, avg=mean(price), sd=sd(price))
monthly <- ddply(gasprices, .(month), summarize, min=min(price), max=max(price), median=median(price), sd=sd(price-avg), var=var(price-avg), sd.avg=sd(avg), price=mean(price), year=year(month[1]))
monthly <- ddply(monthly, .(year), transform, avg.sd=mean(sd))
seasonal <- ddply(gasprices, .(round.year, season), summarize, date=mean(decimal_date(date)), sd=sd(price-avg), var=var(price-avg))

## 
p1 <- ggplot()  + geom_jitter(aes(date, price), data=subset(gasprices, year>=1995), size=2, alpha=.1, colour=I("grey50")) + 
  geom_line(aes(month, price), data=subset(monthly, year(month)>=1995), colour="black", size=1) + 
  ylab("Price per gallon (USD)") + theme_bw()+ xlab("date") + 
  ggtitle("Price of Gasoline in the US, 1995-2014") +
  scale_x_datetime(breaks=date_breaks("4 years"), labels = date_format("%Y"))

p2 <- ggplot( data=subset(monthly, year(month)>=1995) ) + geom_line(aes(month, sd)) + ylab(expression(sigma))  + xlab("date")+ theme_bw() + ggtitle("Standard Deviation of Monthly Prices")

grid.arrange(p1, p2, nrow=1)

##
ggplot()  + geom_boxplot(aes(decimal_date(date), price, group=decimal_date(season.date)), data=subset(gasprices, year>=1995)) +
  ylab("Price per gallon (USD)") + theme_bw()+ xlab("date") + 
  ggtitle("Price of Gasoline in the US, 1995-2014")

# Alternative to boxplot with less chaos
ggplot(data=monthly)  +
  geom_linerange(aes(x=month, ymin=min, ymax=max, y=median)) +
  ylab("Price per gallon (USD)") + theme_bw()+ xlab("date") + 
  ggtitle("Monthly Range of Gasoline Prices in the US, 1995-2014")


##
library(reshape2)
source("code/correctData.R")
model <- smooth.spline(x=decimal_date(gasprices$date), y=gasprices$price, keep.data=TRUE, nknots=230)
gas.pred <- as.data.frame(predict(model))
ggplot(data=gasprices)  +
  geom_jitter(aes(decimal_date(date), price), data=subset(gasprices, year>=1995), size=2, alpha=.1, colour=I("grey50")) + 
  ylab("Price per gallon (USD)") + theme_bw()+ xlab("Time") + 
  ggtitle("Gasoline Prices in the US, 1995-2014") +
  geom_line(data=gas.pred, aes(x=x, y=y), inherit.aes=FALSE)

gas <- gasprices[,c("date", "price")]
gas$date <- decimal_date(gas$date)
gas$pred <- f(gas$date)

gas.full <- gas
gas <- unique(gas[,c(1,3)])

f <- function(x){as.data.frame(predict(model, x=x))$y}
fprime <- function(x){as.data.frame(predict(model, x=x, deriv=1))$y}
f2prime <- function(x){as.data.frame(predict(model, x=x, deriv=2))$y}
# 
# derivs <- data.frame(x=seq(1995, 2014.3, .01))
# derivs$pred <- f(derivs$x)
# derivs$first <- fprime(derivs$x)
# derivs$second <- f2prime(derivs$x)
# derivs <- melt(derivs, id.vars=1, variable.name="Function", value.name="y")
# 
# qplot(data=derivs, x=x, y=y, colour=Function, geom="line")

asp <- .5 # image vs coordinate aspect ratio 
startyear <- 1994
gridlines <- seq(startyear+1, 2014, 1/6)
minor.axis.correction <- 
  c(correct("x", "y", data.frame(x=gridlines, y=f(gridlines)), type="x", f=f, fprime=fprime, f2prime=f2prime, weight=1)$x.correctx, 
    correct("x", "y", data.frame(x=gridlines, y=f(gridlines)), type="x", f=f, fprime=fprime, f2prime=f2prime, weight=.36)$x.correctx
  )

dots <- data.frame(x = minor.axis.correction, 
                   y = rep(-10, length(minor.axis.correction)), 
                   variable = rep(c("Corrected (w=1)", "Corrected (w=0.36)"), 
                                  each=length(gridlines)))


gas.cor <- correct("date", "price", gas, type="x", f=f, fprime=fprime, f2prime=f2prime, weight=1)
gas.cor <- merge(gas.cor, gas.full)
gas.cor <- melt(gas.cor, id.vars=c(2, 4, 5))
gas.cor$variable <- gsub("date.correctx", "Corrected (w=1)", gas.cor$variable)
gas.cor$variable <- gsub("date", "Uncorrected", gas.cor$variable)

gas.cor2 <- correct("date", "price", gas, type="x", f=f, fprime=fprime, f2prime=f2prime, weight=.36)
gas.cor2 <- merge(gas.cor2, gas.full)
gas.cor2 <- melt(gas.cor2, id.vars=c(2, 4, 5))
gas.cor2$variable <- gsub("date.correctx", "Corrected (w=0.36)", gas.cor2$variable)
gas.cor2$variable <- gsub("date", "Uncorrected", gas.cor2$variable)

gas.cor.x <- rbind(gas.cor, subset(gas.cor2, variable!="Uncorrected"))
gas.cor.x <- subset(gas.cor.x, variable!="Uncorrected")
# gas.cor.x$variable <- factor(gas.cor.x$variable, levels=c("Corrected (w=1)", "Corrected (w=0.36)", "Uncorrected"))
gas.cor.x$variable <- factor(gas.cor.x$variable, levels=c("Corrected (w=1)", "Corrected (w=0.36)"))
rm(list=c("gas.cor", "gas.cor2"))

ggplot(data=gas.cor.x) + 
  geom_point(aes(x=value, y=price, group=variable), size=1.5, shape=I(1), colour=I("grey30"), alpha=I(.75))+ 
  facet_wrap(~variable, scales="free_x", ncol=1) +
  geom_line(data=gas.cor.x, aes(x=value, y=price.fit), colour="black", size=1.5) + 
  geom_rug(data=dots, aes(x=x))+
  xlab("Time") +  
  ylab("Price of Gasoline (USD)") + 
  theme_bw()
ggsave(file="Submission2/figure/gas-x-corrected.pdf", width=12, height=7)


gas.cor <- correct("date", "price", gas.full, type="y", f=f, fprime=fprime, f2prime=f2prime, weight=1)
gas.cor <- melt(gas.cor, id.vars=c(1, 3, 5))
gas.cor$variable <- gsub("price.correcty", "Corrected (w=1)", gas.cor$variable)
gas.cor$variable <- gsub("price", "Uncorrected", gas.cor$variable)

gas.cor2 <- correct("date", "price", gas.full, type="y", f=f, fprime=fprime, f2prime=f2prime, weight=.40)
gas.cor2 <- melt(gas.cor2, id.vars=c(1, 3, 5))
gas.cor2$variable <- gsub("price.correcty", "Corrected (w=0.40)", gas.cor2$variable)
gas.cor2$variable <- gsub("price", "Uncorrected", gas.cor2$variable)

gas.cor.y <- rbind(gas.cor, subset(gas.cor2, variable!="Uncorrected"))
gas.cor.y <- subset(gas.cor.y, variable!="Uncorrected")
# gas.cor.y$variable <- factor(gas.cor.y$variable, levels=c("Corrected (w=1)", "Corrected (w=0.40)", "Uncorrected"))
gas.cor.y$variable <- factor(gas.cor.y$variable, levels=c("Corrected (w=1)", "Corrected (w=0.40)"))
rm(list=c("gas.cor", "gas.cor2"))

ggplot(data=gas.cor.y) + 
  geom_point(aes(x=date, y=value, group=variable), size=1.5, shape=I(1), colour=I("grey30"), alpha=I(.75))+ 
  facet_wrap(~variable, scales="free_y", ncol=1) +
  geom_line(data=gas.cor.y, aes(x=date, y=price.fit), colour="black", size=1.5) + 
  xlab("Time") +  
  ylab("Adj. Price of Gasoline (USD)") + 
  theme_bw()
ggsave(file="Submission2/figure/gas-y-corrected.pdf", width=12, height=7)

save.image("Submission2/data/gasprices.RData")
