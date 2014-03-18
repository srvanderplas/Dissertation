# Analysis of data from Amazon Turk experiment conducuted 7/22/13-7/28/13.

library(RMySQL)
library(plyr)
library(lubridate)
library(ggplot2)
library(grid)
library(reshape2)
library(doMC) 

# registerDoMC() 
# mcoptions <- list(preschedule=FALSE, set.seed=FALSE)
# getDoParWorkers() 
# 
# # ----------------- Database Access -------------------------
# con <- dbConnect(MySQL(), group="stat")
# tab <- dbReadTable(con, name="SineIllusionShiny")[-1,]
# dbDisconnect(con)
# Mode <- function(x) {
#   ux <- unique(x)
#   ux[which.max(tabulate(match(x, ux)))]
# }
# 
# diffs <- c(0, .5, .5, .5, .5, .25, .25, .25, .25, .25, .25, .25, .25, .1, .1, .05, .05, .05, .05, .05, .05, .025, .025, .025, .025, .025, .025, .025, .025, .025, .025, .025, .025, .05, .05, .05, .05, .1, .1, .1, .1, .1, .25, .25, .25, .25, .25, .25, .5, .5, .5, .5)
# wopts <- -4 +cumsum(diffs)
# 
# tab$time2 <- ymd_hms(tab$time) # convert to lubridate time
# tab$trialnum <- with(tab, q+skip)
# tab <- ddply(tab, .(iphash, fingerprint, ipid, q), transform, skipdata=max(trialnum)!=trialnum)
#   # binary variable that will only keep the data from the non-skipped question
# tab <- tab[!tab$skipdata,]
# tab <- ddply(tab, .(iphash, fingerprint, ipid), transform, ntrials=length(unique(trialnum)), ntrials.old = length(unique(q+skip))) # add number of trials
# tab <- tab[which(tab$weight>0),] 
#   # remove "initial" observation that no longer corresponds to weight in this iteration of the program. 
# tab <- subset(tab, !grepl("test", userid))
#   # remove "test" runs
# tab <- subset(tab, !grepl("76.84", ipid)) 
#   # remove trial.sum from Auburn/Nebraska City, NE, because it's full of "tests" that weren't marked as such. Whoops.
# tab$weight <- wopts[tab$weight] 
#   # convert javascript vector weight into actual meaningful weight measurement
# 
# tab$training <- (tab$time2>ymd("2013-7-21") & (tab$q+tab$skip)<2) 
#   # two "training" questions added on 7.21.13 that may be biased... test?
# 
# #------------------ Data Reshaping -------------------------
# # dataset of individual questions with each individual observation and time. 
# trial.sequence <- ddply(tab, .(iphash, fingerprint, ipid, q, skip, type, ntrials), transform,
#               weight = weight,
#               start.weight = rep(weight[which.min(time2)],length(time2)), 
#               end.weight = rep(weight[which.max(time2)],length(time2)),
#               time.diff = as.numeric(c(0, diff(time2))),
#               end.time = rep(min(time2)-max(time2),length(time2)),
#               len = rep(length(time2),length(time2)),
#               trial.time = as.numeric(time2-max(time2)),
#               seq = 1:length(time2),
#               post.training = time2>ymd("2013-7-21") & (q+skip)>=2)
# 
# # dataset of individual questions, but without each individual change
# # Includes only trials where the participant interacted with the graph at least twice. 
# trial.sum <- ddply(subset(trial.sequence, len>2 & seq>1), .(iphash, fingerprint, ipid, q, skip, type, ntrials), 
#               function(df){
#                 with(df, 
#                      data.frame(startweight = weight[which.min(time2)], 
#                                 endweight =  weight[which.max(time2)], 
#                                 modeweight =  Mode(weight), 
#                                 len = length(weight), 
#                                 unique.weights = length(unique(weight)),
#                                 time.trial = max(time2)-min(time2),
#                                 training = training[1],
#                                 post.training = min(time2)>ymd("2013-7-21") & (q+skip)[1]>=2))
#               }, .parallel=TRUE)
# 
# trial.sum$fingerid <- as.numeric(factor(trial.sum$fingerprint))
# trial.sequence$fingerid <- as.numeric(factor(trial.sequence$fingerprint))
# 
# write.csv(trial.sequence[,-which(names(trial.sequence)=="userid")], "./data/IndivTrajectory.csv")
# write.csv(trial.sum, "./data/SummaryTable.csv")


# #-----------------------Read in Data------------------------

# use if no access to the SQL database.
trial.sum <- read.csv("./data/SummaryTable.csv", row.names=1, stringsAsFactors=FALSE)
trial.sequence <- read.csv("./data/IndivTrajectory.csv", row.names=1, stringsAsFactors=FALSE)
trial.sequence$time2 <- ymd_hms(trial.sequence$time2)

# what do the different corrections look like outside of the boundary region for outliers?
f <- sin
fprime <- cos
f2prime <- function(x) -sin(x)
dframe <- createSine(40 , len=1, f=f, fprime=fprime, f2prime=f2prime)
wlow <- -2.5
whigh <- 3.5
ycorrlow <- with(dframe, data.frame(x=x, xstart=xstart, xend=xend, 
                                    yend = y + ((1-wlow)*ell/2 + (wlow)*ellx4.u),
                                    corr = "Y", level = "low"))
ycorrhigh <- with(dframe, data.frame(x=x, xstart=xstart, xend=xend, 
                                     y=y, ystart = y - ((1-whigh)*ell/2+(whigh)*ellx4.l), 
                                     yend = y + ((1-whigh)*ell/2 + (whigh)*ellx4.u),
                                     corr = "Y", level = "high"))

xcorrlow <- with(dframe, data.frame(x=correctx(x, fprime, w=wlow), 
                                    xstart=correctx(x, fprime, w=wlow),
                                    xend=correctx(x, fprime, w=wlow), 
                                    y=y, ystart=ystart, yend=yend,
                                    corr = "X", level = "low"))
xcorrhigh <- with(dframe, data.frame(x=correctx(x, fprime, w=whigh), 
                                    xstart=correctx(x, fprime, w=whigh),
                                    xend=correctx(x, fprime, w=whigh), 
                                    y=y, ystart=ystart, yend=yend,
                                    corr = "X", level = "high"))
correction.extremes <- rbind(ycorrlow,ycorrhigh,xcorrlow,xcorrhigh)
qplot(data=correction.extremes, geom="segment", x=xstart, xend=xend, y=ystart, yend=yend) + facet_grid(corr~level)


is.outlier <- function(x){
#   qs <- as.numeric(quantile(x, c(.25, .75)))
#   iqr <- diff(qs)
#   lims <- qs + c(-1, 1)*1.5*iqr
#   !(x>=lims[1] & x <= lims[2]) & !(x>=-1 & x<=2)
  !(x>=wlow & x<=whigh)
}


trial.sum <- ddply(trial.sum, .(startweight), transform, 
                   incl.startwt = startweight<=1 & startweight>=0,
                   incl.trials = ntrials>3,
                   endwt.outlier = is.outlier(endweight)) 

trial.sequence <- ddply(trial.sequence, .(start.weight), transform, 
                        incl.startwt = start.weight<=1 & start.weight>=0,
                        incl.trials = ntrials>3,
                        endwt.outlier = is.outlier(end.weight))

# picture of outlier trial time sequences. Note that in most cases, the value declined monotonically; this potentially indicates applet issues or user issues, i.e. "I'm going to click this button forever" syndrome.
ggplot() + 
  geom_point(data=trial.sum, aes(x=0, y=endweight), alpha=.1) + 
  geom_rug(data=trial.sum, aes(y=endweight), alpha=.1, sides="r") +
  geom_line(data=subset(trial.sequence, len>2 & seq>1 & ntrials>4 & trial.time>-100),
            aes(x=trial.time, y=weight, group=interaction(q+skip, fingerprint), colour=endwt.outlier), alpha=.075) + 
  facet_grid(type~., scales="free_x") + 
  xlab("Time until Trial End") + 
  ylab("Weight") + 
  geom_hline(yintercept=1) + 
  geom_hline(yintercept=0)

noutliers <- sum(trial.sum$endwt.outlier)
lm.data <- subset(trial.sum, incl.startwt & incl.trials & !endwt.outlier)

temp <- rbind.fill(cbind(trial.sum, dataset="full"), cbind(lm.data, dataset="trimmed"))
# not much has changed density wise...
ggplot(data=temp, aes(x=startweight, y=endweight)) + 
  geom_contour(aes(group=dataset, colour=dataset), stat="density2d") + 
  scale_colour_manual("Data", values=c("red", "blue"))+
  xlab("Starting Weight") + ylab("Submitted \"Correct\" Weight") + 
  facet_wrap(~type) + ggtitle("The Effect of Starting Weight on Submitted Weight")
rm("temp")


library(lme4)
library(multcomp)
model.full <- lmer(data=lm.data, 
                   endweight ~ (type-1) + post.training + training + startweight + ((type-1)|fingerprint))
summary(model.full)
# training trials are not significantly different from non-training trials, having training doesn't really matter

model.red <- lmer(data=lm.data, endweight~(type-1) + startweight + ((type-1)|fingerprint))
summary(model.red)
# still have about 2.5x variance in x as in y, but it's not *that* bad.

model <- lmer(data=lm.data, endweight~(type-1)+startweight + (1|fingerprint))
summary(model)

N <- 100
model.mcmc <- mcmcsamp(model, N, saveb=TRUE)

fixef.model <- as.data.frame(t(attr(model.mcmc, "fixef")))
names(fixef.model) <- c("x", "y", "weight")
fixef.model <- melt(fixef.model, id.vars=3)
names(fixef.model) <- c("weight", "type", "intercept")
ranef.model <- as.data.frame(t(attr(model.mcmc, "ranef")))
ranef.model <- cbind(trial = 1:N, ranef.model)
names(ranef.model) <- c("trial", paste("X", 1:125, sep=""))

sim.data <- cbind(fixef.model, rbind(ranef.model, ranef.model))

sim.data <- melt(sim.data, id.vars = 1:4)
sim.data$ub <- sim.data$weight + sim.data$intercept
sim.data$lb <- sim.data$intercept
names(sim.data)[5:6] <- c("user", "ranef")
sim.data <- melt(sim.data, id.vars = 1:6)
names(sim.data)[7:8] <- c("limit", "fixef")
sim.data$eu.level <- sim.data$ranef + sim.data$fixef

sim.data$limit <- relevel(sim.data$limit, ref="lb")
sim.data$limit <- mapvalues(sim.data$limit, from=c("lb", "ub"), to=c("from below", "from above"))
levels(sim.data$type) <- c("Transformation of X axis", "Transformation of Y values")

ints <- ddply(sim.data, .(type, limit), function(x){
  temp <- HPDinterval(t(x$fixef))
  data.frame(lb = temp[1], med = median(x$fixef), ub = temp[2])
})
names(ints) <- c("type", "limit", "xmin", "x",  "xmax")

ints.all <- ddply(ints, .(type), summarise, xmin = min(x), xmax=max(x), ymin=-Inf, ymax=Inf)

ints.users <- ddply(sim.data, .(user, type, limit), function(x){
  temp <- HPDinterval(t(x$eu.level))
  data.frame(lb = temp[1], med = median(x$eu.level), ub = temp[2])
})
names(ints.users) <- c("user", "type", "limit", "xmin", "x", "xmax")

ests <- data.frame(ests = c(fixef(model)[1:2], fixef(model)[1:2]+fixef(model)[3]), type=rep(c("Transformation of X axis", "Transformation of Y values"), 2))

ggplot() + theme_bw() + facet_grid(type~., scales="free") + 
  geom_rect(data=ints.all, aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax), 
            fill="grey75", alpha=.5) +
  geom_histogram(data=subset(sim.data, user=="X1"), aes(x=fixef, y=..density.., group=limit, fill=limit), binwidth=.01, position="identity") + 
  geom_histogram(data=subset(sim.data, user=="X1"), aes(x=fixef, y=..density.., group=limit, fill=limit, colour=limit), binwidth=.01, position="identity", show_guide=FALSE) + 
  geom_vline(xintercept=c(0,1)) +
  geom_line(data=subset(sim.data), 
            aes(x=eu.level, y=..density.., color = limit, group=interaction(user, limit)), 
            alpha=.1, stat='density', trim=TRUE) + 
  ylab("Density") + xlab("Weight w") + 
  scale_colour_manual("Approach", values=c("#B2182B", "#2166AC"),
                      breaks=c("from below", "from above"))  +
  scale_fill_manual("Approach", values=c("#d6604d", "#4393c3"),
                    breaks=c("from below", "from above")) + 
  geom_errorbarh(aes(xmin=xmin, x=x, xmax=xmax, y=32, color=limit), data=ints) + 
  theme(legend.position="bottom") + geom_point(aes(x=ests, y=32), data=ests) +
  xlim(c(-.2, 1.2)) +
  geom_point(data=ints.all, aes(x=xmax, y=0, size=factor("Range")), shape=NA, colour="grey75") +
  guides(colour=guide_legend(order=1), fill=guide_legend(order=1), size=guide_legend("Overall Acceptable Weights", override.aes=list(shape=15, size = 10), order=3))
