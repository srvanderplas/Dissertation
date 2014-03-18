# Analysis of data from Amazon Turk experiment conducted 7/22/13-7/28/13.

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
# tab <- subset(tab, tab$time<ymd("2013-8-15"))
#   # remove trials after August 15, 2013 - not part of this experiment.
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
# trials.removed <- subset(trial.sequence, len<=2 & seq<=1)
# nrow(unique(trials.removed[,c("iphash", "fingerprint", "trialnum")])) 
# # number of trials removed for not interacting with the applet
# 
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

# #----------- Plots and Exploration - raw trial.sum --------------

# plot of directional arrows indicating start and end point by user and trial type.
qplot(data=subset(trial.sum, len>1), 
      x=startweight, xend=endweight, 
      y=fingerprint, yend=fingerprint, geom="segment", 
      arrow=arrow(length = unit(0.1,"cm")), group=q+skip, alpha=I(.2)) + 
  geom_vline(aes(xintercept=0), linetype=2) +
  geom_vline(aes(xintercept=1), linetype=2) +
  facet_grid(post.training~type)

# plot of trial trajectories for each user and trial type (not so useful now that there are a bunch of users)
qplot(data=subset(trial.sequence, len>2 & seq>1 & ntrials>6 & trial.time>-500 & 
                    !training & post.training), 
      x=trial.time, y=weight, group=q+skip, geom="line", 
      colour=factor((q+skip)%%6)) + 
  geom_point(aes(x=0, y=end.weight)) + 
  facet_grid(type~fingerprint, scales="free_x") + 
  xlab("Time until Trial End") + ylab("Weight") + 
  geom_hline(yintercept=1) + geom_hline(yintercept=0)


# plot of trial trajectories for each user
qplot(data=subset(trial.sequence, len>2 & seq>1 & ntrials>25 & 
                    trial.time>-500 & !training), 
      x=trial.time, y=weight, group=q+skip, geom="line", 
      colour=type, alpha=I(.5)) + 
  geom_point(aes(x=0, y=end.weight), alpha=.5) + 
  facet_grid(type~fingerprint, scales="free_x") + 
  xlab("Time until Trial End") + ylab("Weight") + 
  geom_hline(yintercept=1) + geom_hline(yintercept=0)

# plot of trial trajectories, for each trial type across all users. 
ggplot() + 
  geom_point(data=trial.sum, aes(x=0, y=endweight), alpha=.05) + 
  geom_rug(data=trial.sum, aes(y=endweight), alpha=.05, sides="r") +
  geom_line(data=subset(trial.sequence, len>2 & seq>1 & ntrials>4 & trial.time>-100),
            aes(x=trial.time, y=weight, group=interaction(q+skip, fingerprint)), alpha=.1) + 
  facet_grid(type~., scales="free_x") + 
  xlab("Time until Trial End") + 
  ylab("Weight") + 
  geom_hline(yintercept=1) + 
  geom_hline(yintercept=0) +
  xlim(c(-25, 1)) + 
  ylim(c(-1, 2))

# Where (approximately) did weight start?
trial.sum$startweight.cat <- factor(
  sapply(trial.sum$startweight, function(i) sum(i<=quantile(trial.sum$startweight, seq(.2, 1, .2)))),
  labels=paste(round(quantile(trial.sum$startweight, seq(0, .8, .2)), 2), 
               round(quantile(trial.sum$startweight, seq(.2, 1, .2)), 2), sep=" - "))

ggplot() + geom_density(data=trial.sum, aes(x=endweight, group=startweight.cat,
                                       colour=startweight.cat, fill=startweight.cat), alpha=I(.2)) + 
  geom_rug(data=trial.sum, aes(x=endweight), alpha=I(.1)) + 
  facet_grid(startweight.cat~type) + scale_fill_discrete("Starting Weight") + 
  scale_colour_discrete("Starting Weight") + 
  xlab("Final Weight") + ylab("Density") + ggtitle("Density of Final Weight") + xlim(c(-.5, 1.5))

qplot(data=subset(trial.sum, ntrials>10 & q>1), geom="boxplot", 
      x=factor(fingerid), y=endweight) + 
  facet_wrap(~type) + ylim(c(-1, 2)) + ggtitle("Individual boxplots")

ggplot(data=trial.sum, aes(x=startweight, y=endweight)) + 
  geom_polygon(aes(fill=..level.., group=..piece..), stat="density2d", alpha=.5) + 
  xlab("Starting Weight") + 
  ylab("Submitted \"Correct\" Weight") + 
  facet_wrap(~type)

#----------------------Data Cleaning-----------------

# qplot(x=startweight, geom="histogram", data=subset(trial.sum, !training & ntrials>3), binwidth=.1)
# 
# qplot(data=subset(trial.sum, !training & ntrials>3), x=factor(round(startweight, 3)), y=endweight, geom="violin") + 
#   geom_text(aes(x=factor(round(startweight, 3)), y=4, label=..count..), stat="bin")

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
# compute outliers for each possible start weight

trial.sequence <- ddply(trial.sequence, .(start.weight), transform, endwt.outlier = is.outlier(end.weight))

# picture of outlier trial time sequences. Note that in most cases, the value declined monotonically; this potentially indicates applet issues or user issues, i.e. "I'm going to click this button forever" syndrome.
ggplot() + 
  geom_point(data=subset(trial.sum, endwt.outlier), aes(x=0, y=endweight), alpha=.1) + 
  geom_rug(data=subset(trial.sum, endwt.outlier), aes(y=endweight), alpha=.1, sides="r") +
  geom_line(data=subset(trial.sequence, len>2 & seq>1 & ntrials>4 & trial.time>-100 & endwt.outlier),
            aes(x=trial.time, y=weight, group=interaction(q+skip, fingerprint)), alpha=.5) + 
  facet_grid(type~., scales="free_x") + 
  xlab("Time until Trial End") + 
  ylab("Weight") + 
  geom_hline(yintercept=1) + 
  geom_hline(yintercept=0)
# note that x values outside of [0,1] do not preserve the underlying function shape to any degree (i.e. concavity changes, etc.) and y values outside of [0,1] are not at all perceptually reasonable. 

# #--- Plot to compare bivar. density before/after trimming
# temp <- rbind.fill(cbind(trial.sum, dataset="full"), cbind(lm.data, dataset="trimmed"))
# # not much has changed density wise...
# ggplot(data=temp, aes(x=startweight, y=endweight)) + 
#   geom_contour(aes(group=dataset, colour=dataset), stat="density2d") + 
#   scale_colour_manual("Data", values=c("red", "blue"))+
#   xlab("Starting Weight") + ylab("Submitted \"Correct\" Weight") + 
#   facet_wrap(~type) + ggtitle("The Effect of Starting Weight on Submitted Weight")
# rm("temp")

# polygon contour plot of submitted vs starting weight for x and y
# ggplot(data=lm.data, aes(x=startweight, y=endweight)) + geom_polygon(aes(fill=..level.., group=..piece..), stat="density2d", alpha=.5) + xlab("Starting Weight") + ylab("Submitted \"Correct\" Weight") + facet_wrap(~type) + ggtitle("The Effect of Starting Weight on Submitted Weight") + xlim(c(-.2, 1.15))

# stats
nparticipants <- length(unique(trial.sum$fingerprint))
ntrials <- nrow(trial.sum)

nparttrials <- ddply(trial.sum, .(fingerprint), summarise, exclude=ntrials[1]>3)
nparttrials <- sum(nparttrials$exclude)

noutliers <- sum(trial.sum$endwt.outlier)
noutliers.x <- sum(trial.sum$endwt.outlier & trial.sum$type=="x")
noutliers.y <- sum(trial.sum$endwt.outlier & trial.sum$type=="y")


ntrials2 <- nrow(subset(trial.sum, incl.startwt & incl.trials))


lm.data <- subset(trial.sum, incl.startwt & incl.trials & !endwt.outlier, stringsAsFactors=FALSE)
lm.data$type <- relevel(factor(lm.data$type), ref="y")

trials.per.participant <- mean(ddply(lm.data, .(fingerprint), summarise, ntrials=mean(ntrials))$ntrials)


# ------------------- Participant Averages ------------------- 

end.trials <- (lm.data$startweight == 0 | lm.data$startweight == 1)

user.avg <- ddply(subset(lm.data, end.trials), # include trials starting at 0, 1
                  .(fingerprint, type), function(df){
                    avg.0 <- with(subset(df, startweight==0), mean(endweight, na.rm=TRUE))
                    avg.1 <- with(subset(df, startweight==1), mean(endweight, na.rm=TRUE))
                    return(data.frame(fingerprint=df$fingerprint[1], type=df$type[1], 
                                      avg.0=avg.0, avg.1=avg.1, ntrials=df$ntrials[1], 
                                      ntrials.sub = nrow(df)))
                  })

user.avg.all <- subset(user.avg, !is.nan(rowSums(user.avg[,3:4])))
user.avg.all$avg <- rowMeans(user.avg.all[,3:4])
user.avg.all <- ddply(user.avg.all, .(type), transform, overall.avg = mean(avg))

qplot(data=user.avg.all, x=type, y=avg, geom="violin")

ggplot(data=user.avg.all) + 
  geom_segment(aes(x=avg.0, xend=avg.1, y=fingerprint, yend=fingerprint, alpha=ntrials.sub)) + 
  scale_alpha_continuous(range=c(.5, 1)) + 
  geom_point(aes(x=avg, y=fingerprint, alpha=ntrials.sub)) +
  geom_vline(aes(xintercept=overall.avg)) +
  facet_wrap(~type)
               
ggplot(data=user.avg.all) + 
  geom_density(aes(x=avg.0), colour="blue") + 
  geom_density(aes(x=avg.1), colour="red") + 
  geom_vline(aes(xintercept=overall.avg)) + 
  facet_wrap(~type)

ggplot(data=user.avg.all) + 
  geom_density(aes(x=avg, fill=type)) + 
  geom_vline(aes(xintercept=overall.avg)) +
  geom_rug(aes(x=avg, colour=type), alpha=.9) +
  facet_grid(type~.) +
  scale_fill_manual("Transformation", values=c("#d6604d", "#4393c3")) +
  scale_colour_manual("Transformation", values=c("#d6604d", "#4393c3")) + 
  theme_bw() +
  xlim(c(-.2, 1.2))


# ---------------------- Mixed Models ---------------------------
fixed.model <- lm(data=trial.sum, endweight~startweight+type+post.training+training)
summary(fixed.model)
# start weight is significant (p<2e-16), type is not. 
# Whether or not the trial is post-training is also not significant 
#   (though it comes close when training is also included).

library(lme4)
library(multcomp)
model.full <- lmer(data=subset(trial.sum, incl.trials), 
                   endweight ~ (type-1) + post.training + training + startweight + 
                     ((type-1)|fingerprint))
summary(model.full)
# training trials are not significantly different from non-training trials, having training doesn't really matter

model.RdmOutliers <- lmer(data=subset(trial.sum, incl.trials), 
                    endweight ~ (type-1) + startweight + ((type-1)|fingerprint))
summary(model.RdmOutliers)
# removing those extra terms makes things a bit prettier, std error wise.

model.RdmNoOutliers <- lmer(data=subset(trial.sum, incl.trials & !endwt.outlier), 
                            endweight ~ (type-1) + startweight + ((type-1)|fingerprint))
summary(model.RdmNoOutliers)
# without outliers, debug model - for the first time, random effects for type have similar variances

model.Rdm <- lmer(data=lm.data, 
                  endweight ~ (type-1) + startweight + ((type-1)|fingerprint))
summary(model.Rdm)
# without outliers or starting values outside of (0,1) - high influence values without sufficient data points... should remove those as well.

model.x <- lmer(data=subset(lm.data, type=="x"), 
                endweight ~ startweight + (1|fingerprint))
summary(model.x)

model.y <- lmer(data=subset(lm.data, type=="y"),
                endweight ~ startweight + (1|fingerprint))
summary(model.y)

model <- lmer(data=lm.data, endweight~ (type-1) + startweight + (1|fingerprint))
summary(model)
N <- 1000
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
ints <- ddply(sim.data, .(type, limit), function(x){
  temp <- HPDinterval(t(x$fixef))
  data.frame(lb = temp[1], med = median(x$fixef), ub = temp[2])
})
names(ints) <- c("type", "limit", "xmin", "x",  "xmax")

ints.users <- ddply(sim.data, .(user, type, limit), function(x){
  temp <- HPDinterval(t(x$eu.level))
  data.frame(lb = temp[1], med = median(x$eu.level), ub = temp[2])
})
names(ints.users) <- c("user", "type", "limit", "xmin", "x", "xmax")


ggplot() + 
  geom_histogram(data=subset(sim.data, user=="X1"), aes(x=fixef, y=..density.., group=limit, fill=limit), binwidth=.01) +
  geom_line(data=subset(sim.data), 
               aes(x=eu.level, y=..density.., color = limit, group=interaction(user, limit)), alpha=.1, stat='density', trim=TRUE) + 
  facet_grid(type~., scales="free") + 
  ylab("Density") + xlab("Weight") + 
  scale_colour_discrete("Approach") +
  scale_fill_discrete("Approach") + 
  ggtitle("Individual and Group Level Simulations for Optimal Weight Values") + 
  geom_errorbarh(aes(xmin=xmin, x=x, xmax=xmax, y=30, color=limit), data=ints)





# 

# #-----------------------  Change+Direction  --------------------
qplot(data=lm.data, x=jitter(startweight), y=jitter(endweight), 
      xlab="Starting Weight", ylab="Final Weight", geom="point", alpha=I(.25)) + 
  geom_smooth() + ylim(c(-1,2)) + facet_wrap(~type) + theme_bw()
## Much easier to see the *amount* of change here - where amount of change is close to 0,
## correction is (approx) acceptable?


ggplot(data=lm.data, xlab="Starting Weight", ylab="Final Weight") +
  geom_point(aes(x=jitter(startweight), y=jitter(endweight)), alpha=I(.25)) + 
  geom_smooth(aes(x=startweight, y=endweight)) +
  ylim(c(-1,2)) + facet_wrap(~type) + theme_bw()

# #-----------------------Distribution of W-------------------
# logpost <- function(data, par){
#   temp <- sum(dnorm(data$endweight, mean=par[1], sd=par[2], log=TRUE))
#   temp <- temp-max(temp)*.9 # adjust by a constant so small values don't wash out.
# }
# 
# get_posterior_density <- function(data, pars){
#   temp <- sapply(1:nrow(pars), function(i) logpost(data, pars[i,]))
#   temp <- exp(temp)/sum(exp(temp))
#   trial.sum.frame(mean=pars[,1], sd=pars[,2], f=temp)
# }
# 
# #------------------------- Overall Marginals ---------------
# 
# pars <- as.matrix(expand.grid(seq(-.5, 1.2, .005), seq(.4, 1, .005)))
# 
# overall <- ddply(subset(data, ntrials>6), .(type), get_posterior_density, pars=pars)
# overall.mean <- ddply(overall[,-3], .(type, mean), summarise, f=sum(f))
# overall.mean <- ddply(overall.mean, .(type), transform, f=f/sum(f))
# #' Posterior marginal distribution over individual and std. deviation
# qplot(data=overall.mean, x=mean, y=f, geom="line", colour=type, group=type) + 
#   scale_colour_discrete("Function Type") + 
#   xlab("Mean Preferred Weighting") + ylab("Density") + theme_bw()  + theme(legend.position="bottom") 
# # ggsave("figure/fig-OverallMeansW.pdf", width=4, height=4, units="in")
# 
# 
# overall.sd <- ddply(overall[,-2], .(type, sd), summarise, f=sum(f))
# overall.sd <- ddply(overall.sd, .(type), transform, f=f/sum(f))
# 
# #' Posterior marginal distribution over individual and mean
# qplot(data=overall.sd, x=sd, y=f, geom="line", 
#       colour=factor(type), group=type) + 
#   theme_bw() + xlab("Posterior SD") + ylab("Density")
# 
# #' Posterior joint dist of mean, sd over individuals
# #' since stat_density2d won't use weights, ... improvise!
# overall.joint.sample <- sample(1:nrow(overall), size=50000, replace=TRUE, prob=overall$f)
# ggplot(data=overall[overall.joint.sample,], aes(x=mean, y=sd)) + 
#   stat_density2d(n=c(75, 40), geom="density2d", aes(colour=type)) + 
# #   facet_wrap(~type) + scale_colour_discrete(guide="none") + 
#   theme_bw() + 
#   xlab("Mean Preferred Weighting") + ylab("Std. Deviation")
# # ggsave("figure/fig-Joint2dDensityW.pdf", width=4, height=4, units="in")
# 
# #--------------------Individual Distribution of Theta ------
# 
# # get posterior density for each individual
# test <- ddply(subset(trial.sum, ntrials>12), .(fingerprint, fingerid, type), get_posterior_density, pars=pars, .parallel=TRUE)
# 
# test.mean <- ddply(test, .(fingerprint, fingerid, type, mean), summarise, f=sum(f))
# test.mean <- ddply(test.mean, .(fingerprint, fingerid, type), transform, f=f/sum(f))
# 
# participants <- dcast(ddply(trial.sum, .(fingerprint, type), summarise, 
#                             fingerid=unique(fingerid), n=length(type)), 
#                       fingerprint+fingerid~type, value.var="n")
# participants.max <- ddply(test.mean, .(fingerprint, fingerid, type), summarise, x = mean[which.max(f)], y=max(f))
# 
# ipsubset <- subset(participants, rowSums(is.na(participants))==0 & 
#                      rowSums(participants[,3:4]>12, na.rm=TRUE)==2)$fingerprint
# 
# par_labeller <- function(var, value){
#   n <- sapply(value, function(i) sum(subset(participants, fingerid%in%i)[,3:4]))
# #   value <- subset(participants, fingerprint==value)$fingerid
#   value <- paste("Participant ", as.character(value), "\n(n = ", n, ")", sep="")
#   return(value)
# }
# 
# 
# #' Plot 4 individuals who did at least 12 trials of each type 
# qplot(data=subset(test.mean, fingerprint%in%ipsubset), x=mean, y=f, group=type, colour=type, geom="line") + 
#   facet_grid(.~fingerid, labeller=par_labeller) + scale_colour_discrete("Function Type") + theme_bw() + 
#   theme(legend.position="bottom") + xlab("Preferred Weight") + ylab("Density") + 
#   geom_segment(data=subset(participants.max, fingerprint%in%ipsubset), aes(x=x, xend=x, y=0, yend=y, colour=type))
# # ggsave("figure/fig-IndivMeanAllFcnsW.pdf", width=7, height=3.5)      
# 
# 
# 
# #' Posterior mean estimates, including CI information for the individual MEAN 
# #' (i.e. not for any individual observation)
# test.post.indiv<- ddply(test.mean, .(fingerprint, fingerid, type), 
#                         function(x){
#                           ex=sum(x$mean*x$f)
#                           n=sum(trial.sum$fingerprint==x$fingerprint[1] & trial.sum$type==x$type[1])
#                           samp <- matrix(sample(x$mean, n*20, prob=x$f, replace=TRUE), ncol=20)
#                           z <- as.numeric(quantile(rowMeans(samp), c(.025, .5, .975)))
#                           trial.sum.frame(fingerprint=unique(x$fingerprint), type=unique(x$type), lb=z[1], mean = ex, median=z[2], ub=z[3], n=n)
#                         })
# 
# overall.mean.f <- ddply(test.mean, .(type, mean), summarise, f=sum(f))
# overall.mean.f <- ddply(overall.mean.f, .(type), transform, f=f/sum(f))
# 
# overall.mean.bounds <- ddply(overall.mean.f, .(type), function(x){
#   ex=sum(x$mean*x$f)
#   n=length(unique(subset(trial.sum, trial.sum$type==type)$fingerprint))
#   samp <- matrix(sample(x$mean, n*11, prob=x$f, replace=TRUE), ncol=11)
#   sample.mean = mean(samp)                          
#   sdev = sd(rowMeans(samp))
#   lb = as.numeric(quantile(rowMeans(samp), .025))
#   med = as.numeric(quantile(rowMeans(samp), .5))
#   ub = as.numeric(quantile(rowMeans(samp), .975))
#   data.frame(lb=lb, mean=sample.mean, median=med, ub=ub)
# })
# 
# 
# qplot(data=test.post.indiv,  x=lb, xend=ub, y=factor(fingerid), yend=factor(fingerid), geom="segment", colour=type) + 
#   facet_wrap(~type) + geom_point(aes(x=median), colour="black") + 
#   geom_vline(data=overall.mean.bounds, aes(xintercept=lb), linetype=3) + 
#   geom_vline(data=overall.mean.bounds, aes(xintercept=median)) + 
#   geom_vline(data=overall.mean.bounds, aes(xintercept=ub), linetype=3) + 
#   ylab("Participant ID") + xlab("Mean Preferred Weighting") + theme_bw() + theme(legend.position="none") + 
#   scale_colour_discrete("Function Type")
# # ggsave("figure/fig-CIindivMeanW.pdf", width=6, height=6, units="in")
# 
# #' Posterior estimates, including CI information for the individual observations 
# #' (i.e. not for any individual observation)
# indiv.value.bounds <- ddply(test.mean, .(fingerprint, type), function(x){
#   lb=x$mean[which.min(abs(cumsum(x$f)-.025))]
#   med=x$mean[which.min(abs(cumsum(x$f)-.5))]
#   ub=x$mean[which.min(abs(cumsum(x$f)-.975))]
#   data.frame(lb=lb, median=med, ub=ub)
# })
# 
# overall.value.bounds <- ddply(overall.mean.f, .(type), function(x){
#   xnew <- sample(x$mean, length(x$mean), prob=x$f, replace=TRUE)
#   z <- as.numeric(quantile(xnew, c(.025, .5, .975)))
#   data.frame(lb=z[1], median=z[2], ub=z[3])
# })
# # Posterior Distribution for theta without averaging over individuals
# qplot(data=overall.mean.f, x=mean, y=f, geom="line", colour=type) + 
#   xlab("Psychological Lie Factor\nEstimated Distribution for All Individuals") + 
#   theme_bw() + theme(legend.position="bottom") + scale_color_discrete("Function Type") + 
#   ylab("Density")
# 
# qplot(data=indiv.value.bounds,  x=lb, xend=ub, y=fingerprint, yend=fingerprint, geom="segment", colour=type) + 
#   facet_wrap(~type) + geom_point(aes(x=median), colour="black") + 
#   geom_vline(data=overall.value.bounds, aes(xintercept=lb), linetype=3) + 
#   geom_vline(data=overall.value.bounds, aes(xintercept=median)) + 
#   geom_vline(data=overall.value.bounds, aes(xintercept=ub), linetype=3) + 
#   ylab("Participant ID") + xlab("Lie Factor") + theme_bw() + theme(legend.position="bottom") + 
#   scale_colour_discrete("Function Type")
# 
