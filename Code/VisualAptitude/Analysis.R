# setwd("~/Dropbox/GraphicsGroup/TestingVisualAptitude/")

ans <- read.csv(paste(datadir, "VisualGraphicsData.csv", sep=""), stringsAsFactors=F)

key <- ans[1,]
key[,1:19] <- NA
# Visual Search
key[,20:44] <- c(17, 12, 15, 14, 8, 21, 6, 3, 21, 22, 4, 23, 9, 8, 22, 13, 2, 19, 16, 18, 7, 16, 12, 1, 7)
# Lineups 1
key[,45:64] <- c(6, 16, 13, 10, 17, 6, 16, 13, 10, 17, 6, 16, 13, 10, 17, 6, 16, 13, 10, 17)
# Card  Rotation
key[,65:224] <- c("d", "s", "s", "d", "d", "s", "d", "s",  "s", "s", "s", "d", "s", "s", "s", "s", 
                  "s", "d", "d", "d", "s", "s", "s", "d",  "s", "s", "d", "s", "d", "d", "d", "s",
                  "d", "s", "d", "d", "s", "s", "d", "s",  "s", "d", "s", "s", "s", "s", "d", "d", 
                  "s", "d", "s", "d", "d", "s", "s", "s",  "d", "d", "s", "s", "d", "s", "d", "d", 
                  "d", "d", "s", "s", "d", "s", "s", "d",  "s", "d", "d", "s", "d", "d", "s", "s",
                  "s", "s", "d", "d", "s", "s", "d", "d",  "s", "d", "d", "d", "s", "s", "s", "s",
                  "d", "d", "s", "s", "s", "s", "s", "d",  "s", "d", "s", "s", "d", "d", "d", "s",
                  "s", "s", "s", "d", "d", "d", "s", "s",  "d", "s", "s", "d", "s", "d", "d", "d",
                  "s", "s", "d", "d", "d", "d", "d", "s",  "s", "s", "s", "s", "d", "d", "s", "s",
                  "s", "s", "d", "d", "d", "d", "d", "s",  "s", "d", "d", "d", "s", "s", "s", "s")
# Lineups 2
key[,225:244] <- c(15, 15, 15,  4,  6,   15,  6, 15,  4,  6,   6,  6,  4,  4, 15,   4, 15, 15, 15, 15)

# Figure Classification
key[,245:356] <- c(2, 1, 1, 2, 1, 2, 2, 2,  2, 1, 2, 1, 2, 1, 2, 2, 
                   1, 2, 2, 1, 2, 2, 1, 2,  2, 1, 1, 2, 2, 1, 2, 1, 
                   2, 1, 2, 1, 1, 2, 1, 1,  2, 1, 2, 1, 1, 2, 2, 1, 
                   1, 2, 1, 2, 1, 2, 1, 1,  1, 2, 2, 2, 1, 1, 2, 2, 
                   3, 2, 2, 1, 2, 1, 1, 3,  1, 3, 1, 2, 3, 3, 2, 3, 
                   2, 3, 2, 1, 2, 1, 3, 3,  2, 1, 2, 2, 1, 3, 1, 3, 
                   2, 2, 1, 3, 3, 1, 1, 1,  3, 1, 2, 3, 1, 3, 2, 2)

# Lineups 3
key[,357:376] <- c(12, 5, 15, 4, 15, 1, 4, 16, 1, 19, 13, 19, 12, 11, 12, 5, 7, 12, 19, 3)

# Paper Folding
key[,377:396] <- c("a", "d", "b", "d", "b", 
                   "e", "a", "c", "e", "e", 
                   "c", "b", "a", "e", "b", 
                   "a", "e", "d", "d", "c")

library(plyr)
scored <- ddply(ans, .(id), function(j) j==key)
scored[,1:19] <- ans[,1:19]
# scored[,357:376] <- NA # can't find answers for Lineup3?
scored[,20:396] <- apply(scored[,20:396], 2, as.numeric)

library(reshape2)
longform <- melt(scored[,c(1, 20:396)], id.vars=1)
longform$testtype <- gsub("_q[0123456789]+[abcdefgh]?$", "", longform$variable)
longform$testnum <- as.numeric(gsub("[[:alpha:]_]+", "", longform$testtype))
longform$testtype <- gsub("[[:digit:]]", "", longform$testtype)
longform$penalty <- 0
longform$penalty[longform$testtype=="vis_search"] <- 1/23
longform$penalty[longform$testtype=="lineup"] <- 1/19
longform$penalty[longform$testtype=="card_rot"] <- 1
longform$penalty[longform$testtype=="folding"] <- 1/4

longform.sum <- ddply(longform, .(id, testtype, testnum), summarize, 
                      value=ifelse(is.numeric(value), sum(value, na.rm=T)-sum(!value,na.rm=T)*penalty, unique(value)), 
                      pct.answered=sum(!is.na(value))/length(value))
longform.sum$value[longform.sum$testtype=="lineup" & longform.sum$testnum==3] <- NA
ans.summary <- dcast(longform.sum, id~testtype, value.var="value", fun.aggregate = mean, na.rm=TRUE)
ans.summary <- merge(ans[,1:19], ans.summary)
ans.summary2 <- dcast(longform.sum, id~testtype+testnum, value.var = "value")
ans.summary2 <- merge(ans[,1:19], ans.summary2)
pct.ans <- dcast(longform.sum, id~testtype+testnum, value.var="pct.answered")

# 
# qplot(data=ans.summary, x=card_rot, y=lineup, geom="point")
# 
# qplot(data=ans.summary, x=folding, y=lineup, geom="point") 
# 
# qplot(data=ans.summary, x=fig_class, y=lineup, geom="point")
# 
# qplot(data=ans.summary, x=vis_search, y=lineup, geom="point")
# 

lineup.summary <- melt(ans.summary, id.vars=c(1:17, 19, 23))


# qplot(data=lineup.summary, x=value, y=lineup, geom="point") + facet_wrap(~variable, scales="free")

lineup.summary.categorical <- melt(ans.summary[,c(1:3, 12:17, 19, 23)], id.vars=c("id", "lineup"))

