# setwd("~/Dropbox/GraphicsGroup/TestingVisualAptitude/")

# Read in data
ans <- read.csv(paste0(datadir, "VisualGraphicsData.csv"), na.strings=c(" ", ""), stringsAsFactors=F)
# Create anon. IDs that are unique
ans$id <- 1:nrow(ans)
# Fix age categories
ans$age <- gsub("21-24", "21+", gsub("25-28", "21+", ans$age))


ans$stem <- str_detect(ans$major1, "(Engineering)|(Biology)|([cC]hemistry)|(Statistics)|(Genetics)|(Environmental Science)|(Animal Science)|(Horticulture)") 
ans$stem[!is.na(ans$major2)] <- ans$stem[!is.na(ans$major2)] | str_detect(ans$major2[!is.na(ans$major2)], "(Engineering)|(Biology)|([cC]hemistry)|(Statistics)|(Genetics)|(Environmental Science)|(Animal Science)|(Horticulture)")

ans <- ans[,c(1:19, 397, 20:396)]

# Create a key that has the same structure as the answers
key <- ans[1,]
key[,1:20] <- NA
# Visual Search
key[,21:45] <- c(17, 12, 15, 14, 8, 21, 6, 3, 21, 22, 4, 23, 9, 8, 22, 13, 2, 19, 16, 18, 7, 16, 12, 1, 7)
# Lineups 1
key[,46:65] <- c(6, 16, 13, 10, 17, 6, 16, 13, 10, 17, 6, 16, 13, 10, 17, 6, 16, 13, 10, 17)
# Card  Rotation
key[,66:225] <- c("d", "s", "s", "d", "d", "s", "d", "s",  "s", "s", "s", "d", "s", "s", "s", "s", 
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
key[,226:245] <- c(15, 15, 15,  4,  6,   15,  6, 15,  4,  6,   6,  6,  4,  4, 15,   4, 15, 15, 15, 15)

# Figure Classification
key[,246:357] <- c(2, 1, 1, 2, 1, 2, 2, 2,  2, 1, 2, 1, 2, 1, 2, 2, 
                   1, 2, 2, 1, 2, 2, 1, 2,  2, 1, 1, 2, 2, 1, 2, 1, 
                   2, 1, 2, 1, 1, 2, 1, 1,  2, 1, 2, 1, 1, 2, 2, 1, 
                   1, 2, 1, 2, 1, 2, 1, 1,  1, 2, 2, 2, 1, 1, 2, 2, 
                   3, 2, 2, 1, 2, 1, 1, 3,  1, 3, 1, 2, 3, 3, 2, 3, 
                   2, 3, 2, 1, 2, 1, 3, 3,  2, 1, 2, 2, 1, 3, 1, 3, 
                   2, 2, 1, 3, 3, 1, 1, 1,  3, 1, 2, 3, 1, 3, 2, 2)

# Lineups 3
key[,358:377] <- c(12, 19, 14,  4, 13, 
                    4, 14,  7,  8, 19, 
                    8, 19, 12, 18, 12, 
                    7,  7,  6,  6,  5)

# Paper Folding
key[,378:397] <- c("a", "d", "b", "d", "b", 
                   "e", "a", "c", "e", "e", 
                   "c", "b", "a", "e", "b", 
                   "a", "e", "d", "d", "c")


keytype <- key
keytype[,46:65] <- rep(c("boxplot", "density", "scatterplot", "histogram"), each=5)
keytype[,226:245] <- c("boxplot", "boxplot", "violinplot", "violinplot", 
                       "boxplot", "violinplot", "jitter.boxplot", "jitter.boxplot", 
                       "boxplot", "violinplot", "violinplot", "bee.swarm.boxplot", 
                       "bee.swarm.boxplot", "violinplot", "violinplot", "jitter.boxplot",
                       "bee.swarm.boxplot", "jitter.boxplot", "violinplot", "bee.swarm.boxplot")
keytype[,358:377] <- c("qq.line.error", "qq.line", "qq.line.error", "qq.line", 
                   "qq.line", "qq.line.error.rotated", "qq.line.error", "qq.line.error.rotated", 
                   "qq.line.error", "qq.line.error.rotated", "qq.line", "qq.line.error", 
                   "qq.line.error.rotated", "qq.line.error.rotated", "qq.line.error", "qq.line",
                   "qq.line.error", "qq.line", "qq.line.error", "qq.line.error.rotated")

keytype <- keytype[,str_detect(names(keytype), "lineup")]
keytype <- data.frame(question=names(keytype), plot.type=as.character(keytype[1,]))

library(plyr)
# Score each answer by whether it is the same as the key. Make sure NA values stay as NA. 
scored <- ddply(ans, .(id), function(j) as.numeric(j==key) + c(0, NA)[1+is.na(j)])
names(scored) <- names(ans)
# Transfer demographic data
scored[,1:20] <- ans[,1:20]
# scored[,357:376] <- NA # can't find answers for Lineup3?
scored[,21:397] <- apply(scored[,21:397], 2, as.numeric)

# Melt scored data into long format, calculate penalties, and compute scores using ddply
library(reshape2)
library(stringr)
longform <- melt(scored[,c(1, 21:397)], id.vars=1)
longform$testtype <- gsub("_q[0123456789]+[abcdefgh]?$", "", longform$variable)
longform$testnum <- as.numeric(gsub("[[:alpha:]_]+", "", longform$testtype))
longform$testtype <- gsub("[[:digit:]]", "", longform$testtype)
longform$qnum <- as.numeric(str_replace(str_extract(longform$variable, "q\\d{1,2}"), fixed("q"), ""))
longform$penalty <- 0
longform$penalty[longform$testtype=="vis_search"] <- 1/23
longform$penalty[longform$testtype=="lineup"] <- 1/19
longform$penalty[longform$testtype=="card_rot"] <- 1
longform$penalty[longform$testtype=="folding"] <- 1/4
longform$penalty[longform$testtype=="fig_class" & longform$qnum<=8] <- 1
longform$penalty[longform$testtype=="fig_class" & longform$qnum>8] <- 1/2
 
# Compute penalty and positive points for each test type, as well as quantities needed for variance calculation
qrange <- ddply(subset(longform, id==1), .(testtype), summarize, min.score=sum(penalty), max.score=length(testtype), n=length(testtype), k=mean(1/penalty+1))
qrange$var <- with(qrange, n^2/(k-1))

# Compute penalty and positive points for each participant
longform.sum <- ddply(longform, .(id, testtype), summarize, 
                      pos.pts = ifelse(is.numeric(value), sum(value, na.rm=T), unique(value)),
                      neg.pts = ifelse(is.numeric(value), sum((1-value)*penalty, na.rm=T), unique(value)),
                      pct.answered=sum(!is.na(value))/length(value))
longform.sum <- merge(longform.sum, qrange)
# Calculate raw score
longform.sum$value <- with(longform.sum, pos.pts-neg.pts)

# calculate scores by lineup section
lineup.sec.sum <- ddply(longform, .(id, testtype, testnum), summarize, 
                         pos.pts = ifelse(is.numeric(value), sum(value, na.rm=T), unique(value)),
                         neg.pts = ifelse(is.numeric(value), sum((1-value)*penalty, na.rm=T), unique(value)), 
                         total.pos = length(value),
                         total.neg = sum(rep(1, length(value))*penalty))
lineup.sec.sum$score <- with(lineup.sec.sum, (pos.pts-neg.pts))
lineup.sec.sum2 <- ddply(subset(lineup.sec.sum, testtype!="lineup"), .(id, testtype), summarize, score=sum(score), pos.pts=sum(pos.pts), neg.pts=sum(neg.pts), total.pos=sum(total.pos), total.neg=sum(total.neg))
lineup.section.sum <- rbind.fill(lineup.sec.sum2, subset(lineup.sec.sum, testtype=="lineup"))
lineup.section.sum$testnum[is.na(lineup.section.sum$testnum)] <- ""
lineup.section.sum$testnum[nchar(lineup.section.sum$testnum)>0] <- paste0("_", lineup.section.sum$testnum[nchar(lineup.section.sum$testnum)>0])
lineup.section.sum$testtype <- paste0(lineup.section.sum$testtype, lineup.section.sum$testnum)

lineup.section.summary <- merge(ans[,1:20], dcast(lineup.section.sum, id~testtype, value.var="score", na.rm=TRUE))
rm(lineup.sec.sum2, lineup.sec.sum, lineup.section.sum)

# calculate scores by question type
longform2 <- merge(longform, keytype, by.x="variable", by.y="question", all.x=T, all.y=T)
lineup.type.sum <- ddply(longform2, .(id, testtype, plot.type, testnum), summarize, 
                         pos.pts = ifelse(is.numeric(value), sum(value, na.rm=T), unique(value)),
                         neg.pts = ifelse(is.numeric(value), sum((1-value)*penalty, na.rm=T), unique(value)), 
                         total.pos = length(value),
                         total.neg = sum(rep(1, length(value))*penalty))
lineup.type.sum$score <- with(lineup.type.sum, (pos.pts-neg.pts))
lineup.type.sum2 <- ddply(subset(lineup.type.sum, testtype!="lineup"), .(id, testtype), summarize, score=sum(score), pos.pts=sum(pos.pts), neg.pts=sum(neg.pts), total.pos=sum(total.pos), total.neg=sum(total.neg))
lineup.type.sum$score <- with(lineup.type.sum, (pos.pts-neg.pts)/total.pos)
lineup.type.sum2 <- rbind.fill(lineup.type.sum2, subset(lineup.type.sum, testtype=="lineup"))
rm(lineup.type.sum, longform2)

# Save unscaled version
# longform.sum.unscaled <- longform.sum

# Calculate test-wise means (for centering) and center scores
# tmp <- ddply(longform.sum, .(testtype), summarize, test.mean = mean(value, na.rm=T))
# longform.sum <- merge(longform.sum, tmp)
# longform.sum$value <- with(longform.sum, value - test.mean)

# Scale scores
# longform.sum$value <- with(longform.sum, value/sqrt(var))

# Alt. Scaling method by range. 
# longform.sum$value <- with(longform.sum, (pos.pts - neg.pts + min.score)/(min.score+max.score)*100)

# Cast back to wide-ish form with just test totals
ans.summary <- dcast(longform.sum, id~testtype, value.var="value", na.rm=TRUE)
ans.summary <- merge(ans[,1:20], ans.summary)
pct.ans <- dcast(longform.sum, id~testtype, value.var="pct.answered")

# Cast unscaled version
# ans.summary.unscaled <- dcast(longform.sum.unscaled, id~testtype, value.var="value", na.rm=TRUE)
# ans.summary.unscaled <- merge(ans[,1:19], ans.summary.unscaled)

# Create video game factor variable
ans.summary$vidgame_hrs_factor <- ordered(factor(rowSums(sapply(c(-1, 0, 1.99, 4.99), function(i) ans$vidgame_hrs>i)), labels=c("0", "[1, 2)", "[2, 5)", "5+"), levels=1:4))

# Create plot-able data frame with test results
lineup.summary <- melt(ans.summary, id.vars=c(1:20, 24, 26))

# Create plot-able data frame with demographic data
lineup.summary.categorical <- suppressWarnings(melt(ans.summary[,c(1:3, 12:17, 19, 20, 24, 26)], id.vars=c("id", "lineup")))
