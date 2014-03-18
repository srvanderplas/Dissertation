library(lme4)

model2 <- lmer(data=lm.data, endweight~ (type-1) + startweight + (type-1|fingerprint))
## not necessary
# model3 <- lmer(data=lm.data, endweight~ (type-1) + startweight:(type-1) + (type-1|fingerprint))

model2.sim <- simulate(model2, nsim=1000)
res <- llply(model2.sim, function(x) model2 <- lmer(data=lm.data, x~ (type-1) + startweight + (type-1|fingerprint))
)

# extract pieces for confidence intervals
ranefs <- ldply(res, function(x) attr(VarCorr(x)$fingerprint, which="stddev"))
ranefCIs <- ldply(ranefs[,-1], function(x) quantile(x, probs=c(0.025, 0.975)))

sigma <- ldply(res, sigma)
sigmaCI <- quantile(sigma[,2], probs=c(0.025, 0.975))

fixefs <- ldply(res, fixef)
fixefs$upperX <- with(fixefs, typex+startweight)
fixefs$upperY <- with(fixefs, typey+startweight)
fixefsCIs <- ldply(fixefs[,-c(1,4)], function(x) quantile(x, probs=c(0.025, 0.975)))

###
# fixed effects table
fixefsCIs$Transformation <- NA
fixefsCIs$Transformation[grep("[y|Y]$", fixefsCIs$.id)] <- "Y"
fixefsCIs$Transformation[grep("[x|X]$", fixefsCIs$.id)] <- "X"

fixefsCIs$Threshold <- "Lower"
fixefsCIs$Threshold[grep("upper", fixefsCIs$.id)] <- "Upper"

fixefsCIs$Estimate <- c(fixef(model2)[1:2], fixef(model2)[2:1]+fixef(model2)[3])

fixefsCIs <- fixefsCIs[with(fixefsCIs,order(Transformation, Threshold)),]
fixefsCIs$interval <- with(fixefsCIs, sprintf("(%.3f, %.3f)", `2.5%`, `97.5%`))
fixefsCIs <- fixefsCIs[,c(4,5,6,7)]
fixefsCIs$Transformation <- c("X", "", "Y", "")
names(fixefsCIs)[4] <- "95% C.I."
fixeftab <- xtable(fixefsCIs)
print(fixeftab, include.rownames=FALSE, file="fixef.txt")


#####################
# random effects table

ranefCIs <- ldply(ranefs[,-1], function(x) quantile(x, probs=c(0.025, 0.975)))


ranefCIs$Groups <- "Participant"
ranefCIs$Correction <- c("Y", "X")
ranefCIs$Estimate <- attr(VarCorr(model2)$fingerprint, "stddev")
ranefCIs <- data.frame(ranefCIs)

sigmaCIs <- data.frame(t(quantile(sigma[,2], probs=c(0.025, 0.975))))
sigmaCIs$Groups <- "Residual"
sigmaCIs$Correction <- NA
sigmaCIs$Estimate <- sigma(model2)
sigmaCIs$.id <- "Residual"

ranefCIs <- rbind(ranefCIs, sigmaCIs)
ranefCIs$interval <- with(ranefCIs, sprintf("(%.3f, %.3f)", X2.5., X97.5.))
ranefCIs <- ranefCIs[,4:7]
names(ranefCIs)[4] <- "95% C.I."

ranefCIstab <- xtable(ranefCIs)
print(ranefCIstab, include.rownames=FALSE, file="ranef.txt")
