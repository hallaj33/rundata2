rundata <- read.csv("rundata.csv")
head(rundata)
library(statmod)
library(tweedie)
library(MASS)

#Initial Investigation
rundata <- subset(rundata, select = c(runner,avghr,maxhr,avgspeed,course))
runadam <- subset(rundata, runner=="Adam")
runquanadam <- subset(runadam, select = c(avghr,maxhr,avgspeed))
runleah <- subset(rundata, runner=="Leah")
runquanleah <- subset(runleah, select = c(avghr,maxhr,avgspeed))
my_cols <- c("#00AFBB", "#FC4E07")  #runner adam is blue, runner leah is orange

# Correlation panel
panel.cor <- function(x, y){
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(0, 1, 0, 1))
  r <- round(cor(x, y), digits=2)
  txt <- paste0("R = ", r)
  cex.cor <- 0.8/strwidth(txt)
  text(0.5, 0.5, txt, cex = cex.cor * r)
}
# Customize upper panel
upper.panel<-function(x, y){
  points(x,y, pch = 1, col = my_cols[rundata$runner])
}
# Create the plots
pairs(rundata, 
      lower.panel = panel.cor,
      upper.panel = upper.panel)
#Logarithmic relationship between avgspeed~avghr, avgspeed~maxhr for runner leah
#Positive linear relationship betweena avghr~maxhr for both runners

#Sample Mean / Variance by Runner = Adam
sapply(runquanadam,mean) 
sapply(runquanadam,var) 
table(runadam$course)
#Sample Mean / Variance by Runner = Leah
sapply(runquanleah,mean) 
sapply(runquanleah,var) 
table(runleah$course)

#Create breaks
avghrg <- cut(rundata$avghr, breaks = 6 )
maxhrg <- cut(rundata$maxhr, breaks = 6 )
#Log Mean vs Log Variance Plot
vr <- with(rundata, tapply(X = avgspeed, INDEX = list(avghrg,maxhrg,runner,course), FUN = "var" ))
mn <- with(rundata, tapply(X = avgspeed, INDEX = list(avghrg,maxhrg,runner,course), FUN = "mean" ))
vr <- as.vector(vr)
mn <- as.vector(mn)
# plot log(var) against log(mean)
par(mfrow = c(1, 1))
plot(log(vr) ~ log(mn), las = 1,  pch = 19, main = "Log Mean vs Log Variance",
     xlab = "log(group means)", ylab = "log(group variances)")
mf.lm <- lm(log(vr) ~ log(mn) )
coef(mf.lm)
abline(coef(mf.lm), lwd = 2)
# V(mu)>4, however there appears to be an influential point in bottom left that is causing slope to increase.
# As breaks increase, V(mu) trends towards 0 and then into negatives.
# Possibly the best GLM is Tweedie.
# Investigation includes Normal EDM, Gamma, Inverse Gaussian, and Tweedie.

#Normal EDM
par(mfrow=c(2,2))
run.lm <- lm(avgspeed~runner+avghr+course+maxhr, data=rundata)
## Standardized Residuals vs Transformed Fitted Values
plot(rstandard(run.lm) ~ ((fitted(run.lm ))), las=1,
     xlab="Transformed Values", main="Normal Resid vs Fitted", ylab="Standardized residuals")
## Check Linear Predictor
rW <- resid(run.lm , type = "working") # working residuals
lp <- predict(run.lm ) # linear predictor
z <- rW + lp
plot(z~lp, main="Normal Working Resid vs Predictor")
## QQ Plot for residuals
qqnorm( qr1 <- rstandard(run.lm ), las=1, main="Normal QQPlot" ); qqline( qr1 )
## Cook's distance, d
plot( cooks.distance(gamma.glm), main="Normal Cook's Distance", ylab="Cook's distance", las=1, type="h")

#Gamma GLM
gamma.glm <- glm(avgspeed~runner+avghr+course+maxhr, data=rundata, family=Gamma)
summary(gamma.glm)
par(mfrow=c(2,2))
## Standardized Residuals vs Transformed Fitted Values
plot(rstandard(gamma.glm) ~ log((fitted(gamma.glm))), las=1,
     xlab="Transformed Values", main="Gamma Stand Residuals vs Fitted", ylab="Standardized residuals")
## Check Linear Predictor
rW <- resid(gamma.glm, type = "working") # working residuals
lp <- predict(gamma.glm) # linear predictor
z <- rW + lp
plot(z~lp, main="Qamma Working Resid vs Linear Predictor")
## QQ Plot for residuals
qqnorm( qr1 <- qresid(gamma.glm), las=1, main="Gamma QQPlot"); qqline( qr1 )
## Cook's distance, d
plot( cooks.distance(gamma.glm), main="Gamma Cook's Distance", ylab="Cook's distance", las=1, type="h")


#Inverse Gaussian GLM
invg.glm <- glm(avgspeed~runner+avghr+course+maxhr, data=rundata, family=inverse.gaussian)
summary(invg.glm)
par(mfrow=c(2,2))
## Standardized Residuals vs Transformed Fitted Values
plot(rstandard(invg.glm) ~ 1/sqrt((fitted(invg.glm))), las=1,
     xlab="Transformed Values", main="Inv G Stand Residuals vs Fitted", ylab="Standardized residuals")
## Check Linear Predictor
rW <- resid(invg.glm, type = "working") # working residuals
lp <- predict(invg.glm) # linear predictor
z <- rW + lp
plot(z~lp, main="Inv G. Working Resid vs Linear Predictor")
## QQ Plot for residuals
qqnorm( qr1 <- qresid(invg.glm), ylim=c(-5,5), las=1, "Inverse Gaussian QQPlot"); qqline( qr1 )
## Cook's distance, d
plot( cooks.distance(invg.glm), main="Inv G Cook's Distance", ylab="Cook's distance", las=1, type="h")




#Tweedie GLM
par(mfrow = c(1, 1))
out <- tweedie.profile(avgspeed~runner+avghr+course+maxhr+maxhr:course, do.plot=TRUE, data=rundata)
xi.est <- out$xi.max
tweed.glm <- glm(avgspeed~maxhr+runner+avghr+course, data = rundata, 
                 family = tweedie(var.power = xi.est, link.power = 0))
summary(tweed.glm)
confint(tweed.glm)
qres1 <- qresid(tweed.glm)	# quantile  resids,  replication  1
par(mfrow = c(2, 2))
# diagnostic plots replication 1
#Resid vs Fitted
plot(qres1 ~ fitted(tweed.glm), las = 1, main="Tweedie QResid vs Fitted",
     xlab = "Fitted values", ylab = "Quantile residuals" )
#Working vs Predictor
rW <- resid(tweed.glm, type = "working") # working residuals
lp <- predict(tweed.glm) # linear predictor
z <- rW + lp
plot(z~lp, main="Tweedie Working Resid vs Linear Predictor")
#QQPlot
qqnorm(qres1,main = "Tweedie QQPlot", las = 1); qqline(qres1)
#Cook's
plot(cooks.distance(tweed.glm), type = "h", main="Tweedie Cook's Distance", 
     las = 1, ylab = "Tweedie Cook's distance, D")


#Residuals vs Individual Variables
plot(qresid(tweed.glm) ~ factor(rundata$avghr), las = 1, 
     xlab = "avghr", ylab = "Quantile residuals")
plot(qresid(tweed.glm) ~ factor(rundata$maxhr), las = 1, 
     xlab = "maxhr", ylab = "Quantile residuals")
plot(qresid(tweed.glm) ~ factor(rundata$runner), las = 1, 
     xlab = "runner", ylab = "Quantile residuals")
plot(qresid(tweed.glm) ~ factor(rundata$course), las = 1, 
     xlab = "course", ylab = "Quantile residuals")

#Interpretation
exp(coef(tweed.glm))


