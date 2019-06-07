# Importing the data set
shplate <- read.csv("/Users/GuoJia/Documents/NUS/UNIL/Time Series/Data/Shanghai license plate price - Sheet3.csv", header = T)
library("astsa")
library("forecast")
library("rugarch")
library("fGarch")
View(shplate)
class(shplate)

# We decided to work with average price for our purposes
# since we are interested in how the government's policy
# has affected the general affordability of carplates
# without focusing on any other sectors of the society

# Extracting the relevant data of average prices
shplate.avg.org <- shplate[["avg.price"]]
shplate.avg2 <- shplate[["avg.price"]][-c(144:204)]

# comparing between the time series with and without 
# data from 2014 onwards
avgpr.ts <- ts(shplate.avg.org, deltat = 1/12, start = 2002)
avgpr.ts2 <- ts(shplate.avg2, deltat = 1/12, start = 2002)

# PLot of time series of average prices
par(mfrow=c(1,1))
plot(avgpr.ts, main = "Plot of All Avg Carplate Prices in Shanghai, 2002-2018",
     ylab = "Average Prices per month")
# Points that could be potential outliers
points(avgpr.ts[45], x = 2005+8/12, type = "o", col = "red" )
points(avgpr.ts[73], x = 2008+0/12, type = "o", col = "red" )
points(avgpr.ts[107], x = 2010+11/12, type = "o", col = "red")
plot(avgpr.ts2, main = "Plot of Average Prices\nwith 2014-2018 removed", 
     ylab = "Average Prices per month")

# We decided to remove the problematic data in Jan Feb 2008 (obs no. 73) 
# and Dec 2010 (obs no. 107). 
# Both are deemed to be outliers
# The first due to combined bidding in January for both months
# The second in 2010 was due to a system error. 

diff.avgpr.ts2 <- diff(avgpr.ts2)
par(mfrow = c(1,1))
plot(diff.avgpr.ts2, main = "Plot of First Difference of Average Prices", ylab = "Differenced average prices")
points(diff.avgpr.ts2[72], x = 2008+0/12, type = "o", col = "red" )
points(diff.avgpr.ts2[106], x = 2010+10/12, type = "o", col = "red")
# after differencing to remove trends, the outliers are clearly
# marked by the red circles

# hence due to their potential influence on the data, 
# we removed these points to simplify our analysis
shplate <- shplate[-c(73, 107), ]
View(shplate)
shplate.avg22 <- shplate[["avg.price"]][-c(144:204)]
avgpr.ts22 <- ts(shplate.avg22, deltat = 1/12, start = 2002)
plot(avgpr.ts22, main = "Plot of Average Carplate Prices\nwithout years 2014-19 and outliers", ylab = "Average carplate prices per month")
# We can see a clear upwards trend with some jolts (abberations)
# There is also reason to believe that there is some seasonal variation
# to the prices

# Using the stl function to decompose the data into 
# trend and seasonal components
decomppr <- stl(avgpr.ts22, s.window = "periodic")
plot(decomppr, main = "Decomposition of Average Prices Time Series")

# We start by seasonally adjusting the decomposed average prices
seasdecomppr <- seasadj(decomppr)
plot(seasdecomppr,  main = "Plot of Seasonaly Adjusted Prices")

# Analysis with double-differencing
par(mfrow = c(1,3))
plot(diff(diff(seasdecomppr)), main = "Double-Differenced Average Prices")
TSA::acf(diff(diff(seasdecomppr)), main = "ACF of Double-Differenced Avg Prices")
pacf(diff(diff(seasdecomppr)), main = "PACF of Double-Differenced Avg Prices")
# according to Duke People, the lag-1 < -0.5 and the 
# "pattern of excessive changes in sign" shows that the series is over-differenced
# hence double-differencing is not an appropriate detrending method

# Analysis with piecewise linear model
# attempt 3: removing trend with linear model

par(mfrow = c(1,1))
plot(seasdecomppr)

lin_mod <- lm(shplate.avg22 ~ (x<30) *x + (x >= 30 & x < 73)*x + (x>=73)*x)
plot(resid(lin_mod))

par(mfrow = c(1,2))
TSA::acf(resid(lin_mod))
pacf(resid(lin_mod))
# However the model will have too many regressors making it difficult
# to invert

# Using residuals after decomposing for trend using loess 
# and seasonality 

head(decomppr$time.series)
re.pr <- decomppr$time.series[,3]
plot(re.pr, main = "Residuals of Carplate Prices\nTrend & Seasonality Removed")
par(mfrow = c(1,2))
TSA::acf(re.pr)
pacf(re.pr)
# Plots seem to suggest an ARMA model

# function for aic table (consider including BIC and AICc for greater accuracy)
aic_table <- function(data, P,D, Q) {
  table <- matrix(NA, (P + 1), (Q + 1))
  for (p in 0:P) {
    for (q in 0:Q) {
      table[p + 1, q + 1] <- Arima(data, order = c(p, D, q))$aic
    }
  }
  dimnames(table) <- list(paste("<b> AR", 0:P, "</b>", sep = ""), 
                          paste("MA",0:Q, sep = ""))
  table
}
# function for bic table
bic_table <- function(data, P,D, Q) {
  table <- matrix(NA, (P + 1), (Q + 1))
  for (p in 0:P) {
    for (q in 0:Q) {
      table[p + 1, q + 1] <- Arima(data, order = c(p, D, q))$bic
    }
  }
  dimnames(table) <- list(paste("<b> AR", 0:P, "</b>", sep = ""), 
                          paste("MA",0:Q, sep = ""))
  table
}
# function for aicc table
aicc_table <- function(data, P,D, Q) {
  table <- matrix(NA, (P + 1), (Q + 1))
  for (p in 0:P) {
    for (q in 0:Q) {
      table[p + 1, q + 1] <- Arima(data, order = c(p, D, q))$aicc
    }
  }
  dimnames(table) <- list(paste("<b> AR", 0:P, "</b>", sep = ""), 
                          paste("MA",0:Q, sep = ""))
  table
}

# Finding the most parsinomious model using AIC BIC and AICc
pr_aic_table <- aic_table(re.pr, 5,0, 5)
require(knitr)
kable(pr_aic_table, digits = 2)
# suggests ARMA(2,1), ARMA(1, 3) models are most parsimonious
pr_bic_table <- bic_table(re.pr, 5,0, 5)
kable(pr_bic_table, digits = 2)
# suggests ARMA(2,1) model is most parsimonious
pr_aicc_table <- aicc_table(re.pr, 5,0, 5)
kable(pr_aicc_table, digits = 2)
# suggests ARMA(2,1), ARMA(1, 3) models are most parsimonious

# ARMA(2,1) model
model21 <- Arima(re.pr, order = c(2, 0, 1))
model21
MA_roots <- abs(polyroot(c(1, -coef(model21)[3])))
print(MA_roots) 
# MA root = 1.000067 lies on the unit circle, process might not be invertible
AR_roots <- abs(polyroot(c(1, -coef(model21)[1:2])))
print(AR_roots) 
par(mfrow = c(1,3))
plot(resid(model21), ylab = "Residual of ARMA(2,1) model", main = "ARMA(2,1) Model")
TSA::acf(resid(model21), main = "ACF of ARMA(2,1) Residuals")
pacf(resid(model21),  main = "PACF of ARMA(2,1) Residuals")
# both ACF and PACF seems relatively well behaved
# Trying GARCH model
model21.sar <- sarima(re.pr, p = 2, d = 0, q = 1)
plot(resid(model21.sar$fit))
res.sq.21 <- resid(model21.sar$fit)^2
par(mfrow = c(1,2))
TSA::acf(res.sq.21, main = "ACF of Residuals of\nARMA(2,1)-Squared")
pacf(res.sq.21, main = "PACF of Residuals of\nARMA(2,1)-Squared")
fit21 <- garchFit(~arma(2,1) + garch(1,1),re.pr,cond.dist='norm')
fit21
MA_roots <- abs(polyroot(c(1, -coef(fit21)[4])))
print(MA_roots)
AR_roots <- Mod(polyroot(c(1, -coef(fit21)[2:3])))
print(AR_roots)

summary(fit21)
par(mfrow = c(1,2))
plot(fit21, which = 13, main = "QQ-Plot of ARMA(2,1)-GARCH(1,1) Model")
cpgram(residuals(fit21))
# qqplot and cumulative periodogram shows that the model is well-behaved

# Prediction
model <- ugarchspec(variance.model = list(model = "sGARCH", 
                                          garchOrder = c(1,1)), 
                    mean.model = list(armaOrder = c(2, 1), 
                                      include.mean = FALSE),
                    distribution.model = "norm")

# Forecasting with bootstrap using ugarchboot
library("rugarch")
model_fit1 <- ugarchfit(spec = model, data = re.pr)
plot(model_fit1, which = 2) 
set.seed(505)
forc.1 <-ugarchboot(model_fit1, data = NULL, 
                    method =  "Partial", n.ahead = 58)
forc.1
par(mfrow = c(1,1))
plot(forc.1, which = 2)
plot(forc.1, which = 3)
# we can see that the variance of the data  

# we want to add the loess trend back to the data 
loessmod <- loess(shplate$avg.price ~ seq(1:length(shplate$avg.price)))
class(loessmod)
loesspred <- predict(loessmod, n.ahead = 58)
length(loesspred)

sig = sigma(forc.1@forc)
ser = fitted(forc.1@forc)
zs = cbind(t(as.data.frame(forc.1,which = "sigma", 
                           type = "summary")),  
           sig)
zr = cbind(t(as.data.frame(forc.1, which = "series", 
                           type = "summary")), 
           ser)
nrow(zr)
plot(decomppr)
# adding the loess trend and seasonality back into the data
new.predint<- cbind(zr[, 2]+loesspred[144:201]+ decomppr$time.series[,3][1:58],
                         zr[, 3]+loesspred[144:201]+ decomppr$time.series[,3][1:58],
                         zr[, 4]+loesspred[144:201]+ decomppr$time.series[,3][1:58])


predmean <- ts(new.predint[,2], start = 2014, frequency = 12)
plot(predmean, main = "Plot of Series Forecast against Actual Data"
     , ylab = "Average Carplate Price")
# 25% quantile bootstrap error band
predlow25 <- ts(new.predint[,1], start = 2014, frequency = 12)
lines(predlow25, col = "blue")
# 75% quantile bootstrap error band
predhigh25 <- ts(new.predint[,3], start = 2014, frequency = 12)
lines(predhigh25, col = "red")
data14to18 <- ts(shplate$avg.price[144:204], start=2014, frequency = 12)
lines(data14to18, col = "cyan")
?legend
legend(2014, 95000, legend=c("Bootstrap Series Forecast", 
                       "25% quantile bootstrap error", 
                       "75% quantile bootstrap error band", 
                       "Original Carplate Prices from 2014-2018"),
       col=c("black", "blue", "red", "cyan"), lty=1, cex=0.4)


# ARMA(1,3) model
model13 <- Arima(re.pr, order = c(1, 0, 3))
model13
MA_roots <- abs(polyroot(c(1, -coef(model13)[2:4])))
print(MA_roots) 
AR_roots <- abs(polyroot(c(1, -coef(model13)[1])))
print(AR_roots) 
# all roots lie outside the unit circle

plot(resid(model13))
TSA::acf(resid(model13)^2, main = "ACF of Residuals of\nARMA(1,3)-Squared")
pacf(resid(model13)^2, main = "PACF of Residuals of\nARMA(1,3)-Squared")
?qqnorm
qqnorm(resid(model13))
plot(resid(model13),pch = 2)
cpgram(resid(model13))
# model also seems reasonably well behaved

model13.sar <- sarima(re.pr, p = 1, d = 0, q = 3)
res.sq.13 <- resid(model13.sar$fit)^2
TSA::acf(res.sq.13, main = "ACF of Residuals of\nARMA(1,3)-Squared")
pacf(res.sq.13,  main = "PACF of Residuals of\nARMA(1,3)-Squared")
fit13 <- garchFit(~arma(1,3) + garch(1,1), re.pr,cond.dist='norm')

# Because of larger number of parameters, hessian matrix for 
# GARCH model is singular, and hence we are not able to work with 
# the model
