beer = read.csv("project.csv", header = F)
data = beer[,2]
plot.ts(data, main = " Original data - Australian quaterly production of beer")

#Exploratory Analysis

plot.ts(log(data), main = " Log of the data")
plot.ts(diff(log(data)))
plot.ts(diff(log(data), lag = 4))
plot.ts(diff(diff(log(data),lag = 4)), main = " Seasonal difference + First difference of the log data")
ts.beer <- ts(data[1:152], start=c(1956,1), end=c(1993,4), frequency=4)
logbeer = log(ts.beer)

#Month Plot
monthplot(logbeer, main = " Month plot for all quarters ( On the log data)")

#Seasonal Decomposition
plot(stl(logbeer, "periodic"), main = "Seasonal Decomposition of Time Series by Loess")
 
#Fitting a non-linear trend
 library(nls2)
 dsy=logbeer-stl(logbeer,s.window="periodic")$time.series[,1]
 plot(dsy,main = " Non linear trend fitting")
 time=1:152
fit = nls(dsy ~ b0 + b1*time + b2*time^2 , start=list(b0=0, b1=1, b2 =1))
 x<-seq(from=1956,by=0.25,length=152)
 lines(x,fitted(fit))

summary(fit)

Q=factor(rep(1:4,38))
#Checking for Trend, Seasonality 
t = time(logbeer)
fit=lm(logbeer~time(logbeer)+ I(t^2)+Q)
fitr=lm(logbeer~time(logbeer) + I ( t^2)) 
plot.ts(fitr$resid) 
acf(fitr$resid)
 
anova (fit,fitr) 
summary(fit)
summary(fitr)
AIC(fit)
BIC(fit)
AIC(fitr)
BIC(fitr)
plot(fit$resid, main = "Residuals of full model fit")
plot(fitr$resid, main = " Residuals of reduced model fit")
 
library(astsa)
plot(diff(logbeer))
acf2(diff(logbeer))
plot(diff(logbeer,lag = 4))
acf2(diff(logbeer,lag =4))
plot(diff(diff(logbeer,4)))
acf2(diff(diff(logbeer,4)))
abc = diff(diff(logbeer,4))
 
fitt=lm(logbeer~Q) 
x=fitt$resid#deseasonalize the data 
plot.ts(x) 
y=time(x) 
cor.test(x,y, method=c("spearman")) 

#Checking for stationarity
library(lawstat)
acf2(abc)
bartels.test(abc)
Box.test(abc,lag=4, type = "Ljung-Box")
 
#Partitioning the data
time = time(logbeer)
Q = factor(rep(1:4,38))
timer = time[1:132]
train = logbeer[1:132]
Qr = Q[1:132]

#Least Squares Model
model1 = lm( train ~ 1+ timer + I(timer^2)+ I(timer^3)+ Qr, na.action = NULL)
newdata = data.frame(timer = time[133:152], Qr = Q[133:152])
preg = predict(model1, newdata, se.fit = TRUE)

seup=preg$fit+2*preg$se.fit
selow=preg$fit-2*preg$se.fit
time1=time[133:152]
test=logbeer[133:152]
plot(time1,test, col="red",lwd=1:2,main="Actual vs Forecast", type="l", ylim=c(5.8,6.5))
lines(time1,preg$fit, col="green", lwd=1:2)
lines(time1,seup, col="blue", lwd=1:2)
lines(time1,selow, col="blue", lwd=1:2)

#Diagnostics for Least Square
plot(model1$resid)
qqnorm(model1$resid) 
qqline(model1$resid)
acf2(model1$resid)

#Evaluation for Least Squares
mape <- function(y, yhat)
mean(abs((y - yhat)/y))
 
 library(hydroGOF)
 rmse(preg$fit,test)
 mae(preg$fit,test)
 me(preg$fit,test)
 mape(preg$fit,test)
 
#Transformed
 opreg = exp(preg$fit)
 test1= ts.beer[133:152]
 test1 = ts(test1)
 train1 = data[1:156]
 train1 = ts(train1)
 rmse(opreg,test1)
 mae(opreg,test1)
 me(opreg,test1)
 mape(opreg,test1)
 
#Holt Winters Additive
 
train=ts(train, start=1956, frequency=4)
sa=HoltWinters(train, seasonal = c("additive"))
psa=predict(sa,20)
seup_HW=psa+2*sa$SSE
selow_HW=psa-2*sa$SSE

plot(sa, psa, main = "Holt-Winters filtering, log, HW additive")
plot(time1,test, col="red",lwd=1:2,main="Actual vs Forecast", type="l", ylim = c(5.5,6.8))
lines(time1,psa, col="green", lwd=1:2)
lines(time1,seup_HW, col="blue", lwd=1:2)
lines(time1,selow_HW, col="blue", lwd=1:2)

rmse(psa,test)
mae(psa,test)
me(psa,test)
mape(psa,test)

opsa = exp(psa)
opsa = as.numeric(opsa)
rmse(opsa,test1)
mae(opsa,test1)
me(opsa,test1)
mape(opsa,test1)

#Holt Winters Multiplicative

beer <-ts.beer[1:132]
beer=ts(beer, start=1956, frequency=4)
sa1=HoltWinters(beer, seasonal = c("multiplicative"))
psa1=predict(sa1,20)
seup_HW_mult=psa1+2*sa1$SSE
selow_HW_mult=psa1-2*sa1$SSE

plot(sa1, psa1, main = "Holt-Winters filtering, log , HW multiplicative")
plot(time1,test1, col="red",lwd=1:2,main="Actual vs Forecast", type="l")
lines(time1,psa1, col="green", lwd=1:2)
lines(time1,seup_HW_mult, col="blue", lwd=1:2)
lines(time1,selow_HW_mult, col="blue", lwd=1:2)
 psa1 = as.numeric(psa1)
 rmse(psa1,test1)
 mae(psa1,test1)
 me(psa1,test1)
 mape(psa1,test1)
 