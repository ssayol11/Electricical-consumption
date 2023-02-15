#A
serie=ts(read.table("RenewUSA.dat", header=F),start=1990,freq=12)
plot(serie,main="Producció d'Energía renovable en USA", ylab="Trillions of BTU")
abline(v=1990:2019,col=4,lty=3)

boxplot(matrix(serie,nrow=12))

logserie=log(serie)
boxplot(matrix(logserie,nrow=12))

monthplot(logserie)

ts.plot(matrix(logserie,nrow=12))

d12logserie=diff(logserie,lag=12)
plot(d12logserie)

d1d12logserie = diff(d12logserie)
plot(d1d12logserie)
abline(h = mean(d1d12logserie),col = 2)

d1d1d12logserie = diff(d1d12logserie)
plot(d1d1d12logserie)
abline(h = mean(d1d12logserie),col = 2)

var(serie)
var(logserie)
var(d12logserie)
var(d1d12logserie)
var(d1d1d12logserie)

w_serie=d1d12logserie
par(mfrow=c(1,2))
acf(w_serie, ylim=c(-1,1), lag.max=72, col = c(2,rep(1,11)))
pacf(w_serie, ylim=c(-1,1), lag.max=72, col = c(rep(1,11),2))
par(mfrow=c(1,1))

#B
(model1 = arima(w_serie,order=c(0,0,2),seasonal = list(order=c(0,0,2), period=12)))

cat("\nT-ratios:",round(model1$coef/sqrt(diag(model1$var.coef)),2))

(model2 = arima(logserie,order=c(0,1,2),seasonal = list(order=c(0,1,2), period=12)))

cat("\nT-ratios:",round(model2$coef/sqrt(diag(model2$var.coef)),2))

(model3 = arima(w_serie,order=c(3,0,0),seasonal = list(order=c(0,0,2), period=12)))

cat("\nT-ratios:",round(model3$coef/sqrt(diag(model3$var.coef)),2))

(model4 = arima(logserie,order=c(3,1,0),seasonal = list(order=c(0,1,2), period=12)))

cat("\nT-ratios:",round(model4$coef/sqrt(diag(model4$var.coef)),2))

(model5 = arima(logserie,order=c(2,1,0),seasonal = list(order=c(0,1,2), period=12)))

cat("\nT-ratios:",round(model5$coef/sqrt(diag(model5$var.coef)),2))

#C
resi1 = resid(model2)
plot(resi1)
abline(h=0)
abline(h = c(-3,3)*sd(resi1),col = 4, lty=3)

scatter.smooth(sqrt(abs(resi1)), lpars=list(col=2))

qqnorm(resi1)
qqline(resi1,col=2,lwd=2)

hist(resi1,breaks=20, freq=FALSE)
curve(dnorm(x, mean=mean(resi1),sd=sd(resi1)), col=2, add=T)

shapiro.test(resi1)

par(mfrow=c(1,2))
acf(resi1,ylim=c(-1,1),lag.max=72,col=c(2,rep(1,12-1)),lwd=2)
pacf(resi1,ylim=c(-1,1),lag.max=72,col=c(rep(1,12-1),2),lwd=2)
par(mfrow=c(1,1))

tsdiag(model2,gof.lag=72)

Mod(polyroot(c(1,-model2$model$phi)))
Mod(polyroot(c(1,model2$model$theta)))

AIC(model2)
BIC(model2)

resi2 = resid(model4)
plot(resi2)
abline(h=0)
abline(h = c(-3,3)*sd(resi2),col = 4, lty=3)

scatter.smooth(sqrt(abs(resi2)), lpars=list(col=2))

qqnorm(resi2)
qqline(resi2,col=2,lwd=2)

hist(resi2,breaks=20, freq=FALSE)
curve(dnorm(x, mean=mean(resi2),sd=sd(resi2)), col=2, add=T)

shapiro.test(resi2)

par(mfrow=c(1,2))
acf(resi2,ylim=c(-1,1),lag.max=72,col=c(2,rep(1,12-1)),lwd=2)
pacf(resi2,ylim=c(-1,1),lag.max=72,col=c(rep(1,12-1),2),lwd=2)
par(mfrow=c(1,1))

tsdiag(model4,gof.lag=72)

Mod(polyroot(c(1,-model4$model$phi)))
Mod(polyroot(c(1,model4$model$theta)))

AIC(model4)
BIC(model4)

ultim = c(2017,12)
serie2=window(serie,end=ultim)
lnserie2 = log(serie2)

(model6 = arima(lnserie2,order=c(0,1,2),seasonal = list(order=c(0,1,2),period = 12)))

(model2)

pr = predict(model6,n.ahead = 12)
ll=exp(pr$pred-1.96*pr$se)
p=exp(pr$pred)
ul=exp(pr$pred+1.96*pr$se)
ts.plot(serie,ll,ul,p,lty=c(1,2,2,1),col=c(1,4,4,2),xlim=c(2015,2019),type="o")
abline(v=2015:2019,lty=2,col=4)

obs=window(serie,start = 2018)
(RMSE = sqrt(mean((obs-p)^2)))
(MAE = mean(abs(obs-p)))
(RMSPE = sqrt(mean(((obs-p)/p)^2)))
(MAPE = mean(abs((obs-p)/p)))

(model7 = arima(lnserie2,order=c(3,1,0),seasonal = list(order=c(0,1,2),period = 12)))

(model4)

pr = predict(model7,n.ahead = 12)
ll=exp(pr$pred-1.96*pr$se)
p=exp(pr$pred)
ul=exp(pr$pred+1.96*pr$se)
ts.plot(serie,ll,ul,p,lty=c(1,2,2,1),col=c(1,4,4,2),xlim=c(2015,2019),type="o")
abline(v=2015:2019,lty=2,col=4)

obs=window(serie,start = 2018)
(RMSE = sqrt(mean((obs-p)^2)))
(MAE = mean(abs(obs-p)))
(RMSPE = sqrt(mean(((obs-p)/p)^2)))
(MAPE = mean(abs((obs-p)/p)))

#D
pr = predict(model4,n.ahead = 12)
ll=exp(pr$pred-1.96*pr$se)
p=exp(pr$pred)
ul=exp(pr$pred+1.96*pr$se)
ts.plot(serie,ll,ul,p,lty=c(1,2,2,1),col=c(1,4,4,2),xlim=c(2010,2020))
abline(v=2010:2020,lty=2,col=4)

cat(mean(ul-ll))

#E
source("atipics2.R")
(mod.atip = outdetec(model4, dif=c(1,12), crit = 2.9, LS = T))

atipics = mod.atip$atip[order(mod.atip$atip[, 1]),]
meses = c("Ene", "Feb", "Mar", "Abr", "May", "Jun", "Jul", "Ago", "Sep", "Oct", "Nov", "Dic")
data.frame(atipics, Date = paste(meses[(atipics[,1]-1)%%12+1], start(logserie)[1]+ ((atipics[,1]-1)%/%12)), Efect = exp(atipics[,3])*100)

logserie.lin = lineal(logserie, atipics)
plot(logserie.lin,ylab="milers de vehicles", xlab="anys", main="serie linealitzada")

plot(logserie, main="serie")
lines(logserie.lin, col = 2)
legend(x="bottomright",legend=c("serie originial", "serie linealitzada"),lty = c(1, 1), col=c("black", "red"))

(model4.lin = arima(logserie.lin,order=c(3,1,2),seasonal = list(order=c(0,1,2), period=12)))

pr2 = predict(model4.lin,n.ahead = 12)
ll2=exp(pr2$pred-1.96*pr2$se)
p2=exp(pr2$pred)
ul2=exp(pr2$pred+1.96*pr2$se)
ts.plot(exp(logserie.lin),ll2,ul2,p2,lty=c(1,2,2,1),col=c(1,4,4,2),xlim=c(2010,2020))
abline(v=2010:2020,lty=2,col=4)

cat(mean(ul2-ll2))
