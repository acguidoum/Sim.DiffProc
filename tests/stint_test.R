options(prompt="R> ",scipen=16,digits=4,warning=FALSE, message=FALSE)
library(Sim.DiffProc)


f <- expression( exp(w-0.5*t) )
mod1 <- st.int(expr=f,type="ito",M=50,lower=0,upper=1)
print(mod1)
summary(mod1)

plot(mod1)
plot(mod1,type="n")
lines(mod1,col=2,lwd=2)
points(mod1,pch=21,col=5,cex=0.5)
bconfint(mod1)

time(mod1)
mean(mod1)
moment(mod1, center = TRUE , order = 2) ## variance
Median(mod1)
Mode(mod1)
quantile(mod1)
kurtosis(mod1)
skewness(mod1)
cv(mod1)
min(mod1)
max(mod1)
moment(mod1, center= TRUE , order = 4)
moment(mod1, center= FALSE , order = 4)

summary(mod1,at =0.2154)
mean(mod1,at =0.2154)
moment(mod1, center = TRUE , order = 2,at =0.2154) ## variance
Median(mod1,at =0.2154)
Mode(mod1,at =0.2154)
quantile(mod1,at =0.2154)
kurtosis(mod1,at =0.2154)
skewness(mod1,at =0.2154)
cv(mod1,at =0.2154)
min(mod1,at =0.2154)
max(mod1,at =0.2154)
moment(mod1, center= TRUE , order = 4,at =0.2154)
moment(mod1, center= FALSE , order = 4,at =0.2154)
bconfint(mod1,at =0.2154)

##
f <- expression( w )

mod1 <- st.int(expr=f,type="str",M=50,lower=0,upper=1)
print(mod1)
summary(mod1)

plot(mod1)
plot(mod1,type="n")
lines(mod1,col=2,lwd=2)
points(mod1,pch=21,col=5,cex=0.5)
bconfint(mod1)

time(mod1)
mean(mod1)
moment(mod1, center = TRUE , order = 2) ## variance
Median(mod1)
Mode(mod1)
quantile(mod1)
kurtosis(mod1)
skewness(mod1)
cv(mod1)
min(mod1)
max(mod1)
moment(mod1, center= TRUE , order = 4)
moment(mod1, center= FALSE , order = 4)

summary(mod1,at =0.2154)
mean(mod1,at =0.2154)
moment(mod1, center = TRUE , order = 2,at =0.2154) ## variance
Median(mod1,at =0.2154)
Mode(mod1,at =0.2154)
quantile(mod1,at =0.2154)
kurtosis(mod1,at =0.2154)
skewness(mod1,at =0.2154)
cv(mod1,at =0.2154)
min(mod1,at =0.2154)
max(mod1,at =0.2154)
moment(mod1, center= TRUE , order = 4,at =0.2154)
moment(mod1, center= FALSE , order = 4,at =0.2154)
bconfint(mod1,at =0.2154)

###

f <- expression( exp(w-0.5*t) )
mod1 <- st.int(expr=f,type="ito",M=1,lower=0,upper=1)
print(mod1)

plot(mod1)
plot(mod1,type="n")
lines(mod1,col=2,lwd=2)
points(mod1,pch=21,col=5,cex=0.5)
bconfint(mod1)

time(mod1)
mean(mod1)
moment(mod1, center = TRUE , order = 2) ## variance
Median(mod1)
quantile(mod1)
kurtosis(mod1)
skewness(mod1)
cv(mod1)
min(mod1)
max(mod1)
moment(mod1, center= TRUE , order = 4)
moment(mod1, center= FALSE , order = 4)

mean(mod1,at =0.2154)
moment(mod1, center = TRUE , order = 2,at =0.2154) ## variance
Median(mod1,at =0.2154)
quantile(mod1,at =0.2154)
kurtosis(mod1,at =0.2154)
skewness(mod1,at =0.2154)
cv(mod1,at =0.2154)
min(mod1,at =0.2154)
max(mod1,at =0.2154)
moment(mod1, center= TRUE , order = 4,at =0.2154)
moment(mod1, center= FALSE , order = 4,at =0.2154)
bconfint(mod1,at =0.2154)


##
f <- expression( w )

mod1 <- st.int(expr=f,type="str",M=1,lower=0,upper=1)
print(mod1)

plot(mod1)
plot(mod1,type="n")
lines(mod1,col=2,lwd=2)
points(mod1,pch=21,col=5,cex=0.5)
bconfint(mod1)

time(mod1)
mean(mod1)
moment(mod1, center = TRUE , order = 2) ## variance
Median(mod1)
quantile(mod1)
kurtosis(mod1)
skewness(mod1)
cv(mod1)
min(mod1)
max(mod1)
moment(mod1, center= TRUE , order = 4)
moment(mod1, center= FALSE , order = 4)


mean(mod1,at =0.2154)
moment(mod1, center = TRUE , order = 2,at =0.2154) ## variance
Median(mod1,at =0.2154)
quantile(mod1,at =0.2154)
kurtosis(mod1,at =0.2154)
skewness(mod1,at =0.2154)
cv(mod1,at =0.2154)
min(mod1,at =0.2154)
max(mod1,at =0.2154)
moment(mod1, center= TRUE , order = 4,at =0.2154)
moment(mod1, center= FALSE , order = 4,at =0.2154)
bconfint(mod1,at =0.2154)