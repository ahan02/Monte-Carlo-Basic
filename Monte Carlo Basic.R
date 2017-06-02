#####Randomization Analysis of matched-pair experiments##

##null hypothesis of no difference, the mean differential is expected to be 0###
diff.dat<-c() #import the data
n<-length(diff.dat)
t0<-(2*length(diff.dat[diff.dat<0])-n)/n
set.seed(1237)
t<-rep(0,10000)
for (i in 1:10000) {
  m<-rbinom(1,15,0.5)
  t[i]<-(2*m-n)/n
}
#calculate p-value#
pvalue<-sum(abs(t)>=abs(t0))/10000

####Permulation Analysis of comparisons of means studies###
###Compare two group of patients responses, null hypothesis of no differences###
a<-c()#Import data, responses cases were binary, y=1, no=0
b<-c()
c<-c(a,b)

diff.observed<-mean(a)-mean(b)
n<-10000

diff.perm<-rep(0,n)
for (i in 1:n){
  temp<-sample(c,size=length(c),replace=FALSE)
  a.1<-temp[1:10]
  b.1<-temp[11:20]
  diff.perm[i]<-mean(a.1)-mean(b.1)
}

p.value.perm<-sum(abs(diff.perm)>abs(diff.observed))/n
print(p.value.perm)




###Comparison of bias correction from bootstraps and jacknives###

x.10 <- rnorm(10,0,1) #generate random numbers from normal dist#
var.x.10 <- sum((x.10-mean(x.10))^2)/10

#parametric boostrap#
var.x.10.pboot <- rep(0,10000)
for (i in 1:10000) {
  x.10.pboot <- rnorm(10,mean(x.10),var.x.10)
  var.x.10.pboot[i] <- sum((x.10.pboot-mean(x.10.pboot))^2)/10
}
bias.x.10.pboot <- mean(var.x.10.pboot) - var.x.10
var.x.10.pboot.unbiased <- var.x.10-bias.x.10.pboot

#nonparametric bootstrap#

var.x.10.npboot<-rep(0,10000)
for ( i in 1:10000) {
  x.10.npboot<-sample(x.10,10,replace=TRUE)
  var.x.10.npboot[i]<-sum((x.10.npboot-mean(x.10.npboot))^2)/10  
}
bias.x.10.npboot<-mean(var.x.10.npboot)-var.x.10
var.x.10.npboot.unbiased<-var.x.10-bias.x.10.npboot


#Jacknife#

var.x.10.jack<-rep(0,10)
for ( i in 1: 10){
  var.x.10.omit<-sum((x.10[-i]-mean(x.10[-i]))^2/(10-1))
  var.x.10.jack[i]<-var.x.10-(10-1)*(var.x.10-var.x.10.omit)
}
bias.x.10.jack<-mean(var.x.10.jack)-var.x.10
var.x.10.jack.unbiased<-var.x.10-bias.x.10.jack


#bias correction#

var.x.10.npboot.unbiased<-rep(0,100)
var.x.10.jack.unbiased<-rep(0,100)
var.x.10.vec <- rep(0,100)

for (j in 1:100){
  x.10 <- rnorm(10,0,1)
  var.x.10 <- sum((x.10-mean(x.10))^2)/10
  var.x.10.npboot<-rep(0,10000)
  for ( i in 1:10000) {
    x.10.npboot<-sample(x.10,10,replace=TRUE)
    var.x.10.npboot[i]<-sum((x.10.npboot-mean(x.10.npboot))^2)/10    
  }
  bias.x.10.npboot<-mean(var.x.10.npboot)-var.x.10
  
  var.x.10.jack<-rep(0,10)
  for ( i in 1: 10){
    var.x.10.omit<-sum((x.10[-i]-mean(x.10[-i]))^2/(10-1))
    var.x.10.jack[i]<-var.x.10-(10-1)*(var.x.10-var.x.10.omit)    
  }
  bias.x.10.jack<-mean(var.x.10.jack)-var.x.10
  var.x.10.vec[j] <- var.x.10
  var.x.10.jack.unbiased[j]<-var.x.10-bias.x.10.jack
  var.x.10.npboot.unbiased[j]<-var.x.10-bias.x.10.npboot
}

print(paste(mean(var.x.10.vec),mean(var.x.10.npboot.unbiased),mean(var.x.10.jack.unbiased)))
##plot three boostrap##
plot(sort(var.x.10.vec),(1:100)/101,type='s',col='blue',lwd=2,xlab='Estimate',ylab='Frequence')
lines(sort(var.x.10.npboot.unbiased),(1:100)/101,type='s',col='red',lwd=2)
lines(sort(var.x.10.jack.unbiased),(1:100)/101,type='s',col='green',lwd=2)

print(paste(mean(var.x.10.vec),mean(var.x.10.npboot.unbiased),mean(var.x.10.jack.unbiased)))

####Check boostrapping Confidence Intervals for logistic model####
###Given coveriate values for dependent cast/control observations##
##confidence interval for covariate effect##
setwd("Directory")
log.data<- read.table('data.txt',header=TRUE)
log.mod<-glm(Obs~1+Cov,data=log.data)
summary(log.mod)
beta.0<-summary(log.mod)$coef[1,1]
beta.1<-summary(log.mod)$coef[2,1]
sd.beta1<-summary(log.mod)$coef[2,2]

print(paste(beta.1-(1.96*sd.beta1),beta.1,beta.1+(1.96*sd.beta1)))
#paramateric boostrap#
p<-exp(beta.0 +(beta.1*log.data$Cov))/(1+exp(beta.0+(beta.1*log.data$Cov)))
beta.1.boot<-rep(0,10000)
for ( i in 1:10000){
  Obs.1<-rep(0,nrow(log.data))
  for ( j in 1:nrow(log.data)) {
    Obs.1[j]<-rbinom(1,size=1,prob=p[j])
    mod.1<-glm(Obs.1~1+log.data$Cov, data=log.data)
    beta.1.boot[i]<-summary(mod.1)$coef[2,1]
  }
}
print(paste(sort(beta.1.boot)[0.025*10000],mean(beta.1.boot), sort(beta.1.boot)[0.975*10000]))

#nonparamteric bootstrap#

beta.1.npboot<-rep(0,10000)
for ( i in 1:10000){
  ind<-sample(c(1:nrow(log.data)),size=nrow(log.data),replace=TRUE)
  Cov.1<-log.data$Cov[ind]
  Obs.1<-log.data$Obs[ind]
  mod.1<-glm(Obs.1~1+Cov.1, data=log.data)
  beta.1.npboot[i]<-summary(mod.1)$coef[2,1]       
}

print(paste(sort(beta.1.npboot)[0.025*10000],mean(beta.1.npboot), sort(beta.1.npboot)[0.975*10000]))

#Jackknife bootstrap#

beta.1.jack<-rep(0,20)
for ( i in 1:20){
  mod.omitted<-glm(Obs[-i]~1+Cov[-i], data=log.data)
  beta.1.jack[i]<-beta.1 - (20-1)*(beta.1-mod.omitted$coef[2])       
}
print(paste(beta.1-(1.96*sqrt(var(beta.1.jack)/20)),beta.1,beta.1+(1.96*sqrt(var(beta.1.jack)/20))))





