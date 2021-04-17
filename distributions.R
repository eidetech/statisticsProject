## ---- Import av bibliotek og datasett ----
library(metRology)
library(matlib)
library(extraDistr)

## ---- Funksjoner ----

# Nybergs Satterthwaite metode med funksjon:
Satterthwaite<-function(sdX,sdY,nuX,nuY) 
{
  fX<-sdX^2/(nuX+1)
  fY<-sdY^2/(nuY+1)
  en<-(fX+fY)^2
  gX<-fX^2/nuX
  gY<-fY^2/nuY
  dn<-gX+gY
  nZ<-en/dn
  floor(nZ)
}

# Hypotesetester
P_H0_alphacheck = function(P_H0, alpha_hyp){
  if (P_H0 > alpha_hyp){
    print("Vi forkaster ikke H_0 til fordel for H_1: Vi vedder på H_0")
  }else{
    print("Vi forkaster H_0 til fordel for H_1: Vi vedder på H_1")
  }
}

# Plotte hypotesetest
plotHypothesis <- function(range,limit,distr,lower=TRUE){
  xVars=seq(range[1],range[2],(range[2]-range[1])/10000)
  yVars=distr(xVars)
  plot(xVars,yVars,type="l", lwd="2",col="maroon",ylim=c(0,max(yVars)))
  yVars2=yVars*((xVars<=limit)==lower)
  lines(xVars,yVars2,type="h",col="royalblue")   
  lines(xVars,yVars,type="l", lwd="2",col="maroon")
}

# Set work directory:
setwd("C:/Users/marcu/github/statisticsProject")

data0 = read.csv(file="DATA0.csv")
data1 = read.csv(file="DATA1.csv")

data0 = data0[,4]
data1 = data1[,4]

## ---- Prior og posterior for døgnet og ikke døgnet ----
# data0 = Ikke døgnet
# Nøytrale prior parametere
kappa_A_0 = 0
S_A_0 = 0 # stor sum sigma
SS_A_0 = 0
nu_A_0 = -1
C_A_0 = 0

(n_A = length(data0)) # antall elementer
(Sx_A = sum(data0)) # sum av elementer
(Sxx_A = sum(data0^2)) # sum av kvadratet av målingene
(xbar_A = Sx_A/n_A) # gjennomsnitt av n elementer
(xdev_A = data0 - xbar_A) # data0 sine avvik fra gjennomsnittet
(SSx_A = sum(xdev_A^2)) # Sum av kvadratisk avvik

# Posterior hyperparametere
(kappa_A_1 = kappa_A_0 + n_A)
(S_A_1 = S_A_0 + Sx_A)
(nu_A_1 = nu_A_0 + n_A)
(C_A_1 = C_A_0 + Sxx_A)
(m_A_1 = S_A_1/kappa_A_1) 
(SS1_A = C_A_1 - kappa_A_1*m_A_1^2)

(s1_A = sqrt(SS1_A/nu_A_1))

# data1 = Døgnet
# Nøytrale prior parametere
kappa_B_0 = 0
S_B_0 = 0 #sigma
SS_B_0 = 0
nu_B_0 = -1
C_B_0 = 0

(n_B = length(data1)) # antall elementer
(Sx_B = sum(data1)) # sum av elementer
(Sxx_B = sum(data1^2)) # sum av kvadratet av målingene
(xbar_B = Sx_B/n_B) # gjennomsnitt av n elementer
(xdev_B = data1 - xbar_B) # data0 sine avvik fra gjennomsnittet
(SSx_B = sum(xdev_B^2)) # Sum av kvadratisk avvik

# Posterior hyperparametere
(kappa_B_1 = kappa_B_0 + n_B)
(S_B_1 = S_B_0 + Sx_B)
(nu_B_1 = nu_B_0 + n_B)
(C_B_1 = C_B_0 + Sxx_B)
(m_B_1 = S_B_1/kappa_B_1) 
(SS1_B = C_B_1 - kappa_B_1*m_B_1^2)

(s1_B = sqrt(SS1_B/nu_B_1)) # liten stress sigma

## ---- Plotting ----
# mu fordeling for data0 (ikke døgnet)
xVals = seq(3, max(data0), 0.01)
yVals1=dt.scaled(xVals,nu_A_1,m_A_1,s1_A*sqrt(1/kappa_A_1))  # Samme, men med y-verdier innenfor (0,1)
plot(xVals,yVals1,type="l",col="royalblue", lwd="2", xlab="Rundetid (minutt)", ylab="", main="mu fordelinger")

# mu fordeling for data1 (døgnet)
yVals2=dt.scaled(xVals,nu_B_1,m_B_1,s1_B*sqrt(1/kappa_B_1))  # Samme, men med y-verdier innenfor (0,1)
lines(xVals,yVals2,type="l",col="maroon", lwd="2", xlab="Rundetid (minutt)", ylab="", main="mu fordelinger")
legend(4.9, 3.3, legend=c("Ikke døgnet - mu", "Døgnet - mu"),
       col=c("royalblue", "maroon"), lty=1, cex=0.8)

# X+ fordeling for data0 (ikke døgnet)
xVals = seq(0, max(data1), 0.01)
yVals3=dt.scaled(xVals,nu_A_1,m_A_1,s1_A*sqrt(1+1/kappa_A_1))  # Samme, men med y-verdier innenfor (0,1)
plot(xVals,yVals3,type="l",col="royalblue", lwd="2", xlab="Rundetid (minutt)", ylab="", main="X+ fordelinger")

# X+ fordeling for data1 (døgnet)
yVals4=dt.scaled(xVals,nu_B_1,m_B_1,s1_B*sqrt(1+1/kappa_B_1))  # Samme, men med y-verdier innenfor (0,1)
lines(xVals,yVals4,type="l",col="maroon", lwd="2", xlab="Rundetid (minutt)", ylab="", main="X+ fordelinger")
legend(5.1, 0.66, legend=c("Ikke døgnet - X+", "Døgnet - X+"),
       col=c("royalblue", "maroon"), lty=1, cex=0.8)

## ---- mu og X+ for ikke døgnet ----
# mu fordeling for data0 (ikke døgnet)
xVals = seq(2, max(data0), 0.01)
yVals1=dt.scaled(xVals,nu_A_1,m_A_1,s1_A*sqrt(1/kappa_A_1))  # Samme, men med y-verdier innenfor (0,1)
plot(xVals,yVals1,type="l",col="royalblue", lwd="2", xlab="Rundetid (minutt)", ylab="", main="mu og X+ for ikke døgnet")

# X+ fordeling for data0 (ikke døgnet)
yVals3=dt.scaled(xVals,nu_A_1,m_A_1,s1_A*sqrt(1+1/kappa_A_1))  # Samme, men med y-verdier innenfor (0,1)
lines(xVals,yVals3,type="l",col="maroon", lwd="2", xlab="Rundetid (minutt)", ylab="", main="mu og X+ for ikke døgnet")
legend(4.6, 3.2, legend=c("Ikke døgnet - mu", "Ikke døgnet - X+"),
       col=c("royalblue", "maroon"), lty=1, cex=0.8)

## ---- mu og X+ for døgnet ----
# mu fordeling for data1 (døgnet)
xVals = seq(1, max(data1), 0.01)
yVals2=dt.scaled(xVals,nu_B_1,m_B_1,s1_B*sqrt(1/kappa_B_1))  # Samme, men med y-verdier innenfor (0,1)
plot(xVals,yVals2,type="l",col="royalblue", lwd="2", xlab="Rundetid (minutt)", ylab="", main="mu og X+ for døgnet")

# X+ fordeling for data1 (døgnet)
yVals4=dt.scaled(xVals,nu_B_1,m_B_1,s1_B*sqrt(1+1/kappa_B_1))  # Samme, men med y-verdier innenfor (0,1)
lines(xVals,yVals4,type="l",col="maroon", lwd="2", xlab="Rundetid (minutt)", ylab="", main="mu og X+ for døgnet")
legend(5.6, 1.9, legend=c("Døgnet - mu", "Døgnet - X+"),
       col=c("royalblue", "maroon"), lty=1, cex=0.8)

## ---- tau fordelinger ----
# tau fordeling for døgnet:
xVals = seq(0, max(data0), 0.01)
yVals5=dgamma(xVals, nu_B_1/2, SS1_B/2)  # Samme, men med y-verdier innenfor (0,1)
plot(xVals,yVals5,type="l",col="maroon", lwd="2", xlab="Rundetid (minutt)", ylab="", main="tau fordelinger")

# tau fordeling for ikke døgnet:
yVals6=dgamma(xVals, nu_A_1/2, SS1_A/2)  # Samme, men med y-verdier innenfor (0,1)
lines(xVals,yVals6,type="l",col="royalblue", lwd="2", xlab="Rundetid (minutt)", ylab="", main="tau fordelinger")
legend(4.3, 1.2, legend=c("Døgnet", "Ikke døgnet"),
       col=c("royalblue", "maroon"), lty=1, cex=0.8)

## ---- Hypotesetesting ----
## ---- Hypotesetesting - Sammenlikne muA og muB mot hverandre ----
alpha = 0.05

(P_H0_muAB = pnorm(0, m_A_1-m_B_1, sqrt((s1_A^2/kappa_A_1)+(s1_B^2/kappa_B_1))))
# muA <= muB
P_H0_alphacheck(P_H0_muAB, alpha)

## Normalfordeling til å plotte hypotesen.
theDistr = function(x) {
  dnorm(x, m_A_1-m_B_1, sqrt((s1_A^2/kappa_A_1)+(s1_B^2/kappa_B_1)))
}

lower = -2
upper = 2
theLimit = 0
plotHypothesis(c(lower,upper),theLimit,theDistr,TRUE)
# Vi ser en klar overvekt i H0.


## ---- Hypotesetesting - sammenlikne muA/muB mot fast verdi ----
# -- Grunnlag for beregning: --
# mu0 = fast verdi vi vil sjekke mot
# H_1_muA: muA > mu0
# H_0_muA: muA <= mu0
# alpha = 0.05

# P(H_0_muA = P(mu_A <= mu0))
# P_H0_muA = dt.scaled(mu0,nu_A_1,m_A_1,s1_A*sqrt(1/kappa_A_1)))

alpha = 0.05 # Velger signifikans = 0.05
mu0 = 3.8 # Leser av mu fordelingene og ser at ca 3.8 er forventningsverdi.
(P_H0_muA = pt.scaled(mu0,nu_A_1,m_A_1,s1_A*sqrt(1/kappa_A_1)))
(P_H0_muB = pt.scaled(mu0,nu_B_1,m_B_1,s1_B*sqrt(1/kappa_B_1)))

P_H0_alphacheck(P_H0_muA, alpha)
P_H0_alphacheck(P_H0_muB, alpha)

## t fordelinger til å plotte hypotesene.
theDistr1 = function(x) {
  dt.scaled(x,nu_A_1,m_A_1,s1_A*sqrt(1/kappa_A_1))
}
theDistr2 = function(x) {
  dt.scaled(x,nu_B_1,m_B_1,s1_B*sqrt(1/kappa_B_1))
}

lower = 3
upper = 5
theLimit = mu0
plotHypothesis(c(lower,upper),theLimit,theDistr1,TRUE)
plotHypothesis(c(lower,upper),theLimit,theDistr2,TRUE)

