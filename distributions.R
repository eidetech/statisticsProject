## ---- Import av bibliotek og datasett ----
library(metRology)
library(matlib)

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

