rm(list=ls(all=TRUE)) # erasure all objects
system('./LBPA -ind lbpa.dat')  # for model running

source('read.admb.R')
source('read.admbFit.R')

data <-read.rep('For_R.rep')
nages<-10 #this depends on the specie

nages <- length(data$Probability_of_length[,1])
age <- seq(1,nages) 
BinLen <- data$Length_bins
NObsFre <- length(data$Observed_frequencies[,1]) #this numbers of observations depend of the own data
ObsFre <- data$Observed_frequencies
PredFre <- data$Predicted_frequency
CatchLFre<- data$Catch_length_frequency
ProbLen <- data$Probability_of_length

Fcr <- data$`F_Y/R_SSB/R`[,1]
YPR <- data$`F_Y/R_SSB/R`[,2]
SPR <- data$`F_Y/R_SSB/R`[,3]
puntos=data$`F/Ftar_SPR_SPRtar`
results=data$F_L50_slope_a0_cv_Lr_Ftar



#---------FIGURA 1-------------------------------------

#par(mfrow = c(2, 1))

# 1
plot(BinLen,
     PredFre/max(PredFre),type="l",
     ylab="Frequency",
     xlab="Length",
     ylim=c(0, 1),
     main="Model fit",xlim = c(min(BinLen),1.05*max(BinLen)))


for (i in 1:NObsFre) {
  lines(BinLen, ObsFre[i,]/max(ObsFre[i,]), type="p",  cex=0.6, col="black", pch = 16)
}

lines(BinLen,
     PredFre/max(PredFre),type="l",col="red",
     lwd=1,
     ylim=c(0, 1)) 


# 2

plot(ObsFre[1,],PredFre,  cex=0.6, col="black", pch = 16,ylab="Predicted", xlab="Observed", main="Model fit")

for (i in 1:NObsFre) {
  lines(ObsFre[i,],PredFre, type="p", cex=0.6, col="black", pch = 16)
}
lines(PredFre,PredFre, type="l", col="red", lwd=1)

#---------FIGURA 2-------------------------------------

#par(mfrow = c(2, 2))


# 1

L1 <- (data$Length_frequency_of_exploitable_population_current_target_unfished[1,])
L2 <- (data$Length_frequency_of_exploitable_population_current_target_unfished[2,])
L3 <- (data$Length_frequency_of_exploitable_population_current_target_unfished[3,])


plot(BinLen,L3, ylab="Frequency", xlab="Length", type="l", lwd=1, col="green", xlim = c(min(BinLen),1.05*max(BinLen)),
     main="Length frequencies")
polygon(c(BinLen[1], BinLen,BinLen[length(L3)]),c(0,L3,0),col = rgb(0.46, 0.93, 0.78, alpha = 0.5))

lines(BinLen,L2)
polygon(c(BinLen[1], BinLen,BinLen[length(L3)]),c(0,L2,0), col = rgb(1, 0.27, 0, alpha = 0.8))

lines(BinLen,L1,type="l", lwd=2)
legend(x = "topright",
       legend=c("Current","Target","Unfished"),
       bty = "n",
       col=c(1,2,3),
       lwd=3)

# 2
S1 <- (data$Selectivity_and_maturity_at_length[1,])
S2 <- (data$Selectivity_and_maturity_at_length[2,])

plot(BinLen,S1, ylab="Proportion", xlab="Length", type="l", lwd=3,  xlim = c(min(BinLen),1.05*max(BinLen)),
     main="Ogives")
lines(BinLen,S2, lwd=3, col="green", xlim = c(min(BinLen),1.1*max(BinLen)))
legend(x = "topright",legend=c("Selectivity","Maturity"), bty = "n", col=c(1,3), lwd=3)
grid(nx = NULL, ny = NULL, lty = 2, col = "gray",lwd = 1)    

# 2b -------------------------------------------
#par(mfrow = c(1, 1))

plot(BinLen,
     PredFre/max(PredFre),type="l",
     ylab="Frequency",
     xlab="Length",
     ylim=c(0, 1),
     main="Model fit",xlim = c(min(BinLen),1.05*max(BinLen)))


for (i in 1:NObsFre) {
  lines(BinLen, ObsFre[i,]/max(ObsFre[i,]), type="l",  cex=0.6, col="red", pch = 16)
}
lines(BinLen,PredFre/max(PredFre),type="l",lwd=2)
     



# 3
Fcr_est=results[1]
Ftar=results[7]
SPR_est=puntos[2]


plot(Fcr,YPR/max(YPR),type="l", ylab="YPR, SPR", xlab="Fishing mortality", lwd=3, col="green",
     main="Per-recruit analysis")
lines(Fcr,SPR/max(SPR),type="l", ylab="YPR, SPR", xlab="F/Ftar", lwd=3, col="red")
legend(x = "topright",legend=c("YPR","SPR"), bty = "n", col=c(3,2), lwd=3)


lines(Fcr_est,SPR_est,type="p",lwd=10)
abline(h = 0.4, col = "black",lty = 2)
abline(v = Ftar, col = "black",lty = 2)




# -------------FIGURA CONTORNOS DE RIESGO LBPA -------------------------------
par(mfrow = c(1, 1))


files=read.admbFit('LBPA')
attach(files)
library(MASS)
library(MVN)
library(matrixcalc)

U = matrix(est[c(totPar-1, totPar)],nrow = 2,  ncol = 1)
V = cov[c(totPar-1, totPar),c(totPar-1, totPar)]
sd=std[c(totPar-1, totPar)]


p_low=1-pnorm(est[c(totPar-1)],0.4,std[c(totPar-1)])
p_high=pnorm(est[c(totPar)],est[totPar],std[c(totPar)])


set.seed(1)

#Generamos la población de tamaño N=1000
T = mvrnorm(1000,U, V)
x=T[,1]
y=T[,2]
z <- kde2d(x, y, n = 100)

plot(x,y, col='gray',pch = 19,xlab="SPR",ylab="F",main="Uncertainty levels")
     #xlim = c(0,max(x)),ylim=c(0,max(y)))
contour(z, lwd = 1, add = TRUE, nlevels=5, col = "blue")
lines(SPR_est,Fcr_est, type="p", pch = 19, cex=2)

abline(v = 0.4,  col = "black",lty = 2)
abline(h = est[totPar], col = "black",lty = 2)


# Curva probabilidad normal

eje=seq(0,1,by=0.01)
d1=dnorm(eje,U[1],sd[1])

plot(eje,d1/max(d1),type="l", ylab="Density", xlab="SPR",main=paste("p(SPR<0.4)=",round(p_low,3))) 
polygon(c(0, eje, SPR_est), c(0, d1/max(d1), 0), col = rgb(1, 0.27, 0, alpha = 0.2))

abline(v = 0.4,  col = "black",lty = 2)


eje=seq(0,2*U[2],by=0.01)
d2=dnorm(eje,U[2],sd[2])
plot(eje,d2/max(d2),type="l", ylab="Density", xlab="F",main=paste("p(F>Ftar)=",round(p_high,3)))
polygon(c(0, eje, U[2]), c(0, d2/max(d2), 0), col = c("cyan"))
abline(v = est[totPar],  col = "black",lty = 2)

#-----------------------------------------------
#---------------------------------------------------------------------------------------------------------
Nage<-data$Age_Length_s.e_N_Catch_Selectivity_F_Weight_Maturity[,4]
Cage<-data$Age_Length_s.e_N_Catch_Selectivity_F_Weight_Maturity[,5]
Fage<-data$Age_Length_s.e_N_Catch_Selectivity_F_Weight_Maturity[,7]
Nage0<-data$Age_Length_s.e_N_Catch_Selectivity_F_Weight_Maturity[,4]

M=-log(Nage[2]/Nage[1])-Fage[1]

# CALCULO N0
for (i in 2:nages) {
  Nage0[i]=Nage0[i-1]*exp(-M)
}
Nage0[nages]=Nage0[nages]/(1-exp(-M)) 


Cagelength<-ProbLen
Nagelength<-ProbLen
Nage0length<-ProbLen


for (i in 1:nages) {
  Cagelength[i,]<-Cage[i]*ProbLen[i,]/sum(Cage)
  Nagelength[i,]<-Nage[i]*ProbLen[i,]/sum(Cage)
  Nage0length[i,]<-Nage0[i]*ProbLen[i,]/sum(Cage)
}

par(mfrow = c(2,2))


plot(BinLen,colSums(Nage0length), type="l", lwd=3, col="black",
     ylab="Frequency", 
     xlab="Length",
     main="Unfished Population at-length",cex.main = 1.5)


lines(BinLen,Nage0length[1,], 
      type="l",cex.lab = 1.5, 
      lwd=4, 
      col="lightblue", 
      xlim = c(min(BinLen),1.1*max(BinLen)),
      ylim = c(0,max(Nage0length)))


for (i in 2:nages) {
  lines(BinLen,Nage0length[i,], type="l", lwd=4, col="lightblue")
}

#-----------------------------------------------------------
plot(BinLen,colSums(Nagelength), type="l", lwd=3, col="black",
     ylab="Frequency", 
     xlab="Length",
     main="Current population at-length",cex.main = 1.5)


lines(BinLen,Nagelength[1,], 
      type="l", cex.lab = 1.5,
      lwd=4, 
      col="lightblue", 
      xlim = c(min(BinLen),1.1*max(BinLen)),
      ylim = c(0,max(Nagelength)))


for (i in 2:nages) {
  lines(BinLen,Nagelength[i,], type="l", lwd=4, col="lightblue")
}

#--------------------------------------------------------
plot(BinLen,PredFre, type="l", lwd=3, col="black",
     ylab="Frequency", 
     xlab="Length",
     main="Predicted Catch at-length",cex.main = 1.5)


lines(BinLen,Cagelength[1,], 
      type="l", cex.lab = 1.5,
      lwd=4, 
      col="lightblue", 
      xlim = c(min(BinLen),1.1*max(BinLen)),
      ylim = c(0,max(Cagelength)))


for (i in 2:nages) {
  lines(BinLen,Cagelength[i,], type="l", lwd=4, col="lightblue")
}

