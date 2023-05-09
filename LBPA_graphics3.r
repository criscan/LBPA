rm(list=ls(all=TRUE)) # erasure all objects
system('./erase LBPA.rep')  
system('./erase LBPA.std')  


library(knitr)

system('./LBPA -ind lbpa.dat')  # for model running


source('read.admb.R')
source('read.admbFit.R')

data <-read.rep('For_R.rep')

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



#Model fit and residuals-------------------------------------

par(mfrow = c(2, 2))


plot(BinLen, ObsFre[1,], col="gray",
     pch=20,
     ylab="Frequency",
     xlab="Length",
     ylim=c(0, max(ObsFre)),
     xlim = c(min(BinLen),1.05*max(BinLen)),
     cex=1.5,
     cex.lab = 1.,
     cex.axis = 1.)

for (i in 1:NObsFre) {
  lines(BinLen, ObsFre[i,], type="p",col="gray",pch=20,cex=1.5)
}


lines(BinLen,
      PredFre,type="l",col="red",
      lwd=2)


unos=matrix(1,length(ObsFre[,1]),1)
res=(ObsFre-unos%*%PredFre)
res=res/sd(res)


hist(res,xlab="Normalized residuals",prob=T,main="",cex.lab = 1.,
     cex.axis = 1.)
x <- seq(min(res), max(res), length = 200)
lines(x, dnorm(x), col = "red", lwd = 2)
box()



#Theorical length comps-------------------------------------


# 1

L1 <- (data$Length_frequency_of_exploitable_population_current_target_unfished[1,])
L2 <- (data$Length_frequency_of_exploitable_population_current_target_unfished[2,])
L3 <- (data$Length_frequency_of_exploitable_population_current_target_unfished[3,])


#Selectivity and maturity-------------------------------------
S1 <- (data$Selectivity_and_maturity_at_length[1,])
S2 <- (data$Selectivity_and_maturity_at_length[2,])

plot(BinLen,S1, ylab="Proportion", xlab="Length", type="l", lwd=3,  xlim = c(min(BinLen),1.05*max(BinLen)),cex=1.5,cex.lab = 1.0,
     cex.axis = 1.)
lines(BinLen,S2, lwd=3, col="green", xlim = c(min(BinLen),1.1*max(BinLen)))
grid(nx = NULL, ny = NULL, lty = 2, col = "gray",lwd = 1)    

legend(x = "topleft",
       legend=c("Selectivity","Maturity"),
       bty = "n",
       col=c("black","green"),
       lwd=3)

#Per recruit------------------------------------
Fcr_est=results[1]
Ftar=results[7]
SPR_est=puntos[2]
SPR_tar=puntos[3]



plot(Fcr,YPR/max(YPR),type="l", ylab="YPR, SPR", xlab="Fishing mortality", lwd=3, col="green",
     cex=1.5,cex.lab = 1.,cex.axis = 1.)
lines(Fcr,SPR/max(SPR),type="l", ylab="YPR, SPR", lwd=3, col="red")

lines(Fcr_est,SPR_est,type="p",lwd=10)

abline(h = SPR_tar, col = "black",lty = 2)
abline(v = Ftar, col = "black",lty = 2)
legend(x = "topright",
       legend=c("Yield","Biomass"),
       bty = "n",
       col=c("black","green"),
       lwd=3)

par(mfrow = c(2, 2))

plot(c(BinLen,tail(BinLen,1)),c(L3,0), ylab="Frequency", xlab="Length", type="l", 
     ylim = c(0,max(c(L1,L2,L3))),     
     cex=1.5,
     cex.lab = 1.,
     cex.axis = 1.)
polygon(c(BinLen,tail(BinLen,1),BinLen[1]),c(L3,0,L3[1]),col = rgb(0.6, 0.98, 0.6, alpha = 0.3)) 
polygon(c(BinLen,BinLen[1]),c(L2,L2[1]),col = rgb(0.12, 0.56, 1, alpha = 0.3)) 
polygon(c(BinLen,BinLen[1]),c(L1,L1[1]),col = rgb(0.8, 0.2, 0.2, alpha = 0.5))

legend(x = "topleft",
       legend=c("Unfished","Target","Current"),
       bty = "n",
       col=c(rgb(0.6, 0.98, 0.6),rgb(0.12, 0.56, 1),rgb(0.8, 0.2, 0.2)),
       lwd=3)


#Risk and uncertainty------------------------------------

files=read.admbFit('LBPA')
attach(files)
library(MASS)
library(MVN)
library(matrixcalc)

U = matrix(est[c(totPar-1, totPar)],nrow = 2,  ncol = 1)
V = cov[c(totPar-1, totPar),c(totPar-1, totPar)]
sd=std[c(totPar-1, totPar)]


p_low=1-pnorm(SPR_est,SPR_tar,std[c(totPar-1)])
p_high=pnorm(Fcr_est,Ftar,std[c(totPar)])
lispr=SPR_est-1.96*std[c(totPar-1)]
lsspr=SPR_est+1.96*std[c(totPar-1)]
liF=(Fcr_est-1.96*std[c(totPar)])/Ftar
lsF=(Fcr_est+1.96*std[c(totPar)])/Ftar



#Remuestreo N=1000
#set.seed(1)

T = mvrnorm(10000,U, V)
x=T[,1] # SPR
y=T[,2] # F
z <- kde2d(x, y, n = 100)

plot(x[1:1000],y[1:1000], col='gray',pch = 20,xlab="SPR",ylab="Fishing Mortality",
     cex=1.5,cex.lab = 1.0,cex.axis = 1.0)
contour(z, lwd = 1, add = TRUE, nlevels=5, col = "black")
lines(SPR_est,Fcr_est, type="p", pch = 19, cex=1.5)

abline(v = SPR_tar,  col = "black",lty = 2)
abline(h = Ftar, col = "black",lty = 2)



#par(mfrow = c(2, 1))


dspr=density(x)
plot(dspr,ylab="Density", xlab="SPR",main=paste("p(SPR<SPRtar)=",round(p_low,3),"; CI_SPR= [",round(lispr,2),"-",
                                                round(lsspr,2),"]"),cex.main=1.0) 
polygon(dspr, col = rgb(1, 0.27, 0, alpha = 0.2))
abline(v = SPR_tar,  col = "black",lty = 2)


df=density(y)
plot(df, ylab="Density", xlab="F",main=paste("p(F>Ftar)=",round(p_high,3),"; CI_F= [",round(liF,2),"-",round(lsF,2),"]"),cex.main=1.0) 
polygon(df, col = c("cyan"))
abline(v=Ftar,  col = "black",lty = 2)

#Number at-age----------------------------------------


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
     main="Unfished Population at-length",cex.main = 1.)
grid(nx = NULL, ny = NULL, lty = 2, col = "gray",lwd = 1)   

lines(BinLen,Nage0length[1,], 
      type="l",cex.lab = 1.5, 
      lwd=4, 
      col="lightblue", 
      xlim = c(min(BinLen),1.1*max(BinLen)),
      ylim = c(0,max(Nage0length)))
grid(nx = NULL, ny = NULL, lty = 2, col = "gray",lwd = 1)   

for (i in 2:nages) {
  lines(BinLen,Nage0length[i,], type="l", lwd=4, col="lightblue")
}




plot(BinLen,colSums(Nagelength), type="l", lwd=3, col="black",
     ylab="Frequency", 
     xlab="Length",
     main="Current population at-length",cex.main = 1.)


lines(BinLen,Nagelength[1,], 
      type="l", cex.lab = 1.5,
      lwd=4, 
      col="lightblue", 
      xlim = c(min(BinLen),1.1*max(BinLen)),
      ylim = c(0,max(Nagelength)))
grid(nx = NULL, ny = NULL, lty = 2, col = "gray",lwd = 1)   

for (i in 2:nages) {
  lines(BinLen,Nagelength[i,], type="l", lwd=4, col="lightblue")
}

plot(BinLen,PredFre, type="l", lwd=3, col="black",
     ylab="Frequency", 
     xlab="Length",
     main="Predicted Catch at-length",cex.main = 1.)


lines(BinLen,Cagelength[1,], 
      type="l", cex.lab = 1.5,
      lwd=4, 
      col="lightblue", 
      xlim = c(min(BinLen),1.1*max(BinLen)),
      ylim = c(0,max(Cagelength)))


for (i in 2:nages) {
  lines(BinLen,Cagelength[i,], type="l", lwd=4, col="lightblue")
}
grid(nx = NULL, ny = NULL, lty = 2, col = "gray",lwd = 1)  

edad=c(1:nages)
barplot(Cage~edad,col="lightblue",xlab="Age",ylab="Frecuency",
     main="Predicted Catch at-age",cex.main = 1.)
grid(nx = NULL, ny = NULL, lty = 2, col = "gray",lwd = 1)
box()


tabla1 <- matrix(ncol=1, round(data$F_L50_slope_a0_cv_Lr_Ftar[1:6], 2))
rownames(tabla1) <- c("Fishing Mortality (F)", "Selectivity length at 50% (L50)", "Selectivity slope (d)",
                     "Invariant std deviation in length (a0)", "Coeff of variation length at-age (cv)",
                     "Size of recruits (Lr)")


YPRcur=data$YPRcurr_YPRtar[1]
YPRtar=data$YPRcurr_YPRtar[2]

B0=SPR[1]
tabla2 <- matrix(ncol=1, round(c(B0,B0*SPR_est, B0*SPR_tar,SPR_est,SPR_tar,Ftar,Fcr_est/Ftar,YPRcur,YPRtar),2))
rownames(tabla2) <- c("Virginal biomass per-recruit (BPR0)", "Current BPR", "Target BPR","Current spawning potential ratio (SPR)",
                      "Target SPR (SPRtar)", "Target fishing mortality (Ftar)","Overfishing index (F/Ftar)",
                      "Current yield per-recruit (YPRcur)","Target  yield per-recruit (YPRtar)")


knitr::kable(tabla1,"simple",caption = "1: Estimated model parameters")
knitr::kable(tabla2,"simple",caption = "2: Derivates quantities")

