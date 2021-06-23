rm(list=ls(all=TRUE)) # erasure all objects
setwd('G:/Mi unidad/Instituto_Nacional_de_Pesca/Camaron_Pomada/LBPA') # set current directory

#system('./LBPA -ind lbpa4.dat')  # for model running
#system('./LBPA -ind lbpa8.dat')  # for model running
system('./LBPA -ind lbpa12.dat')  # for model running

source('C:/Users/cristian/Desktop/LBPA-master/LBPA-master/read.admb.R')
data <-read.rep('For_R.rep')

nages<-6 #this depends on the specie

nages <- length(data$Probability_of_length[,1])
age <- seq(1,nages) 
BinLen <- data$Length_bins
NObsFre <- length(data$Observed_frequencies[,1]) #this numbers of observations depend of the own data
ObsFre <- data$Observed_frequencies
PredFre <- data$Predicted_frequency
CatchLFre<- data$Catch_length_frequency
ProbLen <- data$Probability_of_length


par(mfrow = c(2, 2))

plot(BinLen,
     PredFre,type="l",
     ylab="Frequency",
     xlab="Length",
     lwd=3,
     ylim=c(0, max(ObsFre)))


for (i in 1:NObsFre) {
  lines(BinLen, ObsFre[i,], type="p", col=2)
}


L1 <- (data$Length_frequency_of_exploitable_population_current_target_unfished[1,])
L2 <- (data$Length_frequency_of_exploitable_population_current_target_unfished[2,])
L3 <- (data$Length_frequency_of_exploitable_population_current_target_unfished[3,])


plot(BinLen,L3, ylab="Frequency", xlab="Length", type="h", lwd=5, col="green", xlim = c(min(BinLen),1.1*max(BinLen)))
lines(BinLen,L2,type="l", lwd=3, col="red")
lines(BinLen,L1,type="l", lwd=3)
legend(x = "topright",
       legend=c("Current","Target","Unfished"),
       bty = "n",
       col=c(1,2,3),
       lwd=3)

Fcr <- data$`F_Y/R_SSB/R`[,1]
YPR <- data$`F_Y/R_SSB/R`[,2]
SPR <- data$`F_Y/R_SSB/R`[,3]

plot(Fcr,YPR/max(YPR),type="l", ylab="YPR, SPR", xlab="Fishing Mortality", lwd=3, col="green")
lines(Fcr,SPR/max(SPR),type="l", ylab="YPR, SPR", xlab="Fishing Mortality", lwd=3, col="red")
legend(x = "topright",legend=c("YPR","SPR"), bty = "n", col=c(3,2), lwd=3)

puntos=data$`F/Ftar_SPR_SPRtar`
results=data$F_L50_slope_a0_cv_Lr_Ftar

Fcr_est=results[1]
Ftar=results[7]
SPR_est=puntos[2]

lines(Fcr_est,SPR_est,type="p",lwd=10)
abline(h = 0.4, col = "black",lty = 2)
abline(v = Ftar, col = "black",lty = 2)


S1 <- (data$Selectivity_and_maturity_at_length[1,])
S2 <- (data$Selectivity_and_maturity_at_length[2,])

plot(BinLen,S1, ylab="Frequency", xlab="Length", type="l", lwd=3,  xlim = c(min(BinLen),1.1*max(BinLen)))
lines(BinLen,S2, ylab="Frequency", xlab="Length", lwd=3, col="green", xlim = c(min(BinLen),1.1*max(BinLen)))
legend(x = "topright",legend=c("Selectivity","Maturity"), bty = "n", col=c(1,3), lwd=3)

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

plot(BinLen,Nage0length[1,], 
     ylab="Frequency", 
     xlab="Length", 
     type="l",cex.lab = 1.5, 
     lwd=4, 
     col="lightblue", 
     xlim = c(min(BinLen),1.1*max(BinLen)),
     ylim = c(0,max(Nage0length)),
     main="Unfished Population at-length",cex.main = 1.5)


for (i in 2:nages) {
  lines(BinLen,Nage0length[i,], type="l", lwd=4, col="lightblue")
}

lines(BinLen,colSums(Nage0length), type="l", lwd=3, col="black")
#-----------------------------------------------------------

plot(BinLen,Nagelength[1,], 
     ylab="Frequency", 
     xlab="Length", 
     type="l", cex.lab = 1.5,
     lwd=4, 
     col="lightblue", 
     xlim = c(min(BinLen),1.1*max(BinLen)),
     ylim = c(0,max(Nagelength)),
     main="Current population at-length",cex.main = 1.5)


for (i in 2:nages) {
  lines(BinLen,Nagelength[i,], type="l", lwd=4, col="lightblue")
}

lines(BinLen,colSums(Nagelength), type="l", lwd=3, col="black")
#--------------------------------------------------------

plot(BinLen,Cagelength[1,], 
     ylab="Frequency", 
     xlab="Length", 
     type="l", cex.lab = 1.5,
     lwd=4, 
     col="lightblue", 
     xlim = c(min(BinLen),1.1*max(BinLen)),
     ylim = c(0,max(Cagelength)),
     main="Predicted Catch at-length",cex.main = 1.5)


for (i in 2:nages) {
  lines(BinLen,Cagelength[i,], type="l", lwd=4, col="lightblue")
}
lines(BinLen,PredFre, type="l", lwd=3, col="black")


#--------------------------------------------------------



