rm(list=ls(all=TRUE)) # erasure all objects
setwd('C:/..............') # set current directory
system('./LBPA')  # for model running
source('C:/................/read.admb.R')
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

tabla <- matrix(ncol=1, round(c(data$F_L50_slope_a0_cv_Lr_Ftar,Fcr_est/Ftar,SPR_est),2))
rownames(tabla) <- c("F curr", "L50", "Slope", "a0", "cv", "Lr", "F Target","F curr/F Target","SPR")
colnames(tabla) <- c("Value")
outcomes<-as.table(tabla)
outcomes

