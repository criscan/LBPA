clear


B=read("lbpa.rep",-1,1,'(A)');
format('v',5)

ubi0=find(B=="Length bins");
Talla=evstr(B(ubi0+2));

ubi1=find(B=="Observed frequencies");
ubi2=find(B=="Predicted frequency");
ubi3=find(B=="Catch length frequency (columns) by age groups (rows)");
ubi4=find(B=="Probability of length (columns) at-age (rows)");
ubi66=find(B=="Age  Length  s.e   N   Catch   Selectivity   F    Weight  Maturity")
ubi5=find(B=="Length frequency of exploitable population : current, target, unfished");
ubi6=find(B=="Selectivity and maturity at-length");
ubi7=find(B=="Model parameters");


Fobs=evstr(B(ubi1+1:ubi2-1));
Fpred=evstr(B(ubi2+1));
frec_age=evstr(B(ubi3+2:ubi4-2));
Prob=evstr(B(ubi4+1:ubi66-3));


Nact=evstr(B(ubi5+2));
Ntar=evstr(B(ubi5+3));
N0=evstr(B(ubi5+4));
Sel=evstr(B(ubi6+2));
Mat=evstr(B(ubi6+3));

figure(1,"BackgroundColor",[1,1,1])

subplot(2,2,1),plot(Talla, Fobs,'r');plot(Talla,Fpred,'k','linewidth',3)
xlabel("Length",'FontSize',4);
ylabel("Frequency",'FontSize',4);

subplot(2,2,2),plot(Talla,N0,'g','linewidth',2),plot(Talla,Ntar,'r','linewidth',1);plot(Talla,Nact,'k','linewidth',2)
xlabel("Length",'FontSize',4);
ylabel("Frequency",'FontSize',4);
legend('Unfished','Target','Current')

subplot(2,2,3),plot(Talla, frec_age','c','linewidth',2);plot(Talla,sum(frec_age,'r'),'k')
xlabel("Length",'FontSize',4);
ylabel("Frequency",'FontSize',4);

subplot(2,2,4),plot(Talla, Sel,'k','linewidth',2);plot(Talla,Mat,':k')
xlabel("Length",'FontSize',4);
ylabel("Proportion",'FontSize',4);
legend('Selectivity','Maturity',4)


 ubi=find(B=="Model parameters")
 pars=evstr(B(ubi+3));
 Ftar=pars(7); 
 
 Fcr=pars(1);
 SPR=evstr(B(ubi+8))(2);
 SPRtar=evstr(B(ubi+8))(3);

 ubi=find(B=="F     Y/R     SSB/R")

 F=[];BPR=[];YPR=[];
 i=2;
 
 while  evstr(B(ubi+i))(1)<=0.9*3*Fcr then

 F=[F;evstr(B(ubi+i))(1)];
 YPR=[YPR; evstr(B(ubi+i))(2)];
 BPR=[BPR; evstr(B(ubi+i))(3)];     
 i=i+1;

 end
 
 figure(2,"BackgroundColor",[1,1,1])
 plot(F,BPR/BPR(1),'k',F,YPR/max(YPR),':r')
 plot(Fcr,SPR,'ok','MarkerSize',10,'MarkerFaceColor','r');
 plot([0 max(F)],[SPRtar SPRtar],':k')
 plot([Ftar Ftar],[0 1.0],':k')
 
 ylabel('SPR, YPR','FontSize',4);
 xlabel('Fishing Mortality','FontSize',4);
 a=gca();
 legend('SPR','YPR')
// a.data_bounds=[0,0;pars(3)*2,1];
 a.y_location = "origin";
 title(['F/Fmsy= ' string(Fcr/Ftar) "  SPR= " string(SPR)],'FontSize',4)



 
//
//
//
////
////
//
//
