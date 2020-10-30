DATA_SECTION

  init_vector biolpar(1,8)
  init_number nages  
  init_number nrep
  init_vector pond(1,nrep)
  init_number nlength
  init_vector len_bins(1,nlength) 
  init_matrix LF_data(1,nrep,1,nlength)

  init_number L50prior // 
  init_number slopeprior // 
  init_number Fcrprior // 
  init_number Loprior //
  init_number s1prior //
  init_number s2prior // 


 // Coef of variation prior
  init_number cv4 // L50
  init_number cv99 // slope
  init_number cv100 // F
  init_number cv1 // Lr
  init_number cv2 // ao
  init_number cv3 // cv

 // Phases
  init_int  f3 // L50
  init_int  f4 // slope
  init_int  f2 // F
  init_int  f7 // Lo
  init_int  f5 // a0
  init_int  f6 // cv

  number logL50ini
  number logslopeini
  number logFcrini
  number logLoini
  number logs1ini
  number logs2ini

 !! logL50ini=log(L50prior);
 !! logslopeini=log(slopeprior);
 !! logFcrini=log(Fcrprior);
 !! logLoini=log(Loprior);
 !! logs1ini=log(s1prior+1E-5);
 !! logs2ini=log(s2prior);


  init_number ratio
  init_number h


  init_number nm // sample size

INITIALIZATION_SECTION

  log_Fcr    logFcrini
  log_L50    logL50ini
  log_rango  logslopeini
  log_alfa   logs1ini
  log_beta   logs2ini
  log_Lo     logLoini
  log_Ftar   -0.5

PARAMETER_SECTION

 init_number log_Fcr(f2) 
 init_number log_L50(f3) 
 init_number log_rango(f4) 
 init_number log_alfa(f5)
 init_number log_beta(f6)
 init_number log_Lo(f7)
 init_number log_Ftar(5)


 vector N0(1,nages)
 vector Ntar(1,nages)
 vector N(1,nages)
 vector Sel_a(1,nages)
 vector Sel(1,nlength)
 vector F(1,nages)
 vector Z(1,nages)
 vector S(1,nages)
 vector mu_edad(1,nages)
 vector sigma_edad(1,nages)
 vector wmed(1,nlength)
 vector msex(1,nlength)
 vector Ps(1,nages)
 vector pred_Ctot_a(1,nages)
 vector pred_Ctot(1,nlength)
 vector likeval(1,10)
 vector edades(1,nages)

 matrix prop_obs(1,nrep,1,nlength)
 vector prop_pred(1,nlength)

 number Linf
 number k
 number Lo
 number M
 number s
 number SPR
 number SPRtar
 number Fref
 number YPR
 number BPR
 
 
 matrix Prob_talla(1,nages,1,nlength) 
 matrix FrecL(1,nages,1,nlength)

 number alfa
 number beta
 number dts
 number B0
 number slope
 

 objective_function_value f

PRELIMINARY_CALCS_SECTION


 Linf=biolpar(1);
 k=biolpar(2);
 M=biolpar(3);
 dts=biolpar(8);
 Ps=0.0;
 

PROCEDURE_SECTION


 Prob_length2age();
 Pop_Dynamic();
 Log_likelihood();
  

FUNCTION Prob_length2age

  int i, j;

 edades.fill_seqadd(1,1);


// genero una clave edad-talla para otros calculos. Se modela desde L(1)
 mu_edad(1)=exp(log_Lo);
 for (i=2;i<=nages;i++)
  {
  mu_edad(i)=Linf*(1-exp(-k))+exp(-k)*mu_edad(i-1);
  }

  sigma_edad=exp(log_alfa)+exp(log_beta)*mu_edad;
  
  Prob_talla = ALK( mu_edad, sigma_edad, len_bins);

  slope=biolpar(7)-biolpar(6);

  wmed=exp(biolpar(4))*pow(len_bins,biolpar(5));
  msex=1./(1+exp(-log(19)*(len_bins-biolpar(6))/(slope)));

//----------------------------------------------------------------------
FUNCTION dvar_matrix ALK(dvar_vector& mu, dvar_vector& sig, dvector& x)
	//RETURN_ARRAYS_INCREMENT();
	int i, j;
	dvariable z1;
	dvariable z2;
	int si,ni; si=mu.indexmin(); ni=mu.indexmax();
	int sj,nj; sj=x.indexmin(); nj=x.indexmax();
	dvar_matrix pdf(si,ni,sj,nj);
	pdf.initialize();
	double xs=0.5*(x[sj+1]-x[sj]);
	for(i=si;i<=ni;i++) //loop over ages
	{
		 for(j=sj;j<=nj;j++) //loop over length bins
		{
			z1=((x(j)-xs)-mu(i))/sig(i);
			z2=((x(j)+xs)-mu(i))/sig(i);
			pdf(i,j)=cumd_norm(z2)-cumd_norm(z1);
		}//end nbins
		pdf(i)/=sum(pdf(i));
	}//end nage
	//RETURN_ARRAYS_DECREMENT();
	return(pdf);
//----------------------------------------------------------------------

FUNCTION Pop_Dynamic

// Selectivity at-length /////////
  Sel=1./(1+exp(-log(19)*(len_bins-exp(log_L50))/(exp(log_rango))));
// Selectivity at-age
  Sel_a=Prob_talla*Sel;


// Survival rate at-age /////////////////////
  F=exp(log_Fcr)*Sel_a;
  Z=F+M;
  S=exp(-1.*Z);


 // Unfished biomass////////////////////////////
  N0(1)=1.0;
  for (int j=2;j<=nages;j++)
  { N0(j)=N0(j-1)*exp(-1.*M);}
    N0(nages)=N0(nages)/(1-exp(-1.*M));
  B0=sum(elem_prod((N0*exp(-dts*M))*Prob_talla,elem_prod(wmed,msex)));

  alfa=4*h/(5*h-1);
  beta=(1-h)/(5*h-1)*B0;


 // LF estimation ////////////////////////////
  N(1)=1.0;
  for (int i=2;i<=nages;i++){
  N(i)=N(i-1)*exp(-Z(i-1));
  }
  N(nages)=N(nages)/(1-exp(-Z(nages)));
  pred_Ctot_a=elem_prod(elem_div(F,Z),elem_prod(1.-S,N));
  pred_Ctot=pred_Ctot_a*Prob_talla;

 // Proportions ////////////////////////////
  for (int j=1;j<=nrep;j++){
  prop_obs(j)=LF_data(j)/sum(LF_data(j));
  }
  prop_pred=pred_Ctot/sum(pred_Ctot);


 // SPR estimation and Ftarg////////////////////////////
  SPR=1/B0*(alfa*sum(elem_prod(elem_prod(N,exp(-dts*Z))*Prob_talla,elem_prod(wmed,msex)))-beta);

  Ntar(1)=1.0;
  Z=M+exp(log_Ftar)*Sel_a;
  
  for (int i=2;i<=nages;i++){
  Ntar(i)=Ntar(i-1)*exp(-Z(i-1));
  }
  Ntar(nages)=Ntar(nages)/(1-exp(-Z(nages)));

  // Per-recruit biomass and yield /////////////////////////
  SPRtar=1/B0*(alfa*sum(elem_prod(elem_prod(Ntar,exp(-dts*Z))*Prob_talla,elem_prod(wmed,msex)))-beta);


FUNCTION Log_likelihood


  s=0;
  for (int j=1;j<=nrep;j++){
  s+=-nm*sum(pond(j)*elem_prod(prop_obs(j),log(prop_pred)));// LF_data
  }


  likeval(1)=s;// LF_data
  likeval(2)=0.5*square((log_Lo-logLoini)/cv1);
  likeval(3)=0.5*square((log_alfa-logs1ini)/cv2);
  likeval(4)=0.5*square((log_beta-logs2ini)/cv3);
  likeval(5)=0.5*square((log_L50-logL50ini)/cv4);
  likeval(6)=0.5*square((log_rango-logslopeini)/cv99);
  likeval(7)=0.5*square((log_Fcr-logFcrini)/cv100);

  if(last_phase()){

  likeval(8)=1000*square((log(SPRtar)-log(ratio)));}
  

  f=sum(likeval);

 


REPORT_SECTION

  for (int i=1;i<=nages;i++){
    FrecL(i)=Prob_talla(i)*pred_Ctot_a(i);
  }


  report << "Length Based Pseudo-Cohort Analysis (LBPA) " << endl;
  report << "Canales C.M, A.E. Punt, & M. Mardones (2021) " << endl;
  report << "cristian.canales.r@pucv.cl " << endl;
  report << "-----------------------------------------------------------------------------" << endl;
  report << "Length bins" << endl;
  report << "---------------------------------------------------------------------------- " << endl;
  report <<  len_bins << endl;
  report << "---------------------------------------------------------------------------- " << endl;
  report << "Observed frequencies" << endl;
  report <<  prop_obs << endl;
  report << "Predicted frequency" << endl;
  report <<  prop_pred << endl;
  report << " " << endl;
  report << "Catch length frequency (columns) by age groups (rows)" << endl;
 report << "-------------------------------------------------------- " <<endl;
  report <<  FrecL  << endl;
  report << " " << endl;
  report << "Probability of length (columns) at-age (rows)" << endl;
 report << "-------------------------------------------------------- " <<endl;
  report <<  Prob_talla  << endl;
  report << " " << endl;
  report << " " << endl;


  report << "Age  Length  s.e   N   Catch   Selectivity   F    Weight  Maturity " << endl;
  report << "----------------------------------------------------- " << endl;
  for (int j=1;j<=nages;j++){ // l
  report << edades(j) <<" "<<mu_edad(j)<<" "<<sigma_edad(j)<<" "<<N(j)<<" "<<pred_Ctot_a(j)<<" "<<Sel_a(j)<<" "<<Sel_a(j)*exp(log_Fcr)<<" "<<Prob_talla(j)*wmed<<" "<<Prob_talla(j)*msex<<endl; 
  }


 report << " "<<endl;
 report << " "<<endl;

    for (int i=1;i<=nages;i++){
    if(mu_edad(i)>=biolpar(1)*2/3){
    Ps(i)=1;
    }};
  

 report << "Length frequency of exploitable population : current, target, unfished" << endl;
 report << "-------------------------------------------------------- " <<endl;
 report << elem_prod(N,Sel_a)*Prob_talla<<endl;
 report << elem_prod(Ntar,Sel_a)*Prob_talla<<endl;
 report << elem_prod(N0,Sel_a)*Prob_talla<<endl;

 report << " "<<endl;
 report << " "<<endl;
 
 report << "Model parameters " << endl;
 report << "----------------------------------------------------- " <<endl;
 report<<"F	 L50	 slope   a0    cv  Lr      Ftar"<<endl;
 report<<exp(log_Fcr)<<" "<<exp(log_L50)<<" "<<exp(log_rango)<<" "<<exp(log_alfa)<<" "<<exp(log_beta)<<" "<<exp(log_Lo)<<" "<<exp(log_Ftar)<<endl;
 report << " "<<endl;
 report << " "<<endl;
 report<<"F/Ftar  SPR"<<endl;
 report<<exp(log_Fcr)/exp(log_Ftar)<<" "<<SPR<<endl;

 
  report << " " << endl;
  report << "Log-likelihood components" << endl;
 report << "----------------------------------------------------- " <<endl;
  report << "Proportions" << endl;
  report << likeval(1) << endl;
  report << "Lr" << endl;
  report << likeval(2) << endl;
  report << "a0" << endl;
  report << likeval(3) << endl;
  report << "cv" << endl;
  report << likeval(4) << endl;
  report << "L50" << endl;
  report << likeval(5) << endl;
  report << "slope" << endl;
  report << likeval(6) << endl;
  report << "F" << endl;
  report << likeval(7) << endl;
  report << "SPRtar" << endl;
  report << likeval(8) << endl;
  report << "Total" << endl;
  report << sum(likeval) << endl;



FINAL_SECTION

 

 ofstream print_pr("per_recruit.txt");

 //Per-recruit Analysis

  print_pr <<"F     Y/R     SSB/R"<<endl;
 //print_pr <<"-------------------"<<endl;

   Fref=0.0;
   while(Fref<=5*M){
   
    F=Fref*Sel_a;
    Z=F+M;
    S=exp(-1.*Z);
 
    N(1)=1.0;
      for (int i=2;i<=nages;i++){
    N(i)=N(i-1)*exp(-Z(i-1));
    }

  N(nages)=N(nages)/(1-exp(-Z(nages)));



  BPR=alfa*sum(elem_prod(elem_prod(N,exp(-dts*Z))*Prob_talla,elem_prod(wmed,msex)))-beta;
  YPR=(alfa*BPR/(beta+BPR))*sum(elem_prod(elem_prod(elem_div(F,Z),elem_prod(1.-S,N))*Prob_talla,wmed));

  print_pr <<Fref<<" "<<YPR<<" "<<BPR<<endl;

  Fref+=0.05;
  }

 

