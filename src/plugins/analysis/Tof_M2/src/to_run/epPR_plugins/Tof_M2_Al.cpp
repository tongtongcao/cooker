#include <Tof_M2.h>
//#include <Det_Ac.h>
#include <iostream>
#include <cmath>
#include "TStyle.h"
#include "TLine.h"
#include <string>
#include <sstream>
#include <fstream>
#include "TVector.h"
// Constructor/destructor
Tof_M2::Tof_M2(TTree *in, TTree *out,TFile *inf_, TFile * outf_, TObject *p):Plugin(in,out,inf_,outf_,p)
{

};

Tof_M2::~Tof_M2()
{

};
Long_t Tof_M2::histos(){
  gStyle->SetOptFit(1101);
  gStyle->SetOptStat(0);
  gStyle->SetNdivisions(510);
  gStyle->SetFitFormat("5.5g");
  gStyle->SetTitleFontSize(0.06);
  gStyle->SetPadTopMargin(0.11);
  gStyle->SetPadLeftMargin(0.14);
  gStyle->SetPadRightMargin(0.08);
  gStyle->SetPadBottomMargin(.14);
  gStyle->SetTitleYOffset(1.);
  gStyle->SetTitleXOffset(0.9);
  gStyle->SetTitleYSize(0.06);
  gStyle->SetTitleXSize(0.075);
  gStyle->SetLabelFont(62,"X");
  gStyle->SetLabelFont(62,"Y");
  gStyle->SetTitleFont(62,"X");
  gStyle->SetTitleFont(62,"Y");
  gStyle->SetLabelSize(0.045,"X");
  gStyle->SetLabelSize(0.045,"Y");

  int adc_c = 400;
  int Ep_c = 4;
  final_PgCAdc =  new TH1F( "final_PgCAdc","final_PgCAdc",70,30,330);
  final_PgCEp =  new TH1F( "final_PgCEp","final_PgCEp",70,0,4);
  final_beta =  new TH1F( "final_beta","final_beta",100,0.4,2);
  final_acAdc =  new TH1F( "final_acAdc","final_acAdc",30,30,4200);
  final_acAdcU =  new TH1F( "final_acAdcU","final_acAdcU",30,30,4200);
  final_acAdcD =  new TH1F( "final_acAdcD","final_acAdcD",30,30,4200);
  final_M2 =  new TH1F( "final_M2","final_M2",100,-10000,30000);
  vert_Mom =  new TH1F( "vert_Mom","vert_Mom",100,0,280);
  hBetaGap =  new TH2D( "hBetaGap","hBetaGap",12,0,12,100,0.4,2);
  theta_pzx = new TH2D("theta_pzx","theta_pzx",100,-0.02,0.02,100,200,300);
  phi_pzx = new TH2D("phi_pzx","phi_pzx",100,-0.004,0.004,100,200,300);
  p_pzx = new TH2D("p_pzx","p_pzx",100,0,280,100,200,300);
  h1adc_acuGap = dH2("Acu vs gap","Acu vs gap",12,-0.5,11.5,80,100,4000);
  h1adc_acGapuGap = dH2("AcuGap vs gap Tof2","AcuGap vs gapTof2",12,-0.5,11.5,12,-0.5,11.5);
  h1adc_acGapdGap = dH2("AcdGap vs gap Tof2","AcdGap vs gap Tof2",12,-0.5,11.5,12,-0.5,11.5);
  h1adc_acGapuGap2 = dH2("AcuGap = Tof2 Gap","AcuGap = Tof2 Gap",12,-0.5,11.5,12,-0.5,11.5);
  h1adc_acGapdGap2 = dH2("AcdGap = Tof2 Gap","AcdGap = Tof2 Gap",12,-0.5,11.5,12,-0.5,11.5);

  for(int i = 0; i < 12; i++){
    string name("ACu_ch");
    name+=(i+1)/10+'0';
    name+=(i+1)%10+'0';
    name+="_ac";
    string title("adc Ch for ac u");
    title+=(i+1)/10+'0';
    title+=(i+1)%10+'0';
    title+="Ch_ADC";
    h1adc_ACu[i] = dH1(name.c_str(),title.c_str(),20,0,4200);
    h1adc_ACuGap[i] = dH2(Form("AcU ADC vs Gap %d",i),Form("AcU ADC vs Gap %d",i),12,0,13,30,0,4000);

    name=string("ACd_ch");
    name+=(i+1)/10+'0';
    name+=(i+1)%10+'0';
    name+="_ac";
    title=string("adc Ch for ac d");
    title+=(i+1)/10+'0';
    title+=(i+1)%10+'0';
    title+="Ch_ADC";
    h1adc_ACd[i] = dH1(name.c_str(),title.c_str(),20,0,4200);

    h1tdc_ACd[i] = dH1(Form("ACD TDC Gap %d",i+1),Form("ACD TDC Gap %d",i+1),30,1200,1800);
    h1tdc_ACu[i] = dH1(Form("ACU TDC Gap %d",i+1),Form("ACU TDC Gap %d",i+1),30,1200,1800);

    h2adc_ACUD[i] = dH2(Form("ACD vs ACU %d",i+1),Form("ACD vs ACU %d",i+1),50,0,4000,50,0,4000);
    h1adc_ACdGap[i] = dH2(Form("AcD ADC vs Gap %d",i),Form("AcD ADC vs Gap %d",i),12,0,13,30,0,4000);

    string name1("PgC_ch");
    name1+=(i+1)/10+'0';
    name1+=(i+1)%10+'0';
    name1+="_PgC";
    string title1("ADC");
    title1+=(i+1)/10+'0';
    title1+=(i+1)%10+'0';

    h1adc_PgC[i] = dH1(name1.c_str(),title1.c_str(),20,5,adc_c);
    h1adc_PgC[i]->GetXaxis()->SetTitle("ADC Channel");
    h1adc_PgC[i]->GetXaxis()->CenterTitle();
    
    string name1kk("PgC_cha");
    name1kk+=(i+1)/10+'0';
    name1kk+=(i+1)%10+'0';
    name1kk+="_PgC";
    string title1kk("E/p");
    title1kk+=(i+1)/10+'0';
    title1kk+=(i+1)%10+'0';

    h1AcAv[i] = dH1(Form("AcAvGap_%d ",i+1),Form("AcAvGap_%d ",i+1),20,0,4000);
    h1AcAv[i]->GetXaxis()->SetTitle("Ac average (GM)");
    h1AcAv[i]->GetXaxis()->CenterTitle();



    h1EP_PgC[i] = dH1(name1kk.c_str(),title1kk.c_str(),20,0,Ep_c);
    h1EP_PgC[i]->GetXaxis()->SetTitle("E/p");
    h1EP_PgC[i]->GetXaxis()->CenterTitle();
    
    h1Mom_PgC[i] = dH1(Form("PgC_Mom_Gap_%d ",i+1),Form("PgC_Mom_Gap_%d ",i+1),60,40,200);
    h1Mom_PgC[i]->GetXaxis()->SetTitle("Track Momentum (MeV)");
    h1Mom_PgC[i]->GetXaxis()->CenterTitle();
    string n1("PgC Gap_");
    n1+=(i+1)/10+'0';
    n1+=(i+1)%10+'0';
    n1+="mom";
    string t1("PgC Gap_");
    t1+=(i+1)/10+'0';
    t1+=(i+1)%10+'0';
    t1+="Mom";
    h1Mom_ADC[i] = dH2(n1.c_str(),t1.c_str(),100,0,500,100,0,300);
    h1Mom_ADC[i]->GetYaxis()->SetTitle("Momentum (MeV)");
    h1Mom_ADC[i]->GetYaxis()->CenterTitle();
    h1Mom_ADC[i]->GetXaxis()->SetTitle("ADC channel");
    h1Mom_ADC[i]->GetXaxis()->CenterTitle();

    string tpn1("PgC Gap_tp_");
    tpn1+=(i+1)/10+'0';
    tpn1+=(i+1)%10+'0';
    tpn1+="momtp";
    string tpt1("PgC Gap_tp");
    tpt1+=(i+1)/10+'0';
    tpt1+=(i+1)%10+'0';
    tpt1+="Momtp";
    h1Mom_ADC_TP[i] = new TProfile(tpn1.c_str(),tpt1.c_str(),100,0,500,0,300);
    h1Mom_ADC_TP[i]->SetLineColor(kRed);
    h1Mom_ADC_TP[i]->SetMarkerColor(kRed);
    h1Mom_ADC_TP[i]->GetYaxis()->SetTitle("Momentum (MeV)");
    h1Mom_ADC_TP[i]->GetYaxis()->CenterTitle();
    h1Mom_ADC_TP[i]->GetXaxis()->SetTitle("ADC channel");
    h1Mom_ADC_TP[i]->GetXaxis()->CenterTitle();

    for(int ii=0;ii<7;ii++){
      string namep("Adc_PGC");
      namep+=(i)/10+'0';
      namep+=(i)%10+'0';
      namep+= "_";
      namep+=(ii)+'0';
      namep+="_pgc";

      string titlep("PgC Gap_");
      titlep+=(i)/10+'0';
      titlep+=(i)%10+'0';
      titlep+= "_";
      titlep+=(ii)+'0';
      titlep+="";

      ADC_PgCh[i][ii]= dH1(namep.c_str(),titlep.c_str(),12,0,adc_c);
      ADC_PgCh[i][ii]->GetXaxis()->SetTitle("ADC Channel");
      ADC_PgCh[i][ii]->GetXaxis()->CenterTitle();

      string namep11("E/p_PGC");
      namep11+=(i)/10+'0';
      namep11+=(i)%10+'0';
      namep11+= "_";
      namep11+=(ii)+'0';
      namep11+="_pgc";

      string titlep11("PgC Gap_");
      titlep11+=(i)/10+'0';
      titlep11+=(i)%10+'0';
      titlep11+= "_";
      titlep11+=(ii)+'0';
      titlep11+="";
      EP_PgCh[i][ii]= dH1(namep11.c_str(),titlep11.c_str(),12,0,Ep_c);
      EP_PgCh[i][ii]->GetXaxis()->SetTitle("ADC/p");
      EP_PgCh[i][ii]->GetXaxis()->CenterTitle();

    }
  }
}
//For tof2
double z_tofc[12] = {160.0, 161.3, 162.1, 162.0, 153.6, 153.9, 153.5, 161.1, 160.8, 160.0, 159.8, 160.1}; //z of each tof
double x_tofc[12] = {139.0,139.0,139.0,139.0,131.6,139.0,131.6,139.0,139.0,139.0,139.0,139.0}; //x of each tof
double x_tof[12], y_tof[12];
double t_tof[12], z_tofb[12];   
double pl_btof[12]; 
double x_btof[12], y_btof[12];  
double x_tof_max[12], x_tof_min[12];
Double_t xv_tof = 80.; Double_t yv_tof = 30.; Double_t zv_tof = 2.; //All Tof has same x and y size
Double_t MEE_tof = 68.7; Double_t Z_A_tof = 0.53768;
Double_t ploss_tof_b[12]; 
Double_t rho_tof = 1.06;

//For Deg

Double_t z_degc[12] = {171.5, 172.2, 172.8, 172.8, 172.9, 173.0, 172.6, 172.2, 171.8, 171.5, 171.2, 171.4}; //z of each tof
Double_t x_degc[12] = {139.0, 139.0, 139.0, 139.0, 131.6, 139.0, 131.6, 139.0, 139.0, 139.0, 139.0, 139.0}; //x of each tof
Double_t x_deg[12],y_deg[12],z_degf[12],z_degb[12],pl_fdeg[12],pl_bdeg[12];
Double_t x_bdeg[12],y_bdeg[12],x_fdeg[12],y_fdeg[12];
Double_t x_deg_max[12], x_deg_min[12];
Double_t xv_deg = 93.; Double_t yv_deg = 41.; Double_t zv_deg =10.; //All Tof has same x and y size
Double_t MEE_deg = 57.4; Double_t Z_A_deg = 0.57034;
Double_t ploss_deg_b[12]; 
Double_t ploss_deg_f[12]; 
Double_t rho_deg =0.89;
Double_t M;
Double_t acGap[12]={0};
Double_t E_PgC[12];
Double_t P_PgC[12];
//Double_t ac_dw[12] ={1155, 735, 735, 735, 945, 735, 945, 735,1155, 945, 735, 945};//875
//Double_t ac_up[12] ={735, 945, 735, 945,1155, 945, 945, 735, 945, 945, 945, 945};//840 

Double_t Tof_M2::pathlength(Double_t z0, Double_t z, Double_t n_z)
{
  Double_t pl = (z -z0)/n_z;
  return pl;
}

Double_t Tof_M2::x_p(Double_t x0, Double_t n_x, Double_t t_p)
{
  Double_t pos_x = x0 + n_x*t_p;
  return pos_x;
}

Double_t Tof_M2::y_p(Double_t y0, Double_t n_y, Double_t t_p)
{
  Double_t pos_y = y0 + n_y*t_p;
  return pos_y;
}

Double_t Tof_M2::dp(Double_t p0, Double_t t_p, Double_t Z_A, Double_t I, Double_t rho)
{
  Double_t E = sqrt(p0*p0 + M*M);
  Double_t K = 0.307075;
  Double_t dE = - rho*t_p*(K/2*Z_A)*(2.*log(2.*M*pow(10,6)/I) + 4.*log(E/M) - 2.);
  Double_t dp_m = E*dE*pow((E*E - M*M),(-0.5));

  return dp_m;
}

Double_t Tof_M2::E_loss(Double_t p0, Double_t p_loss){
  Double_t dE_loss = p0/sqrt(p0*p0+M*M) * p_loss;
}

Double_t E_loss_tof2[12], E_loss_degf[12], E_loss_degb[12];
static Int_t run,event;
static vector<Int_t> gapToF1, gapToF2;
static vector<double> M2, vert_p, Ep, beta;


Long_t Tof_M2::startup()
{
  getBranchObject("testBranch",(TObject **) &ptr);
  getBranchObject("tofBranch",(TObject **) &ptof);
  getBranchObject("hitPgC",(TObject **) &treePgC);
  getBranchObject("hitAc",(TObject **) &treeAc);
  out->Branch("run",&run,"run/I");
  out->Branch("event",&event,"event/I");
  out->Branch("gapToF1",&gapToF1);
  out->Branch("gapToF2",&gapToF2);
  out->Branch("M2",&M2);
  out->Branch("vert_p",&vert_p);
  out->Branch("beta",&beta);
  out->Branch("Ep",&Ep);

  //  getBranchObject("RawPgcInfo",(TObject **) &treeRaw); 
  par = "kpm";
  //if(par=="kpm") M = 0.511;
//  else if(par=="kpm") M = 0.;

  string walk_parm =  Form("/data/trek/E36/gautam/crnRoot/cooker/inputFiles/tofM2/walkparm/ke3_walkparm.dat",par.Data());
  ifstream walkpar(walk_parm.c_str());
  if (!walkpar.is_open()){
    cout<<"cannot open the file"<<endl;
    exit(0);
  }
  walkpar>>fixed>>off>>u>>d>>AI>>AO>>BI>>BO>>u0>>d0>>AIo>>AOo>>BIo>>BOo>>p_max>>p_min;
  walkpar.close();
 

  int ix =0;
  string inputfile_A=Form("/data/trek/E36/gautam/crnRoot/cooker/inputFiles/tofM2/offset/meanA_kmu2_1.dat",par.Data());
  ifstream infileA(inputfile_A.c_str());
  if (!infileA.is_open()){
    cout<<"cannot open the file"<<endl;
    exit(0);
  }
  for(int i=0;i<12;i++){
    for (int k=-1;k<=1;k++)
      {
        ix  = (k+i)%12;
        if(ix>-1){
          infileA>>fixed>>i>>ix>>mean_a[i][ix];
	}
      }
  }

  string inputfile_B = Form("/data/trek/E36/gautam/crnRoot/cooker/inputFiles/tofM2/offset/meanB_kmu2_1.dat",par.Data());
  ifstream infileB(inputfile_B.c_str());
  if (!infileB.is_open()){
    cout<<"cannot open the file"<<endl;
    exit(0);
  }

  for(int i=0;i<12;i++){
    for (int k=-1;k<=1;k++)
      {
        ix  = (k+i)%12;
        if(ix>-1){
          infileB>>fixed>>i>>ix>>mean_b[i][ix];
        }
      }
  }
  infileB.close();
  infileA.close();

  string inputfile_C = Form("/data/trek/E36/gautam/crnRoot/cooker/inputFiles/tofM2/pgc_gains/PgC_Gain_align_peak_epPR.dat",par.Data());
  ifstream infileC(inputfile_C.c_str());
  if (!infileC.is_open()){
    cout<<"cannot open the file"<<endl;
    exit(0);
  }

  for(int i=0;i<12;i++){
    for (int k=0;k<7;k++)
      {
	infileC>>fixed>>i>>k>>peak_val[i][k];
	cout<<peak_val[i][k]<<" ";
	cout<<"     ";
      }

  }
  infileC.close();
 
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(1101);


  return 0;
};

Long_t Tof_M2::process()
{
 run=treeAc->run;
 event=treeAc->event;
 gapToF1.clear();
 gapToF2.clear();
 M2.clear();
 vert_p.clear();
  beta.clear();
  Ep.clear();
  gapToF1.push_back(ptr->fgapNumTof1);
  gapToF2.push_back(ptr->fgapNumTof2);

  fTof1TimeMean =0.;
  fTof2TimeMeanA =0.;
  fTof2TimeMeanB =0.;
  F_beta =0.;
  F_M2=0.;
  acAdcf = 0;
  acAdcu = 0;
  acAdcd = 0;
  acAdcfU[ptr->fgapNumTof2]=0;  
  acAdcfD[ptr->fgapNumTof2]=0;  
  AcAv[ptr->fgapNumTof2]=0;
  sum_beta_A[ptr->fgapNumTof2]=0.;
  sum_beta_B[ptr->fgapNumTof2]=0.;
  sum_M2_A[ptr->fgapNumTof2]=0.;
  sum_M2_B[ptr->fgapNumTof2]=0.;
  ADC_PgC[ptr->fgapNumTof2]=0;
  ADC_PgCf=0;
  Ep_PgCf=0;
  if(ptr->nTracks==1 &&  ptr->fVertSP <226. && 
     ptr->fVertSP>50 && ptr->fChi2PerNDF<10 && 
     0<ptr->fChi2PerNDF && ptr->fFlagToFMatch==kTRUE &&
     (ptr->fFlagTargetToAC==0 || ptr->fFlagTargetToAC==1 )){
    fTof1TimeMean = off + (ptof->timeUToF1 -  u/sqrt(ptof->adcUToF1 -u0) + ptof->timeDToF1 -  d/sqrt(ptof->adcDToF1 -d0))/2.;
    fTof2TimeMeanA = (ptof->timeAIToF2 - AI/sqrt(ptof->adcAIToF2 -AIo) + ptof->timeAOToF2 - AO/sqrt(ptof->adcAOToF2-AIo) )/2.;
    fTof2TimeMeanB = (ptof->timeBIToF2 - BI/sqrt(ptof->adcBIToF2 - BOo) + ptof->timeBOToF2 - BO/sqrt(ptof->adcBOToF2-BOo) )/2.;
  

    for(int q = -1;q<=1;q++){
      pndex = (ptr->fgapNumTof2 + q)%12; //index for tracking variables
      ppndex = (ptof->gapNumToF2 + q)%13;//index for cooked variables
      if( ppndex==-2 || ppndex==-1)continue;
      if((pndex == ptr->fgapNumTof1) &&  (ppndex == ptof->gapNumToF1) ) {
	time_tof1A[ppndex] = fTof1TimeMean;
	time_tof1B[ppndex] = fTof1TimeMean;
	fTof2TimeMeanA1[ppndex] = fTof2TimeMeanA;
	fTof2TimeMeanB1[ppndex] = fTof2TimeMeanB;
	if(fTof2TimeMeanA1[ppndex]>0 && fTof2TimeMeanA1[ppndex]<50 && ptr->fTof2PathTime>0 && ptr->fTof2PathTime<30){
	  fTof2PathTimeA1[pndex] = ptr->fTof2PathTime;
	  fTof1PathTime1A[pndex] = ptr->fTof1PathTime;
	  TOF_cal_A[ptr->fgapNumTof2][pndex] = (fTof2PathTimeA1[pndex] - fTof1PathTime1A[pndex]);
	  delta_corr_A[ptr->fgapNumTof2][pndex] = TOF_cal_A[ptr->fgapNumTof2][pndex] +(fTof2TimeMeanA1[ppndex] - time_tof1A[ppndex] - TOF_cal_A[ptr->fgapNumTof2][pndex] - mean_a[ptr->fgapNumTof2][pndex]) ;
	  beta_A[ptr->fgapNumTof2][pndex] = (ptr->fTof2PathLocation - ptr->fTof1PathLocation)/(3.*100.*delta_corr_A[ptr->fgapNumTof2][pndex]);
	  sum_beta_A[ptr->fgapNumTof2]+=beta_A[ptr->fgapNumTof2][pndex] ;
	  sum_M2_A[ptr->fgapNumTof2]+=ptr->fMwpc4SP*ptr->fMwpc4SP*(1./pow(beta_A[ptr->fgapNumTof2][pndex],2) - 1.);
	}
	if(fTof2TimeMeanB1[ppndex]>0 && fTof2TimeMeanB1[ppndex]<50 && ptr->fTof2PathTime>0 && ptr->fTof2PathTime<30){
	  fTof2PathTimeB1[pndex] = ptr->fTof2PathTime;
	  fTof1PathTime1B[pndex] = ptr->fTof1PathTime;
	  TOF_cal_B[ptr->fgapNumTof2][pndex] = (fTof2PathTimeB1[pndex] - fTof1PathTime1B[pndex]);
	  if(pndex>-1 && ppndex>0){
	    delta_corr_B[ptr->fgapNumTof2][pndex] = TOF_cal_B[ptr->fgapNumTof2][pndex] +  (fTof2TimeMeanB1[ppndex] - time_tof1B[ppndex] - TOF_cal_B[ptr->fgapNumTof2][pndex] - mean_b[ptr->fgapNumTof2][pndex]);
	    beta_B[ptr->fgapNumTof2][pndex] = (ptr->fTof2PathLocation - ptr->fTof1PathLocation)/(3.*100.*delta_corr_B[ptr->fgapNumTof2][pndex]);
	    sum_beta_B[ptr->fgapNumTof2]+=beta_B[ptr->fgapNumTof2][pndex];
	    sum_M2_B[ptr->fgapNumTof2]+=ptr->fMwpc4SP*ptr->fMwpc4SP*(1./pow(beta_B[ptr->fgapNumTof2][pndex],2) - 1.);
	  }
	}

      }
    }
    F_beta+=sum_beta_A[ptr->fgapNumTof2]+sum_beta_B[ptr->fgapNumTof2];
    F_M2+=sum_M2_A[ptr->fgapNumTof2]+sum_M2_B[ptr->fgapNumTof2];

    if( F_beta>0.4){
      fTofbeta = F_beta;
      fTofM2 = F_M2;

      if(fTofM2<3200){
	final_beta->Fill(fTofbeta);
	final_M2->Fill(fTofM2);
	theta_pzx->Fill(ptr->theta,ptr->p_zx);
	phi_pzx->Fill(ptr->phi,ptr->p_zx);
	p_pzx->Fill(ptr->fMwpc4SP,ptr->p_zx);
      
	vert_Mom->Fill(ptr->fVertSP);          
	final_beta->Fill(fTofbeta);
	hBetaGap->Fill(ptr->fgapNumTof2,fTofbeta);
         vert_p.push_back(ptr->fVertSP);
         M2.push_back(fTofM2);
         beta.push_back(fTofbeta);
         

	//Tof2 
	tof2_mom[ptr->fgapNumTof2] = ptr->fTof2SP;
	x_tof_max[ptr->fgapNumTof2] = x_tofc[ptr->fgapNumTof2] + xv_tof/2.;
	x_tof_min[ptr->fgapNumTof2] = x_tofc[ptr->fgapNumTof2] - xv_tof/2.;
	z_tofb[ptr->fgapNumTof2] = z_tofc[ptr->fgapNumTof2] + zv_tof/2.; //Back z

	pl_btof[ptr->fgapNumTof2] =  pathlength(ptr->fTof2SZ/10., z_tofb[ptr->fgapNumTof2], ptr->fTof2SNz);//Pathlength back
	x_btof[ptr->fgapNumTof2] = x_p(ptr->fTof2SX/10.,ptr->fTof2SNx, pl_btof[ptr->fgapNumTof2]);//x back 
	y_btof[ptr->fgapNumTof2] = y_p(ptr->fTof2SY/10.,ptr->fTof2SNy, pl_btof[ptr->fgapNumTof2]);//y back
	if(x_tof_min[ptr->fgapNumTof2]<x_btof[ptr->fgapNumTof2] && x_btof[ptr->fgapNumTof2]<x_tof_max[ptr->fgapNumTof2] && -yv_tof/2<y_btof[ptr->fgapNumTof2] && y_btof[ptr->fgapNumTof2]<yv_deg/2 )
	  {
	    ploss_tof_b[ptr->fgapNumTof2] = dp(ptr->fTof2SP, pl_btof[ptr->fgapNumTof2], Z_A_tof, MEE_tof, rho_tof);// front plane
	  }

	//deg
	x_deg_max[ptr->fgapNumTof2] = x_degc[ptr->fgapNumTof2] + xv_deg/2.;
	x_deg_min[ptr->fgapNumTof2] = x_degc[ptr->fgapNumTof2] - xv_deg/2.;

	z_degb[ptr->fgapNumTof2] = z_degc[ptr->fgapNumTof2] + zv_deg/2.; //Back z
	z_degf[ptr->fgapNumTof2] = z_degc[ptr->fgapNumTof2] - zv_deg/2.; //Back z

	pl_bdeg[ptr->fgapNumTof2] =  pathlength(ptr->fTof2SZ/10. + z_degc[ptr->fgapNumTof2] - z_tofc[ptr->fgapNumTof2], z_degb[ptr->fgapNumTof2], ptr->fTof2SNz);//Pathlength back
	pl_fdeg[ptr->fgapNumTof2] =  pathlength(ptr->fTof2SZ/10. + z_degc[ptr->fgapNumTof2] - z_tofc[ptr->fgapNumTof2], z_degf[ptr->fgapNumTof2], ptr->fTof2SNz);//Pathlength back
	x_bdeg[ptr->fgapNumTof2] = x_p(ptr->fTof2SX/10., ptr->fTof2SNx, pl_bdeg[ptr->fgapNumTof2]);//x back 
	y_bdeg[ptr->fgapNumTof2] = y_p(ptr->fTof2SY/10., ptr->fTof2SNy, pl_bdeg[ptr->fgapNumTof2]);//y back

	x_fdeg[ptr->fgapNumTof2] = x_p(ptr->fTof2SX/10., ptr->fTof2SNx, pl_fdeg[ptr->fgapNumTof2]);//x back 
	y_fdeg[ptr->fgapNumTof2] = y_p(ptr->fTof2SY/10., ptr->fTof2SNy, pl_fdeg[ptr->fgapNumTof2]);//y back

	if(x_deg_min[ptr->fgapNumTof2]<x_fdeg[ptr->fgapNumTof2]  && x_fdeg[ptr->fgapNumTof2]<x_deg_max[ptr->fgapNumTof2] && -yv_deg/2.<y_fdeg[ptr->fgapNumTof2] && y_fdeg[ptr->fgapNumTof2]<yv_deg/2.)
	  {
	    ploss_deg_f[ptr->fgapNumTof2] = dp(ptr->fTof2SP - TMath::Abs(ploss_tof_b[ptr->fgapNumTof2]), pl_fdeg[ptr->fgapNumTof2], Z_A_deg, MEE_deg,rho_deg);// front plane
	  }
	if(x_deg_min[ptr->fgapNumTof2]<x_bdeg[ptr->fgapNumTof2] && x_bdeg[ptr->fgapNumTof2]<x_deg_max[ptr->fgapNumTof2] && -yv_deg/2<y_bdeg[ptr->fgapNumTof2] && y_bdeg[ptr->fgapNumTof2]<yv_deg/2.)
	  {
	    ploss_deg_b[ptr->fgapNumTof2] = dp(ptr->fTof2SP - TMath::Abs(ploss_tof_b[ptr->fgapNumTof2]) - TMath::Abs(ploss_deg_f[ptr->fgapNumTof2]), pl_bdeg[ptr->fgapNumTof2], Z_A_deg, MEE_deg, rho_deg);// front plane
	  }
	ploss_tof_b[ptr->fgapNumTof2] = fabs(ploss_tof_b[ptr->fgapNumTof2]);
	ploss_deg_b[ptr->fgapNumTof2]= fabs(ploss_deg_b[ptr->fgapNumTof2]);
	ploss_deg_f[ptr->fgapNumTof2] = fabs(ploss_deg_f[ptr->fgapNumTof2]);
	//deg

	P_PgC[ptr->fgapNumTof2] = ptr->fTof2SP - TMath::Abs(ploss_tof_b[ptr->fgapNumTof2]) - 2.*TMath::Abs(ploss_deg_f[ptr->fgapNumTof2]);
	double fac_mom = 0.511/78.5;
	P_PgC[ptr->fgapNumTof2] = P_PgC[ptr->fgapNumTof2];
	h1Mom_PgC[ptr->fgapNumTof2]->Fill(tof2_mom[ptr->fgapNumTof2]);

	int sf_m =1 ;
	int cut_val =150;
	for(int i=0;i<12;i++){
	  ADC_PgC[i]=0;// ADC_PgCf=0;	
	  int count_j=0;

	  for(int ii=0;ii<7;ii++){
	    Array_ADC[i][ii]=treePgC->adcPgC[i][ii];    
	    if(i==11 && ii==5) Array_ADC[i][ii]=0;        
	    Array_ADC[i][ii]=Array_ADC[i][ii];
	  }
	  if(Array_ADC[i][0]>cut_val && Array_ADC[i][1]<cut_val && Array_ADC[i][2]<cut_val && Array_ADC[i][3]<cut_val && Array_ADC[i][4]<cut_val && Array_ADC[i][5]<cut_val && Array_ADC[i][6]<cut_val)  
	    {
	      count_j++;
	    
	      sf = sf_m/peak_val[i][0]; 
              array_ADC[i][0]= Array_ADC[i][0]*sf;

	      EP_PgCh[i][0]->Fill(array_ADC[i][0]/tof2_mom[i]);
	      ADC_PgCh[i][0]->Fill(array_ADC[i][0]);
	      ADC_PgC[i] +=array_ADC[i][0]; 

	    }
	  else if(Array_ADC[i][0]<cut_val && Array_ADC[i][1]>cut_val && Array_ADC[i][2]<cut_val &&Array_ADC[i][3]<cut_val &&Array_ADC[i][4]<cut_val &&Array_ADC[i][5]<cut_val &&Array_ADC[i][6]<cut_val )  
	    {
	      sf = sf_m/peak_val[i][1]; 
	      array_ADC[i][1]= Array_ADC[i][1]*sf;
	      EP_PgCh[i][1]->Fill(array_ADC[i][1]/tof2_mom[i]);
	      ADC_PgCh[i][1]->Fill(array_ADC[i][1]);
	      ADC_PgC[i] +=array_ADC[i][1]; 
	    }

	  else if(Array_ADC[i][0]<cut_val && Array_ADC[i][1]<cut_val && Array_ADC[i][2]>cut_val &&Array_ADC[i][3]<cut_val &&Array_ADC[i][4]<cut_val &&Array_ADC[i][5]<cut_val &&Array_ADC[i][6]<cut_val )  
	    {

	      sf = sf_m/peak_val[i][2]; 
	      array_ADC[i][2]= Array_ADC[i][2]*sf;
	      EP_PgCh[i][2]->Fill(array_ADC[i][2]/tof2_mom[i]);
	      ADC_PgCh[i][2]->Fill(array_ADC[i][2]);
	      ADC_PgC[i] +=array_ADC[i][2]; 
	    }

	  else if(Array_ADC[i][0]<cut_val && Array_ADC[i][1]<cut_val && Array_ADC[i][2]<cut_val &&Array_ADC[i][3]>cut_val &&Array_ADC[i][4]<cut_val &&Array_ADC[i][5]<cut_val &&Array_ADC[i][6]<cut_val )  
	    {
	      sf = sf_m/peak_val[i][3]; 
	      array_ADC[i][3]= Array_ADC[i][3]*sf;
	      EP_PgCh[i][3]->Fill(array_ADC[i][3]/tof2_mom[i]);
	      ADC_PgCh[i][3]->Fill(array_ADC[i][3]);
	      ADC_PgC[i] +=array_ADC[i][3]; 
	    }
  
	  else if(Array_ADC[i][0]<cut_val && Array_ADC[i][1]<cut_val && Array_ADC[i][2]<cut_val &&Array_ADC[i][3]<cut_val &&Array_ADC[i][4]>cut_val &&Array_ADC[i][5]<cut_val &&Array_ADC[i][6]<cut_val )  
	    {
	      sf = sf_m/peak_val[i][4]; 
              array_ADC[i][4]= Array_ADC[i][4]*sf;
              EP_PgCh[i][4]->Fill(array_ADC[i][4]/tof2_mom[i]);
	      ADC_PgCh[i][4]->Fill(array_ADC[i][4]);
	      ADC_PgC[i] +=array_ADC[i][4]; 
	    }

  
	  else if(Array_ADC[i][0]<cut_val && Array_ADC[i][1]<cut_val && Array_ADC[i][2]<cut_val &&Array_ADC[i][3]<cut_val &&Array_ADC[i][4]<cut_val &&Array_ADC[i][5]>cut_val &&Array_ADC[i][6]<cut_val )  
	    {
	      sf = sf_m/peak_val[i][5]; 
	      array_ADC[i][5]= Array_ADC[i][5]*sf;
	      EP_PgCh[i][5]->Fill(array_ADC[i][5]/tof2_mom[i]);
	      ADC_PgCh[i][5]->Fill(array_ADC[i][5]);
	      ADC_PgC[i] +=array_ADC[i][5]; 
	    }

	  else if(Array_ADC[i][0]<cut_val && Array_ADC[i][1]<cut_val && Array_ADC[i][2]<cut_val &&Array_ADC[i][3]<cut_val &&Array_ADC[i][4]<cut_val &&Array_ADC[i][5]<cut_val &&Array_ADC[i][6]>cut_val )  
	    {
	      sf = sf_m/peak_val[i][6]; 
	      array_ADC[i][6]= Array_ADC[i][6]*sf;
	      EP_PgCh[i][6]->Fill(array_ADC[i][6]/tof2_mom[i]);
	      ADC_PgCh[i][6]->Fill(array_ADC[i][6]);
	      ADC_PgC[i] +=array_ADC[i][6];   
	    }   

	  ////////////
	  else if(Array_ADC[2][6]>200 && array_ADC[2][6]/tof2_mom[2]>0.4) {
            sf = sf_m/peak_val[2][6];
            array_ADC[2][6]= Array_ADC[2][6]*sf;
            EP_PgCh[2][6]->Fill(array_ADC[2][6]/tof2_mom[2]);
            ADC_PgCh[2][6]->Fill(array_ADC[2][6]);
            ADC_PgC[2] +=array_ADC[2][6];
          }


          if(i==2 && ADC_PgC[2]>30){
	    if(ADC_PgC[2]/tof2_mom[2]>0.1) h1EP_PgC[2]->Fill(ADC_PgC[2]/tof2_mom[2]);
            h1adc_PgC[2]->Fill(ADC_PgC[2]);
            h1Mom_ADC[2]->Fill(ADC_PgC[2],tof2_mom[2]);
            h1Mom_ADC_TP[2]->Fill(ADC_PgC[2],tof2_mom[2]);
            ADC_PgCf += ADC_PgC[2];
	    if(ADC_PgC[2]/tof2_mom[2]>0.1) Ep_PgCf += ADC_PgC[2]/tof2_mom[2];

          }
          else if (i!=2){
	    if(ADC_PgC[i]/tof2_mom[i]>0.1) h1EP_PgC[i]->Fill(ADC_PgC[i]/tof2_mom[i]);
	    h1adc_PgC[i]->Fill(ADC_PgC[i]);
	    h1Mom_ADC[i]->Fill(ADC_PgC[i],tof2_mom[i]);
	    h1Mom_ADC_TP[i]->Fill(ADC_PgC[i],tof2_mom[i]);
	    ADC_PgCf += ADC_PgC[i]; //
	    if(ADC_PgC[i]/tof2_mom[i]>0.1)  Ep_PgCf += ADC_PgC[i]/tof2_mom[i];

	  }
	}


	if(ADC_PgCf>14)final_PgCAdc->Fill(ADC_PgCf);
	if(Ep_PgCf>0.13)final_PgCEp->Fill(Ep_PgCf);
         Ep.push_back(Ep_PgCf);
	int timeh = 2500;
	int timel = 500;
	int acGapU,acGapD;
	int acGapU1,acGapD1;
              //3 needed to be finalized step1
              //int ac_cutvld[12] = {120,195,170,170,100,80,110,130,310,80,130,170};
              //int ac_cutvlu[12] = {360,400,455,455,240,210,300,310,300,310,360,350};
              //4th gap peak final before gain alignment
              //int ac_cutvld[12] = {120,195,170,170,100,80,110,130,310,80,130,170};
              //int ac_cutvlu[12] = {360,400,455,327,240,210,300,310,300,310,360,350};	
              //4th aligned
              //int ac_cutvld[12] = {120,195,170,170,100,80,110,130,310,80,130,170};
             // int ac_cutvlu[12] = {360,400,455,365,240,210,300,310,300,310,360,350};	
              //Final 1
              //int ac_cutvld[12] = {120,220,170,170,100,80,110,130,120,80,130,170};
              //int ac_cutvlu[12] = {360,400,455,365,240,315,410,340,300,310,360,350};	
             //Final 2
             // int ac_cutvld[12] = {120,220,170,170,100,80,110,130,120,80,130,170};
             // int ac_cutvlu[12] = {360,400,455,365,240,315,460,460,300,310,360,350};	
             Double_t ac_dw[12] ={1155, 525, 945, 945, 525, 525, 735, 735,1155, 735, 735, 735};//875
	     Double_t ac_up[12] ={945,1155, 945, 525, 525,1155, 945,1365,1155, 945, 735, 735};//840 
  
              int ac_cutvld[12] = {120,250,120,170,100,70,110,10,70,90,130,140};
              int ac_cutvlu[12] = {310,400,455,360,240,315,450,490,260,130,360,350};	
	    
	if (ADC_PgCf>20 && Ep_PgCf>0.1) {
	  for(int i = 0; i < 12; i++){

	   ac_dsf[i]=817./ac_dw[i];
	    ac_usf[i]=817./ac_up[i]; 


	    if(i==0 && treeAc->adcACU[i]>ac_cutvlu[i]  && treeAc->adcACD[i]<ac_cutvlu[i] && (treeAc->timeACU[i]>timel && treeAc->timeACU[i]<timeh)){
	      acGapU=i;
	      if(i==ptr->fgapNumTof2 && (ptr->fgapNumTof2== ptr->fgapNumTof1+1||ptr->fgapNumTof2== ptr->fgapNumTof1||ptr->fgapNumTof2== ptr->fgapNumTof1-1) && ((treeAc->timeACU[i]>timel && treeAc->timeACU[i]<timeh))){
	        ac_usf[i] =ac_usf[i];
		h1adc_ACu[i]->Fill(treeAc->adcACU[i]*ac_usf[i]);
		h1tdc_ACu[i]->Fill(treeAc->timeACU[i]);
		acAdcu+= treeAc->adcACU[i]*ac_usf[i];
		acGapU1=i;
	      }
	    }
	    else if(i==1 && treeAc->adcACU[i]>ac_cutvlu[i] && treeAc->adcACD[i]<ac_cutvlu[i] && (treeAc->timeACU[i]>timel && treeAc->timeACU[i]<timeh)){
	      acGapU=i;
	      if(i==ptr->fgapNumTof2 && (ptr->fgapNumTof2== ptr->fgapNumTof1+1||ptr->fgapNumTof2== ptr->fgapNumTof1||ptr->fgapNumTof2== ptr->fgapNumTof1-1)){
	        ac_usf[i] = ac_usf[i]*817./525.*817./945.;
         	h1adc_ACu[i]->Fill(treeAc->adcACU[i]*ac_usf[i]);
		h1tdc_ACu[i]->Fill(treeAc->timeACU[i]);
		acAdcu+= treeAc->adcACU[i]*ac_usf[i];
		acGapU1=i;
	      }
	    }


	    else if(i==2 && treeAc->adcACU[i]>ac_cutvlu[i] && treeAc->adcACD[i]<ac_cutvlu[i] && (treeAc->timeACU[i]>timel && treeAc->timeACU[i]<timeh)){
	      acGapU=i;
	      if( i==ptr->fgapNumTof2 && (ptr->fgapNumTof2== ptr->fgapNumTof1+1||ptr->fgapNumTof2== ptr->fgapNumTof1||ptr->fgapNumTof2== ptr->fgapNumTof1-1)){
	        ac_usf[i] = ac_usf[i]*817./945.*817./945.;
        	h1adc_ACu[i]->Fill(treeAc->adcACU[i]*ac_usf[i]);
		h1tdc_ACu[i]->Fill(treeAc->timeACU[i]);
		acAdcu+= treeAc->adcACU[i]*ac_usf[i];

		acGapU1=i;
	      }
	    }
	    else if(i==3 && treeAc->adcACU[i]>ac_cutvlu[i] && treeAc->adcACD[i]<ac_cutvlu[i] && (treeAc->timeACU[i]>timel && treeAc->timeACU[i]<timeh)){
	      acGapU=i;
	      if(i==ptr->fgapNumTof2 && (ptr->fgapNumTof2== ptr->fgapNumTof1+1||ptr->fgapNumTof2== ptr->fgapNumTof1||ptr->fgapNumTof2== ptr->fgapNumTof1-1)){
	        ac_usf[i] = ac_usf[i]*817./945.*817./945.*817./945.*817./525.;
         	h1adc_ACu[i]->Fill(treeAc->adcACU[i]*ac_usf[i]);
		h1tdc_ACu[i]->Fill(treeAc->timeACU[i]);
		acAdcu+= treeAc->adcACU[i]*ac_usf[i];
		acGapU1=i;
	      }
	    }

	    else if(i==4 && treeAc->adcACU[i]>ac_cutvlu[i] && treeAc->adcACD[i]<ac_cutvlu[i] && (treeAc->timeACU[i]>timel && treeAc->timeACU[i]<timeh)){
	      acGapU=i;
	      if( i==ptr->fgapNumTof2 && (ptr->fgapNumTof2== ptr->fgapNumTof1+1||ptr->fgapNumTof2== ptr->fgapNumTof1||ptr->fgapNumTof2== ptr->fgapNumTof1-1)){
	        ac_usf[i] = ac_usf[i]*817./525.;
		h1adc_ACu[i]->Fill(treeAc->adcACU[i]*ac_usf[i]);
		h1tdc_ACu[i]->Fill(treeAc->timeACU[i]);
		acAdcu+= treeAc->adcACU[i]*ac_usf[i];
		acGapU1=i;
	      }

	    }


	    else if(i==5 && treeAc->adcACU[i]>ac_cutvlu[i] && treeAc->adcACD[i]<ac_cutvlu[i] && (treeAc->timeACU[i]>timel && treeAc->timeACU[i]<timeh)){
	      acGapU=i;
	      if(i==ptr->fgapNumTof2 && (ptr->fgapNumTof2== ptr->fgapNumTof1+1||ptr->fgapNumTof2== ptr->fgapNumTof1||ptr->fgapNumTof2== ptr->fgapNumTof1-1)){
	        ac_usf[i] = ac_usf[i]*817./525.*817./1155.*817./525.*817./945.;
		h1adc_ACu[i]->Fill(treeAc->adcACU[i]*ac_usf[i]);
		h1tdc_ACu[i]->Fill(treeAc->timeACU[i]);
		acAdcu+= treeAc->adcACU[i]*ac_usf[i];
		acGapU1=i;
	      }

	    }


	    else if(i==6 && treeAc->adcACU[i]>ac_cutvlu[i] && treeAc->adcACD[i]<ac_cutvlu[i] && (treeAc->timeACU[i]>timel && treeAc->timeACU[i]<timeh)){
	      acGapU=i;
	      if(i==ptr->fgapNumTof2 && (ptr->fgapNumTof2== ptr->fgapNumTof1+1||ptr->fgapNumTof2== ptr->fgapNumTof1||ptr->fgapNumTof2== ptr->fgapNumTof1-1)){
	        ac_usf[i] = ac_usf[i];
		h1adc_ACu[i]->Fill(treeAc->adcACU[i]*ac_usf[i]);
		h1tdc_ACu[i]->Fill(treeAc->timeACU[i]);
		acAdcu+= treeAc->adcACU[i]*ac_usf[i];
		acGapU1=i;
	      }

	    }

	    else if(i==7 && treeAc->adcACU[i]>ac_cutvlu[i] && treeAc->adcACD[i]<ac_cutvlu[i] && (treeAc->timeACU[i]>timel && treeAc->timeACU[i]<timeh)){
	      acGapU=i;
	      if(i==ptr->fgapNumTof2 && (ptr->fgapNumTof2== ptr->fgapNumTof1+1||ptr->fgapNumTof2== ptr->fgapNumTof1||ptr->fgapNumTof2== ptr->fgapNumTof1-1)){
	        ac_usf[i] = ac_usf[i]*817./525.*817./945.*817./945.;
		h1adc_ACu[i]->Fill(treeAc->adcACU[i]*ac_usf[i]);
		h1tdc_ACu[i]->Fill(treeAc->timeACU[i]);
		acAdcu+= treeAc->adcACU[i]*ac_usf[i];
		acGapU1=i;
	      }
	    }

	    else if(i==8 && treeAc->adcACU[i]>ac_cutvlu[i] && treeAc->adcACD[i]<ac_cutvlu[i] && (treeAc->timeACU[i]>timel && treeAc->timeACU[i]<timeh)){
	      acGapU=i;
	      if(i==ptr->fgapNumTof2 && (ptr->fgapNumTof2== ptr->fgapNumTof1+1||ptr->fgapNumTof2== ptr->fgapNumTof1||ptr->fgapNumTof2== ptr->fgapNumTof1-1)){
	        ac_usf[i] =ac_usf[i]*817./945.;
		h1adc_ACu[i]->Fill(treeAc->adcACU[i]*ac_usf[i]);
		h1tdc_ACu[i]->Fill(treeAc->timeACU[i]);
		acAdcu+= treeAc->adcACU[i]*ac_usf[i];
		acGapU1=i;
	      }
	    }

	    else if(i==9 && treeAc->adcACU[i]>ac_cutvlu[i] && treeAc->adcACD[i]<ac_cutvlu[i] && (treeAc->timeACU[i]>timel && treeAc->timeACU[i]<timeh)){
	      acGapU=i;
	      if(i==ptr->fgapNumTof2 && (ptr->fgapNumTof2== ptr->fgapNumTof1+1||ptr->fgapNumTof2== ptr->fgapNumTof1||ptr->fgapNumTof2== ptr->fgapNumTof1-1)){
	        ac_usf[i] = ac_usf[i];

		h1adc_ACu[i]->Fill(treeAc->adcACU[i]*ac_usf[i]);
		h1tdc_ACu[i]->Fill(treeAc->timeACU[i]);
		acAdcu+= treeAc->adcACU[i]*ac_usf[i];
		acGapU1=i;
	      }
	    }

	    else if(i==10 && treeAc->adcACU[i]>ac_cutvlu[i] && treeAc->adcACD[i]<ac_cutvlu[i] && (treeAc->timeACU[i]>timel && treeAc->timeACU[i]<timeh)){
	      acGapU=i;
	      if(i==ptr->fgapNumTof2 && (ptr->fgapNumTof2== ptr->fgapNumTof1+1||ptr->fgapNumTof2== ptr->fgapNumTof1||ptr->fgapNumTof2== ptr->fgapNumTof1-1)){
	        ac_usf[i] = ac_usf[i]*817./945.;
		h1adc_ACu[i]->Fill(treeAc->adcACU[i]*ac_usf[i]);
		h1tdc_ACu[i]->Fill(treeAc->timeACU[i]);
		acAdcu+= treeAc->adcACU[i]*ac_usf[i];
		acGapU1=i;
	      }
	    }
	    else if(i==11 && treeAc->adcACU[i]>ac_cutvlu[i] && treeAc->adcACD[i]<ac_cutvlu[i] && (treeAc->timeACU[i]>timel && treeAc->timeACU[i]<timeh)){
	      acGapU=i;
	      if(i==ptr->fgapNumTof2 && (ptr->fgapNumTof2== ptr->fgapNumTof1+1||ptr->fgapNumTof2== ptr->fgapNumTof1||ptr->fgapNumTof2== ptr->fgapNumTof1-1)){
	        ac_usf[i] = ac_usf[i];
		h1adc_ACu[i]->Fill(treeAc->adcACU[i]*ac_usf[i]);
		h1tdc_ACu[i]->Fill(treeAc->timeACU[i]);
		acAdcu+= treeAc->adcACU[i]*ac_usf[i];
		acGapU1=i;
	      }
	    }
	    else if(treeAc->adcACU[i]>ac_cutvlu[i] && treeAc->adcACD[i]<ac_cutvlu[i] && (treeAc->timeACU[i]>timel && treeAc->timeACU[i]<timeh)){
	      acGapU=i;
	      if(i==ptr->fgapNumTof2 && (ptr->fgapNumTof2== ptr->fgapNumTof1+1||ptr->fgapNumTof2== ptr->fgapNumTof1||ptr->fgapNumTof2== ptr->fgapNumTof1-1)){ 
		h1adc_ACu[i]->Fill(treeAc->adcACU[i]*ac_usf[i]);
		h1tdc_ACu[i]->Fill(treeAc->timeACU[i]);
		acAdcu+= treeAc->adcACU[i]*ac_usf[i]; 
		acGapU1=i;
	      }
	    }

	    //Downstream

	    else if(i==0 && treeAc->adcACU[i]<ac_cutvld[i] && treeAc->adcACD[i]>ac_cutvld[i] && (treeAc->timeACD[i]>timel && treeAc->timeACD[i]<timeh)){
	      acGapD=i;
	      if( i==ptr->fgapNumTof2 && (ptr->fgapNumTof2== ptr->fgapNumTof1-1 || ptr->fgapNumTof2== ptr->fgapNumTof1||ptr->fgapNumTof2== ptr->fgapNumTof1+1)){ 
                //ac_dsf[i] = ac_dsf[i];
         	h1adc_ACd[i]->Fill(treeAc->adcACD[i]*ac_dsf[i]);
		h1tdc_ACd[i]->Fill(treeAc->timeACD[i]);
		acAdcd+= treeAc->adcACD[i]*ac_dsf[i];
		acGapD1=i;
	      }
	    }
	    else if(i==1 && treeAc->adcACU[i]<ac_cutvld[i] && treeAc->adcACD[i]>ac_cutvld[i] && (treeAc->timeACD[i]>timel && treeAc->timeACD[i]<timeh)){
	      acGapD=i;
	      if( i==ptr->fgapNumTof2 && (ptr->fgapNumTof2== ptr->fgapNumTof1-1 || ptr->fgapNumTof2== ptr->fgapNumTof1||ptr->fgapNumTof2== ptr->fgapNumTof1+1)){ 
                ac_dsf[i] = ac_dsf[i]*817./945.;
         	h1adc_ACd[i]->Fill(treeAc->adcACD[i]*ac_dsf[i]);
		h1tdc_ACd[i]->Fill(treeAc->timeACD[i]);
		acAdcd+= treeAc->adcACD[i]*ac_dsf[i];
		acGapD1=i;
	      }
	    }


	    else if(i==2 && treeAc->adcACU[i]<ac_cutvld[i] && treeAc->adcACD[i]>ac_cutvld[i] && (treeAc->timeACD[i]>timel && treeAc->timeACD[i]<timeh)){
	      acGapD=i;
	      if(i==ptr->fgapNumTof2 && (ptr->fgapNumTof2== ptr->fgapNumTof1-1 || ptr->fgapNumTof2== ptr->fgapNumTof1||ptr->fgapNumTof2== ptr->fgapNumTof1+1)){ 
		h1adc_ACd[i]->Fill(treeAc->adcACD[i]*ac_dsf[i]);
		h1tdc_ACd[i]->Fill(treeAc->timeACD[i]);
		acAdcd+= treeAc->adcACD[i]*ac_dsf[i];
		acGapD1=i;
	      }
	    }


	    else if(i==3 && treeAc->adcACU[i]<ac_cutvld[i] && treeAc->adcACD[i]>ac_cutvld[i] && (treeAc->timeACD[i]>timel && treeAc->timeACD[i]<timeh)){
	      acGapD=i;
	      if(i==ptr->fgapNumTof2 && (ptr->fgapNumTof2== ptr->fgapNumTof1-1 || ptr->fgapNumTof2== ptr->fgapNumTof1||ptr->fgapNumTof2== ptr->fgapNumTof1+1)){ 
		h1adc_ACd[i]->Fill(treeAc->adcACD[i]*ac_dsf[i]);
		h1tdc_ACd[i]->Fill(treeAc->timeACD[i]);
		acAdcd+= treeAc->adcACD[i]*ac_dsf[i];
		acGapD1=i;
	      }
	    }
	    else if(i==4 && treeAc->adcACU[i]<ac_cutvld[i] && treeAc->adcACD[i]>ac_cutvld[i] && (treeAc->timeACD[i]>timel && treeAc->timeACD[i]<timeh)){
	      acGapD=i;
	      if(i==ptr->fgapNumTof2 && (ptr->fgapNumTof2== ptr->fgapNumTof1-1 || ptr->fgapNumTof2== ptr->fgapNumTof1||ptr->fgapNumTof2== ptr->fgapNumTof1+1)){ 
		h1adc_ACd[i]->Fill(treeAc->adcACD[i]*ac_dsf[i]);
		h1tdc_ACd[i]->Fill(treeAc->timeACD[i]);
		acAdcd+= treeAc->adcACD[i]*ac_dsf[i];
		acGapD1=i;
	      }
	    }


	    else if(i==5 && treeAc->adcACU[i]<ac_cutvld[i] && treeAc->adcACD[i]>ac_cutvld[i] && (treeAc->timeACD[i]>timel && treeAc->timeACD[i]<timeh)){
	      acGapD=i;
	      if(i==ptr->fgapNumTof2 && (ptr->fgapNumTof2== ptr->fgapNumTof1-1 || ptr->fgapNumTof2== ptr->fgapNumTof1||ptr->fgapNumTof2== ptr->fgapNumTof1+1)){ 
		 ac_dsf[i] = ac_dsf[i]*817./945.*817./945.*857./1050.;
                h1adc_ACd[i]->Fill(treeAc->adcACD[i]*ac_dsf[i]);
		h1tdc_ACd[i]->Fill(treeAc->timeACD[i]);
		acAdcd+= treeAc->adcACD[i]*ac_dsf[i];
		acGapD1=i;
	      }
	    }


	    else if(i==6 && treeAc->adcACU[i]<ac_cutvld[i] && treeAc->adcACD[i]>ac_cutvld[i] && (treeAc->timeACD[i]>timel && treeAc->timeACD[i]<timeh)){
	      acGapD=i;
	      if(i==ptr->fgapNumTof2 && (ptr->fgapNumTof2== ptr->fgapNumTof1-1 || ptr->fgapNumTof2== ptr->fgapNumTof1||ptr->fgapNumTof2== ptr->fgapNumTof1+1)){ 
		 ac_dsf[i] = ac_dsf[i]*817./945.;
		h1adc_ACd[i]->Fill(treeAc->adcACD[i]*ac_dsf[i]);
		h1tdc_ACd[i]->Fill(treeAc->timeACD[i]);
		acAdcd+= treeAc->adcACD[i]*ac_dsf[i];
		acGapD1=i;
	      }
	    }


	    else if(i==7 && treeAc->adcACU[i]<ac_cutvld[i] && treeAc->adcACD[i]>ac_cutvld[i] && (treeAc->timeACD[i]>timel && treeAc->timeACD[i]<timeh)){
	      acGapD=i;
	      if(i==ptr->fgapNumTof2 && (ptr->fgapNumTof2== ptr->fgapNumTof1-1 || ptr->fgapNumTof2== ptr->fgapNumTof1||ptr->fgapNumTof2== ptr->fgapNumTof1+1)){ 
		 ac_dsf[i] = ac_dsf[i]*817./945.;
		h1adc_ACd[i]->Fill(treeAc->adcACD[i]*ac_dsf[i]);
		h1tdc_ACd[i]->Fill(treeAc->timeACD[i]);
		acAdcd+= treeAc->adcACD[i]*ac_dsf[i];
		acGapD1=i;
	      }
	    }

	    else if(i==8 && treeAc->adcACU[i]<ac_cutvld[i] && treeAc->adcACD[i]>ac_cutvld[i] && (treeAc->timeACD[i]>timel && treeAc->timeACD[i]<timeh)){
	      acGapD=i;
	      if(i==ptr->fgapNumTof2 && (ptr->fgapNumTof2== ptr->fgapNumTof1-1 || ptr->fgapNumTof2== ptr->fgapNumTof1||ptr->fgapNumTof2== ptr->fgapNumTof1+1)){ 
		ac_dsf[i] =ac_dsf[i];
		h1adc_ACd[i]->Fill(treeAc->adcACD[i]*ac_dsf[i]);
		h1tdc_ACd[i]->Fill(treeAc->timeACD[i]);
		acAdcd+= treeAc->adcACD[i]*ac_dsf[i];
		acGapD1=i;
	      }
	    }


	    else if(i==9 && treeAc->adcACU[i]<ac_cutvld[i] && treeAc->adcACD[i]>ac_cutvld[i] && (treeAc->timeACD[i]>timel && treeAc->timeACD[i]<timeh)){
	      acGapD=i;
	      if(i==ptr->fgapNumTof2 && (ptr->fgapNumTof2== ptr->fgapNumTof1-1 || ptr->fgapNumTof2== ptr->fgapNumTof1||ptr->fgapNumTof2== ptr->fgapNumTof1+1)){ 
		 ac_dsf[i] = ac_dsf[i]*817./1155.;
		h1adc_ACd[i]->Fill(treeAc->adcACD[i]*ac_dsf[i]);
		h1tdc_ACd[i]->Fill(treeAc->timeACD[i]);
		acAdcd+= treeAc->adcACD[i]*ac_dsf[i];
		acGapD1=i;
	      }
	    }


	    else if(i==10 && treeAc->adcACU[i]<ac_cutvld[i] && treeAc->adcACD[i]>ac_cutvld[i] && (treeAc->timeACD[i]>timel && treeAc->timeACD[i]<timeh)){
	      acGapD=i;
	      if(i==ptr->fgapNumTof2 && (ptr->fgapNumTof2== ptr->fgapNumTof1-1 || ptr->fgapNumTof2== ptr->fgapNumTof1||ptr->fgapNumTof2== ptr->fgapNumTof1+1)){ 
		 //ac_dsf[i] = ac_dsf[i]*817./945.;
		 ac_dsf[i] = ac_dsf[i]*817./945. * 857./1050.;
		h1adc_ACd[i]->Fill(treeAc->adcACD[i]*ac_dsf[i]);
		h1tdc_ACd[i]->Fill(treeAc->timeACD[i]);
		acAdcd+= treeAc->adcACD[i]*ac_dsf[i];
		acGapD1=i;
	      }
	    }

	    else if(i==11 && treeAc->adcACU[i]<ac_cutvld[i] && treeAc->adcACD[i]>ac_cutvld[i] && (treeAc->timeACD[i]>timel && treeAc->timeACD[i]<timeh)){
	      acGapD=i;
	      if(i==ptr->fgapNumTof2 && (ptr->fgapNumTof2== ptr->fgapNumTof1-1 || ptr->fgapNumTof2== ptr->fgapNumTof1||ptr->fgapNumTof2== ptr->fgapNumTof1+1)){ 
		h1adc_ACd[i]->Fill(treeAc->adcACD[i]*ac_dsf[i]*857./1050.);
		h1tdc_ACd[i]->Fill(treeAc->timeACD[i]);
		acAdcd+= treeAc->adcACD[i]*ac_dsf[i]*857./1050.;
		acGapD1=i;
	      }
	    }

	    else if(treeAc->adcACU[i]<ac_cutvld[i] && treeAc->adcACD[i]>ac_cutvld[i] && (treeAc->timeACD[i]>timel && treeAc->timeACD[i]<timeh)){
	      acGapD=i;
	      if(i==ptr->fgapNumTof2 && (ptr->fgapNumTof2== ptr->fgapNumTof1-1 || ptr->fgapNumTof2== ptr->fgapNumTof1||ptr->fgapNumTof2== ptr->fgapNumTof1+1)){ 
		h1adc_ACd[i]->Fill(treeAc->adcACD[i]*ac_dsf[i]);
		h1tdc_ACd[i]->Fill(treeAc->timeACD[i]);
		acAdcd+= treeAc->adcACD[i]*ac_dsf[i];
		acGapD1=i;
	      }
	    }
	    else continue;
	  }
	  h1adc_acGapuGap->Fill(ptr->fgapNumTof2,acGapU);
	  h1adc_acGapuGap2->Fill(ptr->fgapNumTof2,acGapU1);
	  h1adc_acGapdGap->Fill(ptr->fgapNumTof2,acGapD);
	  h1adc_acGapdGap2->Fill(ptr->fgapNumTof2,acGapD1);
	  h1adc_acuGap->Fill(ptr->fgapNumTof2,acAdcu);     
	  if(acAdcu>1)final_acAdcU->Fill(acAdcu);
	  if(acAdcd>1)final_acAdcD->Fill(acAdcd);
	  if((acAdcu+acAdcd)>1)final_acAdc->Fill((acAdcu+acAdcd));
	}    

      }
    }
  }

  return 0;
};

Long_t Tof_M2::done()
{
  gStyle->SetOptFit(1101);
  //gStyle->SetOptStat(0);
  
  TCanvas* c1_com = new TCanvas("c1_com", " AC com ", 1400,900);
  c1_com->Divide(4,3);
  TCanvas* c1 = new TCanvas("c1", " AC ADC Ped Sub. Upstream ", 1400,900);
  c1->Divide(4,3);
  TCanvas* c2 = new TCanvas("c2", " AC ADC Ped Sub. Downstream ", 1400,900);
  c2->Divide(4,3);
  TCanvas* c1tdc = new TCanvas("c1tdc", " AC TDC Ped Sub. Upstream ", 1400,900);
  c1tdc->Divide(4,3);
  TCanvas* c2tdc = new TCanvas("c2tdc", " AC TDC Ped Sub. Downstream ", 1400,900);
  c2tdc->Divide(4,3);



  TCanvas* c1_uGap = new TCanvas("c1_uGap", " AC ADC Up vs gap ", 1400,900);
  c1_uGap->Divide(4,3);
  TCanvas* c1_dGap = new TCanvas("c1_dGap", " AC ADC Down vs gap ", 1400,900);
  c1_dGap->Divide(4,3);
  TCanvas* c1_ac = new TCanvas("c1_ac", " AC ADC Ped Sub. Upstream ", 100, 100, 800, 650);
  c1_ac->Divide(3,3);
  c1_ac->cd(1);  
  final_acAdcU->Draw("HIST");
  final_acAdcU->GetXaxis()->SetTitle("ADC sum Upstream");
  final_acAdcU->GetXaxis()->CenterTitle();
  c1_ac->cd(2);  
  final_acAdcD->Draw("HIST");
  final_acAdcD->GetXaxis()->SetTitle("ADC sum Downstream");
  final_acAdcD->GetXaxis()->CenterTitle();
  c1_ac->cd(3);  
  final_acAdc->Draw("HIST");
  final_acAdc->GetXaxis()->SetTitle("ADC sum ");
  final_acAdc->GetXaxis()->CenterTitle();
  c1_ac->cd(4);  
  h1adc_acuGap->Draw("COLZ"); 
  h1adc_acuGap->GetXaxis()->SetTitle("ADC sum Upstream ");
  h1adc_acuGap->GetXaxis()->CenterTitle();
  h1adc_acuGap->GetYaxis()->SetTitle("Tof2 Gap Number ");
  h1adc_acuGap->GetYaxis()->CenterTitle();
  c1_ac->cd(5);  
  h1adc_acGapuGap->Draw("COLZ");
  h1adc_acGapuGap->GetYaxis()->SetTitle("AC gap upstream ");
  h1adc_acGapuGap->GetYaxis()->CenterTitle();
  h1adc_acGapuGap->GetXaxis()->SetTitle("Tof2 Gap Number ");
  h1adc_acGapuGap->GetXaxis()->CenterTitle();
  c1_ac->cd(6);  
  h1adc_acGapdGap->Draw("COLZ");
  h1adc_acGapdGap->GetYaxis()->SetTitle("AC gap Downstream ");
  h1adc_acGapdGap->GetYaxis()->CenterTitle();
  h1adc_acGapdGap->GetXaxis()->SetTitle("Tof2 Gap Number ");
  h1adc_acGapdGap->GetXaxis()->CenterTitle();
  c1_ac->cd(7);  
  h1adc_acGapuGap2->Draw("COLZ");
  h1adc_acGapuGap2->GetYaxis()->SetTitle("AC gap upstream (");
  h1adc_acGapuGap2->GetYaxis()->CenterTitle();
  h1adc_acGapuGap2->GetXaxis()->SetTitle("Tof2 Gap Number ");
  h1adc_acGapuGap2->GetXaxis()->CenterTitle();
  c1_ac->cd(8);  
  h1adc_acGapdGap2->Draw("COLZ");
  h1adc_acGapdGap2->GetYaxis()->SetTitle("AC gap Downstream ");
  h1adc_acGapdGap2->GetYaxis()->CenterTitle();
  h1adc_acGapdGap2->GetXaxis()->SetTitle("Tof2 Gap Number ");
  h1adc_acGapdGap2->GetXaxis()->CenterTitle();

  TCanvas* c3 = new TCanvas("c3", " all In one final", 1900, 1000);
  c3->Divide(3,2);
  TCanvas* c4 = new TCanvas("c4", " PGC ADC 12 gaps ", 100, 100, 800, 650);
  c4->Divide(3,4);
  TCanvas* c1_pPgC = new TCanvas("c1_pPgC", " p PgC", 100, 100, 800, 650);
  c1_pPgC->Divide(3,4);
  TCanvas* c5 = new TCanvas("c5", " Ene Adc sum", 100, 100, 800, 650);
  c5->Divide(3,4); 
  TCanvas* c5tp = new TCanvas("c5tp", " Ene Adc sum", 100, 100, 800, 650);
  c5tp->Divide(3,4); 
  TCanvas* cEp1 = new TCanvas("cEp1", "cEp1", 100, 100, 800, 650);
  cEp1->Divide(3,4); 

  TCanvas* c1_AcUD = new TCanvas("c1_AcUD", "c1_AcUD", 1400,900);
  c1_AcUD->Divide(3,4); 


 
  c3->cd(1);
  vert_Mom->Draw("HIST");
  vert_Mom->GetXaxis()->SetTitle("Vert Momemntum (MeV)");
  vert_Mom->GetXaxis()->CenterTitle();

  c3->cd(2);
  final_beta->Draw("HIST");
  final_beta->GetXaxis()->SetTitle("#beta");
  final_beta->GetXaxis()->CenterTitle();
  c3->cd(3);
  final_M2->Draw("HIST");
  final_M2->GetXaxis()->SetTitle("M^{2}(MeV/c^{2})^{2}");
  final_M2->GetXaxis()->CenterTitle();
  int max_M2 = final_M2->GetMaximum();
  TLine *line = new TLine(3000,0,3000,max_M2);
  line->SetLineColor(kRed);
  //  line->Draw();
  c3->cd(4);
  final_PgCAdc->Draw();
  final_PgCAdc->GetXaxis()->SetTitle("ADC sum");
  final_PgCAdc->GetXaxis()->CenterTitle();
  c3->cd(5);
  final_PgCEp->Draw("HIST");
  final_PgCEp->GetXaxis()->SetTitle("E/p");
  final_PgCEp->GetXaxis()->CenterTitle();
  c3->cd(6);
  hBetaGap->Draw("colz");
  hBetaGap->GetXaxis()->SetTitle("Gap Number");
  hBetaGap->GetXaxis()->CenterTitle();
  hBetaGap->GetYaxis()->SetTitle("#beta");
  hBetaGap->GetYaxis()->CenterTitle();
  TCanvas *c11[12];
  TCanvas *c12[12];
 FILE *AC_U = fopen(Form("/data/trek/E36/gautam/crnRoot/cooker/inputFiles/tofM2/pgc_gains/AcU_peak_%s.dat",par.Data()),"w");
  FILE *AC_D = fopen(Form("/data/trek/E36/gautam/crnRoot/cooker/inputFiles/tofM2/pgc_gains/AcD_peak_%s.dat",par.Data()),"w");
  for(int i = 0; i < 12; i++){
    c1_pPgC->cd(i+1);
    h1Mom_PgC[i]->Draw("HIST");

    c1_com->cd(i+1);
    h1AcAv[i]->Draw("HIST"); 
    c1->cd(i+1);
    h1adc_ACu[i]->Draw("HIST"); 
    Ac_U[i] =  h1adc_ACu[i]->GetBinCenter(h1adc_ACu[i]->GetMaximumBin());;
  fprintf(AC_U,"%4.0f,",Ac_U[i]);  
    c2->cd(i+1);
    h1adc_ACd[i]->Draw("HIST"); 
    Ac_D[i] =  h1adc_ACd[i]->GetBinCenter(h1adc_ACd[i]->GetMaximumBin());;
   fprintf(AC_D,"%4.0f,", Ac_D[i]);

    c1tdc->cd(i+1);
    h1tdc_ACu[i]->Draw("HIST"); 
  
    c2tdc->cd(i+1);
    h1tdc_ACd[i]->Draw("HIST"); 

    c1_AcUD->cd(i+1);
    h2adc_ACUD[i]->Draw("colz");
    h2adc_ACUD[i]->GetXaxis()->SetTitle(Form("ADC Ac Up Gap %d",i+1));
    h2adc_ACUD[i]->GetXaxis()->CenterTitle();
    h2adc_ACUD[i]->GetYaxis()->SetTitle(Form("ADC Ac Down Gap %d",i+1));
    h2adc_ACUD[i]->GetYaxis()->CenterTitle();


    c4->cd(i+1);
    h1adc_PgC[i]->Draw("HIST");  

    c1_uGap->cd(i+1);
    h1adc_ACuGap[i]->Draw("COLZ");  
    c1_dGap->cd(i+1);
    h1adc_ACdGap[i]->Draw("COLZ");  
 
    c11[i] = new TCanvas(Form("c11[%d]",i));
    c12[i] = new TCanvas(Form("c12[%d]",i));
    c11[i]->Divide(4,2);
    c12[i]->Divide(4,2);
    for(int ii = 0; ii < 7; ii++){
      c11[i]->cd(ii+1);
      ADC_PgCh[i][ii]->Draw("HIST");
      c12[i]->cd(ii+1);
      EP_PgCh[i][ii]->Draw("HIST");
    }
    c5tp->cd(i+1);
    h1Mom_ADC_TP[i]->Draw("cloz");
    h1Mom_ADC_TP[i]->Fit("pol1","","",0,300);
    c5->cd(i+1);
    h1Mom_ADC[i]->Draw("cloz");
    c11[i]->Write(Form("PgC_ADC_%d",i));
    c12[i]->Write(Form("PgC1_Ep_%d",i));
     
    cEp1->cd(i+1);
    h1EP_PgC[i]->Draw("HIST");
 
  }
  c1_uGap->Write("AC_adcUsub_Gap");
  c1_dGap->Write("AC_adcDsub_Gap");
  c1->Write("AC_adcsub");
  c1_ac->Write("ac_f");
  c3->Write("final_var");
  c4->Write("PgC ADC 12 Gaps");
  c1_pPgC->Write("PgC Mom 12 Gaps");
  c5tp->Write("tp ADC_vs_Mom");
  c5->Write("ADC_vs_Mom");
  cEp1->Write("E1_vs_p");
  c1_com->Write("accom12"); 
  c1->Write("acadcu"); 
  c2->Write("acadcd"); 
  c1_AcUD->Write("acUD2D");
  c1tdc->Write("timeU");
  c2tdc->Write("timeD");
  fclose(AC_U);
  fclose(AC_D);


  return 0;
};

// Cooker command line option function
Long_t Tof_M2::cmdline(char *cmd)
{

  return 0;
};


extern "C"
{
  Plugin *factory(TTree *in, TTree *out,TFile *inf_, TFile * outf_, TObject *p)
  {
    return (Plugin *) new Tof_M2(in,out,inf_,outf_,p);
  }
}

ClassImp(Tof_M2);
