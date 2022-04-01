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
  //        
  final_PgCAdc =  new TH1F( "final_PgCAdc","final_PgCAdc",100,0,3000);
  final_beta =  new TH1F( "final_beta","final_beta",100,0,2);
  final_acAdc =  new TH1F( "final_acAdc","final_acAdc",100,-200,4500);
  final_acAdcU =  new TH1F( "final_acAdcU","final_acAdcU",100,-200,4500);
  final_acAdcD =  new TH1F( "final_acAdcD","final_acAdcD",100,-200,4500);
  final_M2 =  new TH1F( "final_M2","final_M2",100,-10000,30000);
  for(int i = 0; i < 12; i++){
    string name("ACu_ch");
    name+=(i+1)/10+'0';
    name+=(i+1)%10+'0';
    name+="_ac";
    string title("adc Ch for ac u");
    title+=(i+1)/10+'0';
    title+=(i+1)%10+'0';
    title+="Ch_ADC";
    h1adc_ACu[i] = dH1(name.c_str(),title.c_str(),200,0,4500);

    name=string("ACd_ch");
    name+=(i+1)/10+'0';
    name+=(i+1)%10+'0';
    name+="_ac";
    title=string("adc Ch for ac d");
    title+=(i+1)/10+'0';
    title+=(i+1)%10+'0';
    title+="Ch_ADC";
    h1adc_ACd[i] = dH1(name.c_str(),title.c_str(),200,0,4500);

    string name1("PgC_ch");
    name1+=(i+1)/10+'0';
    name1+=(i+1)%10+'0';
    name1+="_PgC";
    string title1("ADC");
    title1+=(i+1)/10+'0';
    title1+=(i+1)%10+'0';

    h1adc_PgC[i] = dH1(name1.c_str(),title1.c_str(),100,0,1200);
    h1adc_PgC[i]->SetTitleFont(62);
    h1adc_PgC[i]->SetTitleSize(1.2);
    h1adc_PgC[i]->GetXaxis()->SetTitle("ADC Channel");
    h1adc_PgC[i]->GetXaxis()->CenterTitle();
    h1adc_PgC[i]->GetXaxis()->SetTitleSize(0.075);
    h1adc_PgC[i]->SetTitleOffset(0.6,"X");
    h1adc_PgC[i]->SetTitleFont(62,"X");
    h1adc_PgC[i]->SetLabelFont(62,"X");
    h1adc_PgC[i]->SetLabelFont(62,"Y");
    h1adc_PgC[i]->SetTitleFont(62,"X");
    h1adc_PgC[i]->SetTitleFont(62,"Y");
    h1adc_PgC[i]->SetLabelSize(0.055,"X");
    h1adc_PgC[i]->SetLabelSize(0.055,"Y");

for(int ii=0;ii<7;ii++){
      string namep("Adc_PGC");
      namep+=(i+1)/10+'0';
      namep+=(i+1)%10+'0';
      namep+= "_";
      namep+=(ii+1)+'0';
      namep+="_pgc";

      string titlep("PgC Gap_");
      titlep+=(i+1)/10+'0';
      titlep+=(i+1)%10+'0';
      titlep+= "_";
      titlep+=(ii+1)+'0';
      titlep+="";

      ADC_PgCh[i][ii]= dH1(namep.c_str(),titlep.c_str(),100,0,1200);
      ADC_PgCh[i][ii]->SetTitleFont(62);
      ADC_PgCh[i][ii]->SetTitleSize(1.2);
      ADC_PgCh[i][ii]->GetXaxis()->SetTitle("ADC Channel");
      ADC_PgCh[i][ii]->GetXaxis()->CenterTitle();
      ADC_PgCh[i][ii]->GetXaxis()->SetTitleSize(0.075);
      ADC_PgCh[i][ii]->SetTitleOffset(0.6,"X");
      ADC_PgCh[i][ii]->SetTitleFont(62,"X");
      ADC_PgCh[i][ii]->SetLabelFont(62,"X");
      ADC_PgCh[i][ii]->SetLabelFont(62,"Y");
      ADC_PgCh[i][ii]->SetTitleFont(62,"X");
      ADC_PgCh[i][ii]->SetTitleFont(62,"Y");
      ADC_PgCh[i][ii]->SetLabelSize(0.055,"X");
      ADC_PgCh[i][ii]->SetLabelSize(0.055,"Y");

    }



  }
}

Long_t Tof_M2::startup()
{
  getBranchObject("testBranch",(TObject **) &ptr);
  getBranchObject("tofBranch",(TObject **) &ptof);
  getBranchObject("hitAc",(TObject **) &treeAc);
   getBranchObject("hitPgC",(TObject **) &treePgC);
  getBranchObject("RawPgcInfo",(TObject **) &treeRaw); 

  TString par = "ke3";

  string walk_parm =  Form("/data/trek/E36/gautam/crnRoot/cooker/inputFiles/tofM2/walkparm/%s_walkparm.dat",par.Data());
  ifstream walkpar(walk_parm.c_str());
  if (!walkpar.is_open()){
    cout<<"cannot open the file"<<endl;
    exit(0);
  }
  walkpar>>fixed>>off>>u>>d>>AI>>AO>>BI>>BO>>u0>>d0>>AIo>>AOo>>BIo>>BOo>>p_max>>p_min;
  walkpar.close();
 

  int ix =0;
  string inputfile_A=Form("/data/trek/E36/gautam/crnRoot/cooker/inputFiles/tofM2/offset/meanA_%s_1.dat",par.Data());
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

  string inputfile_B = Form("/data/trek/E36/gautam/crnRoot/cooker/inputFiles/tofM2/offset/meanB_%s_1.dat",par.Data());
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

  for (int i=0; i<12; i++)
    {
      int j=0; 
      ifstream file(Form("/data/trek/E36/gautam/crnRoot/cooker/datfiles/PgC_3994/PgC_ped3994_%d.dat", i+1), ios::in);


      string line;
      if ( file.is_open() )
	{
	  while ( ! file.eof() )
	    {
	      getline(file,line);
	      if(line.size()>0)
		{
		  istringstream vars(line);
           
		  vars>>pgc_ped[i][j];
		  cout<<i<<"   "<<j<<"  PgC: "<<endl;
		  j++;

		}
	    }
	}
      file.close();

	}

    
      gStyle->SetOptStat(0);
      return 0;
    };

  Long_t Tof_M2::process()
  {
    fTof1TimeMean =0.;
    fTof2TimeMeanA =0.;
    fTof2TimeMeanB =0.;
    F_beta =0.;
    F_M2=0.;
    acAdcf = 0;
    acAdcfU=0;  
    acAdcfD=0;  
    sum_beta_A[ptr->fgapNumTof2]=0.;
    sum_beta_B[ptr->fgapNumTof2]=0.;
    sum_M2_A[ptr->fgapNumTof2]=0.;
    sum_M2_B[ptr->fgapNumTof2]=0.;
    ADC_PgC[ptr->fgapNumTof2]=0;
    ADC_PgCf=0;

    if(ptr->nTracks==1 &&  ptr->fVertSP <300. && 
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


      if(F_beta>0){
	fTofbeta = F_beta;
	fTofM2 = F_M2;
	final_beta->Fill(fTofbeta);
	final_M2->Fill(fTofM2);
	if(fTofM2<3000){
          



	  for(int i = 0; i < 12; i++){
	    h1adc_ACu[i]->Fill(treeAc->adcACU[i]);
	    h1adc_ACd[i]->Fill(treeAc->adcACD[i]);
	    acAdcfU+= treeAc->adcACU[i]; 
	    acAdcfD+= treeAc->adcACD[i]; 
	  }
	  acAdcf = acAdcfU+acAdcfD;         

	  final_acAdcU->Fill(acAdcfU);
	  final_acAdcD->Fill(acAdcfD);
	  final_acAdc->Fill(acAdcf); 
    
      for(int i=0;i<12;i++){
      for(int ii=0;ii<7;ii++){
         
       ADC_PgCh[i][ii]->Fill(treePgC->adcPgC[i][ii]);  
       ADC_PgC[i] += treePgC->adcPgC[i][ii];
     //if(treePgC->adcPgC[i][ii]>200)   cout<<" i:  "<<i<<" ii: "<<ii<<"  AdC:  "<<treePgC->adcPgC[i][ii]<<endl;
         }
      h1adc_PgC[i]->Fill(ADC_PgC[i]);
      ADC_PgCf+= ADC_PgC[i];
    }
    final_PgCAdc->Fill(ADC_PgCf);


	}
      }
    }
    return 0;
  };

  Long_t Tof_M2::done()
  {
    TCanvas* c1 = new TCanvas("c1", " AC ADC Ped Sub. Upstream ", 100, 100, 800, 650);
    c1->Divide(3,4);
    TCanvas* c2 = new TCanvas("c2", " AC ADC Ped Sub. Downstream ", 100, 100, 800, 650);
    c2->Divide(3,4);
    TCanvas* c3 = new TCanvas("c3", " all In one final", 100, 100, 800, 650);
    c3->Divide(3,2);
    TCanvas* c4 = new TCanvas("c4", " PGC ADC 12 gaps ", 100, 100, 800, 650);
    c4->Divide(3,4);
  
    c3->cd(1);
    final_beta->Draw("HIST");
    final_beta->GetXaxis()->SetTitle("#beta");
    final_beta->GetXaxis()->CenterTitle();
    c3->cd(2);
    final_M2->Draw("HIST");
    final_M2->GetXaxis()->SetTitle("M^{2}(MeV/c^{2})^{2}");
    final_M2->GetXaxis()->CenterTitle();
    int max_M2 = final_M2->GetMaximum();
    TLine *line = new TLine(3000,0,3000,max_M2);
    line->SetLineColor(kRed);
    line->Draw();
    c3->cd(3);
    final_acAdcU->Draw("HIST");
    final_beta->GetXaxis()->SetTitle("Adc sum (Upstream Ac)");
    final_beta->GetXaxis()->CenterTitle();
    c3->cd(4);
    final_acAdcD->Draw("HIST");
    final_acAdcD->GetXaxis()->SetTitle("Adc sum (Downstream Ac)");
    final_acAdcD->GetXaxis()->CenterTitle();

    c3->cd(5);
    final_acAdc->Draw("HIST");
    final_acAdc->GetXaxis()->SetTitle("Adc sum (Ac)");
    final_acAdc->GetXaxis()->CenterTitle();
    c3->cd(6);
    final_PgCAdc->Draw("HIST");
    final_PgCAdc->GetXaxis()->SetTitle("Adc sum (PgC)");
    final_PgCAdc->GetXaxis()->CenterTitle();
    gPad->SetLogy();
TCanvas *c11[12];
    for(int i = 0; i < 12; i++){
      c1->cd(i+1);
      string htitle("AC ADC gap: ");
      htitle+=(i+1)/10+'0';
      htitle+=(i+1)%10+'0';
      htitle+="Gap";

      c1->cd(i+1);
      h1adc_ACu[i]->SetTitle(htitle.c_str());
      h1adc_ACu[i]->SetLineColor(kMagenta-7);
      h1adc_ACu[i]->GetXaxis()->SetTitle("ADC Ch");
      gPad->SetLogy();
      h1adc_ACu[i]->Draw("HIST");
      c2->cd(i+1);
      h1adc_ACd[i]->SetTitle(htitle.c_str());
      h1adc_ACd[i]->SetLineColor(kMagenta-7);
      h1adc_ACd[i]->GetXaxis()->SetTitle("ADC Ch");
      gPad->SetLogy();
      h1adc_ACd[i]->Draw("HIST");
      c4->cd(i+1);
      h1adc_PgC[i]->Draw("HIST");  
      gPad->SetLogy();

     c11[i] = new TCanvas(Form("c11[%d]",i));
     c11[i]->Divide(4,2);
     for(int ii = 0; ii < 7; ii++){
     c11[i]->cd(ii+1);
     ADC_PgCh[i][ii]->Draw("HIST");
    gPad->SetLogy();
     
     }
   c11[i]->Write(Form("PgC_ADC_%d",i));
    }
    
    c1->Write("AC_adcUsub");
    c2->Write("AC_adcDsub");
    c3->Write("final_var");
    c4->Write("PgC ADC 12 Gaps");
   


    return 0;
  };

  // Cooker command line option function
  Long_t Tof_M2::cmdline(char *cmd)
  {
    // Add command line handling flags here

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
