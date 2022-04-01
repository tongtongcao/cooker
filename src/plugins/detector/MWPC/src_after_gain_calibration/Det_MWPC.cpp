#include <Det_MWPC.h>

#include<cmath>

#include <fstream>
#include <iostream>
#include <TStyle.h>
#include <TCanvas.h>
#include <sstream>
#include "C2X_Strip_Trans.h"
using namespace std;

static TH2D* h2AdcVSSt_c2yl;
static TH2D* h2AdcVSSt_c2yr;
static TH2D* h2AdcVSSt_c3yl;
static TH2D* h2AdcVSSt_c3yr;
static TH2D* h2AdcVSSt_c4yl;
static TH2D* h2AdcVSSt_c4yr;
static TH2D* h2AdcVSSt_c2xl;
static TH2D* h2AdcVSSt_c2xr;
static TH2D* h2AdcVSSt_c3xl;
static TH2D* h2AdcVSSt_c3xr;
static TH2D* h2AdcVSSt_c4xl;
static TH2D* h2AdcVSSt_c4xr;

static MwpcAdcSpectra* AdcSpecC2XL[56];
static MwpcAdcSpectra* AdcSpecC2XR[56];
static MwpcAdcSpectra* AdcSpecC2YL[16];
static MwpcAdcSpectra* AdcSpecC2YR[16];
static MwpcAdcSpectra* AdcSpecC3XL[64];
static MwpcAdcSpectra* AdcSpecC3XR[64];
static MwpcAdcSpectra* AdcSpecC3YL[16];
static MwpcAdcSpectra* AdcSpecC3YR[16];
static MwpcAdcSpectra* AdcSpecC4XL[72];
static MwpcAdcSpectra* AdcSpecC4XR[72];
static MwpcAdcSpectra* AdcSpecC4YL[16];
static MwpcAdcSpectra* AdcSpecC4YR[16];

// Ch16 and Ch21 (starting from 0) of c4xr have problem for ADC readout 
TH1D *h1AdcSt16_c4xr;
TH1D *h1AdcSt21_c4xr;
TH1D *h1AdcSt18_c4xr;
TH1D *h1AdcSt16_c4xr_excl0;
TH1D *h1AdcSt21_c4xr_excl0;
TH1D *h1AdcSt18_c4xr_excl0;

Det_MWPC::Det_MWPC(TTree *in_,TTree *out_,TFile *inf_, TFile * outf_,TObject *p_):Plugin(in_,out_,inf_,outf_,p_)
{
  // Set defaults for various options
  treeCali=0;
};


Det_MWPC::~Det_MWPC()
{
};

Long_t Det_MWPC::histos(){
  h2AdcVSSt_c2yl=dH2("MwpcC2YL","Adc VS Strip for MwpcC2YL;Strip;ADC",16,-0.5,15.5,150,0,5000);
  h2AdcVSSt_c3yl=dH2("MwpcC3YL","Adc VS Strip for MwpcC3YL;Strip;ADC",16,-0.5,15.5,150,0,5000);
  h2AdcVSSt_c4yl=dH2("MwpcC4YL","Adc VS Strip for MwpcC4YL;Strip;ADC",16,-0.5,15.5,150,0,5000);
  h2AdcVSSt_c2yr=dH2("MwpcC2YR","Adc VS Strip for MwpcC2YR;Strip;ADC",16,-0.5,15.5,150,0,5000);
  h2AdcVSSt_c3yr=dH2("MwpcC3YR","Adc VS Strip for MwpcC3YR;Strip;ADC",16,-0.5,15.5,150,0,5000);
  h2AdcVSSt_c4yr=dH2("MwpcC4YR","Adc VS Strip for MwpcC4YR;Strip;ADC",16,-0.5,15.5,150,0,5000);
  h2AdcVSSt_c2xl=dH2("MwpcC2XL","Adc VS Strip for MwpcC2XL;Strip;ADC",56,-0.5,55.5,150,0,5000);
  h2AdcVSSt_c2xr=dH2("MwpcC2XR","Adc VS Strip for MwpcC2XR;Strip;ADC",56,-0.5,55.5,150,0,5000);
  h2AdcVSSt_c3xl=dH2("MwpcC3XL","Adc VS Strip for MwpcC3XL;Strip;ADC",64,-0.5,63.5,150,0,5000);
  h2AdcVSSt_c3xr=dH2("MwpcC3XR","Adc VS Strip for MwpcC3XR;Strip;ADC",64,-0.5,63.5,150,0,5000);
  h2AdcVSSt_c4xl=dH2("MwpcC4XL","Adc VS Strip for MwpcC4XL;Strip;ADC",72,-0.5,71.5,150,0,5000);
  h2AdcVSSt_c4xr=dH2("MwpcC4XR","Adc VS Strip for MwpcC4XR;Strip;ADC",72,-0.5,71.5,150,0,5000);

  h1AdcSt16_c4xr=dH1("MwpcC4XR_St16", "Adc Spectra for Strip 16 of MwpcC4XR; ADC; Counts",200,0,200);
  h1AdcSt21_c4xr=dH1("MwpcC4XR_St21", "Adc Spectra for Strip 16 of MwpcC4XR; ADC; Counts",200,0,200);
  h1AdcSt16_c4xr_excl0=dH1("MwpcC4XR_St16_excl0", "Adc Spectra for Strip 16 of MwpcC4XR excluding 0; ADC; Counts",199,1,200);
  h1AdcSt21_c4xr_excl0=dH1("MwpcC4XR_St21_excl0", "Adc Spectra for Strip 16 of MwpcC4XR excluding 0; ADC; Counts",199,1,200);

  int binNum=200,lowerValue=0,upperValue=4200;
  ostringstream name;
  for(int iChy=0;iChy<16;iChy++){
    name<<"MwpcC2YL_St"<<iChy;
    AdcSpecC2YL[iChy]=new MwpcAdcSpectra(binNum,lowerValue,upperValue,name.str().c_str());
    name.str("");
    name<<"MwpcC3YL_St"<<iChy;
    AdcSpecC3YL[iChy]=new MwpcAdcSpectra(binNum,lowerValue,upperValue,name.str().c_str());
    name.str("");
    name<<"MwpcC4YL_St"<<iChy;
    AdcSpecC4YL[iChy]=new MwpcAdcSpectra(binNum,lowerValue,upperValue,name.str().c_str());
    name.str("");
  }

  for(int iChx2=0;iChx2<56;iChx2++){
    name<<"MwpcC2XL_St"<<iChx2;
    AdcSpecC2XL[iChx2]=new MwpcAdcSpectra(binNum,lowerValue,upperValue,name.str().c_str());
    name.str("");
  }

  for(int iChx3=0;iChx3<64;iChx3++){
    name<<"MwpcC3XL_St"<<iChx3;
    AdcSpecC3XL[iChx3]=new MwpcAdcSpectra(binNum,lowerValue,upperValue,name.str().c_str());
    name.str("");
  }

  for(int iChx4=0;iChx4<72;iChx4++){
    name<<"MwpcC4XL_St"<<iChx4;
    AdcSpecC4XL[iChx4]=new MwpcAdcSpectra(binNum,lowerValue,upperValue,name.str().c_str());
    name.str("");
  }


  for(int iChy=0;iChy<16;iChy++){
    name<<"MwpcC2YR_St"<<iChy;
    AdcSpecC2YR[iChy]=new MwpcAdcSpectra(binNum,lowerValue,upperValue,name.str().c_str());
    name.str("");
    name<<"MwpcC3YR_St"<<iChy;
    AdcSpecC3YR[iChy]=new MwpcAdcSpectra(binNum,lowerValue,upperValue,name.str().c_str());
    name.str("");
    name<<"MwpcC4YR_St"<<iChy;
    AdcSpecC4YR[iChy]=new MwpcAdcSpectra(binNum,lowerValue,upperValue,name.str().c_str());
    name.str("");
  }

  for(int iChx2=0;iChx2<56;iChx2++){
    name<<"MwpcC2XR_St"<<iChx2;
    AdcSpecC2XR[iChx2]=new MwpcAdcSpectra(binNum,lowerValue,upperValue,name.str().c_str());
    name.str("");
  }

  for(int iChx3=0;iChx3<64;iChx3++){
    name<<"MwpcC3XR_St"<<iChx3;
    AdcSpecC3XR[iChx3]=new MwpcAdcSpectra(binNum,lowerValue,upperValue,name.str().c_str());
    name.str("");
  }
  
  for(int iChx4=0;iChx4<72;iChx4++){
    name<<"MwpcC4XR_St"<<iChx4;
    AdcSpecC4XR[iChx4]=new MwpcAdcSpectra(binNum,lowerValue,upperValue,name.str().c_str());
    name.str("");
  }

  return 0;
}

Long_t Det_MWPC::startup()
{
  // return 0;
  getBranchObject("RawMwpcInfo",(TObject **) &treeRaw);
  getBranchObject("RawBeamInfo",(TObject **) &treeBeam);
  gStyle->SetOptStat(0);
  // cout<<treeBeam->TDC_Trig[0][0]<<endl;

  return 0;
};


Long_t Det_MWPC::process()
 
{
  // int cylmax, cylmaxbin, cylbinl, cylbinr 
  for(int i=0;i<16;i++){
    h2AdcVSSt_c2yl->Fill(i,treeRaw->ADC_C2Y_L[i]);
    h2AdcVSSt_c3yl->Fill(i,treeRaw->ADC_C3Y_L[i]);
    h2AdcVSSt_c4yl->Fill(i,treeRaw->ADC_C4Y_L[i]);
    h2AdcVSSt_c2yr->Fill(i,treeRaw->ADC_C2Y_R[i]);
    h2AdcVSSt_c3yr->Fill(i,treeRaw->ADC_C3Y_R[i]);
    h2AdcVSSt_c4yr->Fill(i,treeRaw->ADC_C4Y_R[i]);
  }

  for(int iChx2=0;iChx2<56;iChx2++){
    h2AdcVSSt_c2xl->Fill(iChx2,treeRaw->ADC_C2X_L[C2XStTrans[iChx2]]);
    h2AdcVSSt_c2xr->Fill(iChx2,treeRaw->ADC_C2X_R[C2XStTrans[iChx2]]);
  }
  for(int iChx3=0;iChx3<64;iChx3++){
    h2AdcVSSt_c3xl->Fill(iChx3,treeRaw->ADC_C3X_L[iChx3]);
    h2AdcVSSt_c3xr->Fill(iChx3,treeRaw->ADC_C3X_R[iChx3]);
  }
  for(int iChx4=0;iChx4<72;iChx4++){
    h2AdcVSSt_c4xl->Fill(iChx4,treeRaw->ADC_C4X_L[iChx4]);
    h2AdcVSSt_c4xr->Fill(iChx4,treeRaw->ADC_C4X_R[iChx4]);
  }

  h1AdcSt16_c4xr->Fill(treeRaw->ADC_C4X_L[16]);
  h1AdcSt21_c4xr->Fill(treeRaw->ADC_C4X_L[21]);
  h1AdcSt16_c4xr_excl0->Fill(treeRaw->ADC_C4X_L[16]);
  h1AdcSt21_c4xr_excl0->Fill(treeRaw->ADC_C4X_L[21]);

  for(int iChy=0;iChy<16;iChy++){
    if(treeRaw->ADC_C2Y_L[iChy]>0)
      AdcSpecC2YL[iChy]->addData(treeRaw->ADC_C2Y_L[iChy]);
    if(treeRaw->ADC_C2Y_R[iChy]>0)
      AdcSpecC2YR[iChy]->addData(treeRaw->ADC_C2Y_R[iChy]);
    if(treeRaw->ADC_C3Y_L[iChy]>0)
      AdcSpecC3YL[iChy]->addData(treeRaw->ADC_C3Y_L[iChy]);
    if(treeRaw->ADC_C3Y_R[iChy]>0)
      AdcSpecC3YR[iChy]->addData(treeRaw->ADC_C3Y_R[iChy]);
    if(treeRaw->ADC_C4Y_L[iChy]>0)
      AdcSpecC4YL[iChy]->addData(treeRaw->ADC_C4Y_L[iChy]);
    if(treeRaw->ADC_C4Y_R[iChy]>0)
      AdcSpecC4YR[iChy]->addData(treeRaw->ADC_C4Y_R[iChy]);
  }
  for(int iChx2=0;iChx2<56;iChx2++){
    if(treeRaw->ADC_C2X_L[C2XStTrans[iChx2]]>0)
      AdcSpecC2XL[iChx2]->addData(treeRaw->ADC_C2X_L[C2XStTrans[iChx2]]); // Strip transformation
    if(treeRaw->ADC_C2X_R[C2XStTrans[iChx2]]>0)
      AdcSpecC2XR[iChx2]->addData(treeRaw->ADC_C2X_R[C2XStTrans[iChx2]]); // Strip transformation
  }
  for(int iChx3=0;iChx3<64;iChx3++){
    if(treeRaw->ADC_C3X_L[iChx3]>0)
      AdcSpecC3XL[iChx3]->addData(treeRaw->ADC_C3X_L[iChx3]);
    if(treeRaw->ADC_C3X_R[iChx3]>0)
      AdcSpecC3XR[iChx3]->addData(treeRaw->ADC_C3X_R[iChx3]);
  }
  for(int iChx4=0;iChx4<72;iChx4++){
    if(treeRaw->ADC_C4X_L[iChx4]>0)
      AdcSpecC4XL[iChx4]->addData(treeRaw->ADC_C4X_L[iChx4]);
    if(treeRaw->ADC_C4X_R[iChx4]>0)
      AdcSpecC4XR[iChx4]->addData(treeRaw->ADC_C4X_R[iChx4]);
  }

  return 0;
};


Long_t Det_MWPC::done()
{
  TCanvas *canv=new TCanvas("mwpcADCVsSt","mwpcADCvsSt",1200,900);
  canv->Divide(4,3);
  canv->cd(1);
  gPad->SetLogz();
  h2AdcVSSt_c2xl->Draw("colz");
  canv->cd(2);
  gPad->SetLogz();
  h2AdcVSSt_c2xr->Draw("colz");
  canv->cd(3);
  gPad->SetLogz();
  h2AdcVSSt_c2yl->Draw("colz");
  canv->cd(4);
  gPad->SetLogz();
  h2AdcVSSt_c2yr->Draw("colz");
  canv->cd(5);
  gPad->SetLogz();
  h2AdcVSSt_c3xl->Draw("colz");
  canv->cd(6);
  gPad->SetLogz();
  h2AdcVSSt_c3xr->Draw("colz");
  canv->cd(7);
  gPad->SetLogz();
  h2AdcVSSt_c3yl->Draw("colz");
  canv->cd(8);
  gPad->SetLogz();
  h2AdcVSSt_c3yr->Draw("colz");
  canv->cd(9);
  gPad->SetLogz();
  h2AdcVSSt_c4xl->Draw("colz");
  canv->cd(10);
  gPad->SetLogz();
  h2AdcVSSt_c4xr->Draw("colz");
  canv->cd(11);
  gPad->SetLogz();
  h2AdcVSSt_c4yl->Draw("colz");
  canv->cd(12);
  gPad->SetLogz();
  h2AdcVSSt_c4yr->Draw("colz");
  canv->SaveAs("mwpcADCVsSt.pdf");

  TCanvas *canv_ch=new TCanvas("mwpcADCC4XR_St16_St21","mwpcADCC4XR_St16_St21",1000,600);
  canv_ch->Divide(2,2);
  canv_ch->cd(1);
  h1AdcSt16_c4xr->Draw();
  canv_ch->cd(2);
  h1AdcSt21_c4xr->Draw();
  canv_ch->cd(3);
  h1AdcSt16_c4xr_excl0->Draw();
  canv_ch->cd(4);
  h1AdcSt21_c4xr_excl0->Draw();
  canv_ch->SaveAs("mwpcADCC4XR_St16_St21.pdf");

  for(int iChy=0;iChy<16;iChy++){
    AdcSpecC2YL[iChy]->savePlot();
    AdcSpecC3YL[iChy]->savePlot();
    AdcSpecC4YL[iChy]->savePlot();
    AdcSpecC2YR[iChy]->savePlot();
    AdcSpecC3YR[iChy]->savePlot();
    AdcSpecC4YR[iChy]->savePlot();
  }

  for(int iChx2=0;iChx2<56;iChx2++){
    AdcSpecC2XL[iChx2]->savePlot();
    AdcSpecC2XR[iChx2]->savePlot();
  }

  for(int iChx3=0;iChx3<64;iChx3++){
    AdcSpecC3XL[iChx3]->savePlot();
    AdcSpecC3XR[iChx3]->savePlot();
  }

  for(int iChx4=0;iChx4<72;iChx4++){
    AdcSpecC4XL[iChx4]->savePlot();
    AdcSpecC4XR[iChx4]->savePlot();
  }

  return 0;
};

Long_t Det_MWPC::cmdline(char *cmd)
{
  //add cmdline hanling here

  return 0; // 0 = all ok
};

/////////////////////////////////////////////////////////////////////////////////////////////////// 

extern "C"{
  Plugin *factory(TTree *in_, TTree *out_, TFile *inf_, TFile * outf_,TObject *p_)
  {
      return (Plugin *) new Det_MWPC(in_,out_,inf_,outf_,p_);
  }
}


ClassImp(Det_MWPC);
