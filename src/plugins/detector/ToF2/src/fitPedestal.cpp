#include <Det_ToF2.h>

#include<cmath>

#include <fstream>
#include <iostream>
#include <TStyle.h>
#include "TCanvas.h"
using namespace std;



Long_t Det_ToF2::histos_ped()
{
  h1Pedestal=new TH1D("Pedestal","Pedestal for TOF2;TOF2 Paddle [AO1-12,AI,BO,BI]",48,0.5,48.5);
  for(int iTof=0;iTof<12;iTof++){
    string name("Tof2AO_");
    name+=iTof/10+'0';
    name+=iTof%10+'0';
    adcTof2AO[iTof]=new TrekAdcPedestal(500,0,5000,name);
  }
  for(int iTof=0;iTof<12;iTof++){
    string name("Tof2AI_");
    name+=iTof/10+'0';
    name+=iTof%10+'0';
    adcTof2AI[iTof]=new TrekAdcPedestal(500,0,5000,name);
  }
  for(int iTof=0;iTof<12;iTof++){
    string name("Tof2BO_");
    name+=iTof/10+'0';
    name+=iTof%10+'0';
    adcTof2BO[iTof]=new TrekAdcPedestal(500,0,5000,name);
  }
  for(int iTof=0;iTof<12;iTof++){
    string name("Tof2BI_");
    name+=iTof/10+'0';
    name+=iTof%10+'0';
    adcTof2BI[iTof]=new TrekAdcPedestal(500,0,5000,name);
  }

  return 0;

}

Long_t Det_ToF2::startup_ped()
{
  getBranchObject("RawTof2Info",(TObject **) &treeRaw);
  gStyle->SetOptStat(0);
  return 0;
};


Long_t Det_ToF2::process_ped()
{
  for(int iTof=0;iTof<12;iTof++){
    if(treeRaw->ADC_TOF2AO[iTof]>0)
      adcTof2AO[iTof]->addData(treeRaw->ADC_TOF2AO[iTof]);
    if(treeRaw->ADC_TOF2AI[iTof]>0)
      adcTof2AI[iTof]->addData(treeRaw->ADC_TOF2AI[iTof]);
    if(treeRaw->ADC_TOF2BO[iTof]>0)
      adcTof2BO[iTof]->addData(treeRaw->ADC_TOF2BO[iTof]);
    if(treeRaw->ADC_TOF2BI[iTof]>0)
      adcTof2BI[iTof]->addData(treeRaw->ADC_TOF2BI[iTof]);

  }
  //  cout<<treeBeam->TDC_Trig[0][0]<<endl;
  return 0;
};


Long_t Det_ToF2::done_ped()
{

  double pedestal[48];
  TCanvas* c1=new TCanvas("AO","AO",800,600);
  c1->Divide(4,3);
  TCanvas* c2=new TCanvas("AI","AI",800,600);
  c2->Divide(4,3);
  TCanvas* c3=new TCanvas("BO","BO",800,600);
  c3->Divide(4,3);
  TCanvas* c4=new TCanvas("BI","BI",800,600);
  c4->Divide(4,3);

  for(int iTof=0;iTof<12;iTof++){
    cout<<iTof+1<<"&";
    c1->cd(iTof+1);
    gPad->SetLogy();
    gPad->SetLogx();
    double ped2AO=adcTof2AO[iTof]->getPedestal(true);
    double err2AO=adcTof2AO[iTof]->getError();
    h1Pedestal->SetBinContent(1+iTof,ped2AO);
    h1Pedestal->SetBinError(1+iTof,err2AO);
    pedestal[iTof]=ped2AO;
    cout<<ped2AO<<"&";

    c2->cd(iTof+1);
    gPad->SetLogy();
    gPad->SetLogx();
    double ped2AI=adcTof2AI[iTof]->getPedestal(true);
    double err2AI=adcTof2AI[iTof]->getError();
    h1Pedestal->SetBinContent(12+1+iTof,ped2AI);
    h1Pedestal->SetBinError(12+1+iTof,err2AI);
    pedestal[iTof+12]=ped2AI;
    cout<<ped2AI<<"&";

    c3->cd(iTof+1);
    gPad->SetLogy();
    gPad->SetLogx();
    double ped2BO=adcTof2BO[iTof]->getPedestal(true);
    double err2BO=adcTof2BO[iTof]->getError();
    h1Pedestal->SetBinContent(24+1+iTof,ped2BO);
    h1Pedestal->SetBinError(24+1+iTof,err2BO);
    pedestal[iTof+24]=ped2BO;
    cout<<ped2BO<<"&";

    c4->cd(iTof+1);
    gPad->SetLogy();
    gPad->SetLogx();
    double ped2BI=adcTof2BI[iTof]->getPedestal(true);
    double err2BI=adcTof2BI[iTof]->getError();
    h1Pedestal->SetBinContent(36+1+iTof,ped2BI);
    h1Pedestal->SetBinError(36+1+iTof,err2BI);
    pedestal[iTof+36]=ped2BI;
    cout<<ped2BI<<"\\"<<endl;

  }
  c1->Write();
  c2->Write();
  c3->Write();
  c4->Write();
  cout<<"vector<double> listPedestal(48,0);"<<endl;
  for(int iTof=0;iTof<48;iTof++){
    cout<<"listPedestal["<<iTof<<"]="<<pedestal[iTof]<<";"<<endl;
  }
  return 0;
};
