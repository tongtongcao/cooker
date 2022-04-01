#include <Det_MWPC.h>

#include<cmath>

#include <fstream>
#include <iostream>
#include <TStyle.h>
#include <sstream>
#include "TLine.h"
#include "TCanvas.h"
#include "C2X_Strip_Trans.h"
#include "TFile.h"
#include "TTree.h"
using namespace std;

static MwpcPedestal* pedC2XL[56];
static MwpcPedestal* pedC2XR[56];
static MwpcPedestal* pedC2YL[16];
static MwpcPedestal* pedC2YR[16];
static MwpcPedestal* pedC3XL[64];
static MwpcPedestal* pedC3XR[64];
static MwpcPedestal* pedC3YL[16];
static MwpcPedestal* pedC3YR[16];
static MwpcPedestal* pedC4XL[72];
static MwpcPedestal* pedC4XR[72];
static MwpcPedestal* pedC4YL[16];
static MwpcPedestal* pedC4YR[16];

static TH1D* h1PedestalC2XL;
static TH1D* h1PedestalC2XR;
static TH1D* h1PedestalC2YL;
static TH1D* h1PedestalC2YR;
static TH1D* h1PedestalC3XL;
static TH1D* h1PedestalC3XR;
static TH1D* h1PedestalC3YL;
static TH1D* h1PedestalC3YR;
static TH1D* h1PedestalC4XL;
static TH1D* h1PedestalC4XR;
static TH1D* h1PedestalC4YL;
static TH1D* h1PedestalC4YR;

Long_t Det_MWPC::histos_ped()
{
  h1PedestalC2XL=new TH1D("PedestalC2XL","Pedestal for MWPCC2XL;Strip;Pedestal",56,-0.5,55.5);
  h1PedestalC2XR=new TH1D("PedestalC2XR","Pedestal for MWPCC2XR;Strip;Pedestal",56,-0.5,55.5);
  h1PedestalC2YL=new TH1D("PedestalC2YL","Pedestal for MWPCC2YL;Strip;Pedestal",16,-0.5,15.5);
  h1PedestalC2YR=new TH1D("PedestalC2YR","Pedestal for MWPCC2YR;Strip;Pedestal",16,-0.5,15.5);
  h1PedestalC3XL=new TH1D("PedestalC3XL","Pedestal for MWPCC3XL;Strip;Pedestal",64,-0.5,63.5);
  h1PedestalC3XR=new TH1D("PedestalC3XR","Pedestal for MWPCC3XR;Strip;Pedestal",64,-0.5,63.5);
  h1PedestalC3YL=new TH1D("PedestalC3YL","Pedestal for MWPCC3YL;Strip;Pedestal",16,-0.5,15.5);
  h1PedestalC3YR=new TH1D("PedestalC3YR","Pedestal for MWPCC3YR;Strip;Pedestal",16,-0.5,15.5);
  h1PedestalC4XL=new TH1D("PedestalC4XL","Pedestal for MWPCC4XL;Strip;Pedestal",72,-0.5,71.5);
  h1PedestalC4XR=new TH1D("PedestalC4XR","Pedestal for MWPCC4XR;Strip;Pedestal",72,-0.5,71.5);
  h1PedestalC4YL=new TH1D("PedestalC4YL","Pedestal for MWPCC4YL;Strip;Pedestal",16,-0.5,15.5);
  h1PedestalC4YR=new TH1D("PedestalC4YR","Pedestal for MWPCC4YR;Strip;Pedestal",16,-0.5,15.5);

  int binNum=100,lowerValue=0,upperValue=500;
  ostringstream name;
  for(int iChy=0;iChy<16;iChy++){
    name<<"MwpcC2YL_St"<<iChy;
    pedC2YL[iChy]=new MwpcPedestal(binNum,lowerValue,upperValue,name.str().c_str());
    name.str("");
    name<<"MwpcC3YL_St"<<iChy;
    pedC3YL[iChy]=new MwpcPedestal(binNum,lowerValue,upperValue,name.str().c_str());
    name.str("");
    name<<"MwpcC4YL_St"<<iChy;
    pedC4YL[iChy]=new MwpcPedestal(binNum,lowerValue,upperValue,name.str().c_str());
    name.str("");
  }

  for(int iChx2=0;iChx2<56;iChx2++){
    name<<"MwpcC2XL_St"<<iChx2;
    pedC2XL[iChx2]=new MwpcPedestal(binNum,lowerValue,upperValue,name.str().c_str());
    name.str("");
  }

  for(int iChx3=0;iChx3<64;iChx3++){
    name<<"MwpcC3XL_St"<<iChx3;
    pedC3XL[iChx3]=new MwpcPedestal(binNum,lowerValue,upperValue,name.str().c_str());
    name.str("");
  }

  for(int iChx4=0;iChx4<72;iChx4++){
    name<<"MwpcC4XL_St"<<iChx4;
    pedC4XL[iChx4]=new MwpcPedestal(binNum,lowerValue,upperValue,name.str().c_str());
    name.str("");
  }


  for(int iChy=0;iChy<16;iChy++){
    name<<"MwpcC2YR_St"<<iChy;
    pedC2YR[iChy]=new MwpcPedestal(binNum,lowerValue,upperValue,name.str().c_str());
    name.str("");
    name<<"MwpcC3YR_St"<<iChy;
    pedC3YR[iChy]=new MwpcPedestal(binNum,lowerValue,upperValue,name.str().c_str());
    name.str("");
    name<<"MwpcC4YR_St"<<iChy;
    pedC4YR[iChy]=new MwpcPedestal(binNum,lowerValue,upperValue,name.str().c_str());
    name.str("");
  }

  for(int iChx2=0;iChx2<56;iChx2++){
    name<<"MwpcC2XR_St"<<iChx2;
    pedC2XR[iChx2]=new MwpcPedestal(binNum,lowerValue,upperValue,name.str().c_str());
    name.str("");
  }

  for(int iChx3=0;iChx3<64;iChx3++){
    name<<"MwpcC3XR_St"<<iChx3;
    pedC3XR[iChx3]=new MwpcPedestal(binNum,lowerValue,upperValue,name.str().c_str());
    name.str("");
  }
  
  for(int iChx4=0;iChx4<72;iChx4++){
    name<<"MwpcC4XR_St"<<iChx4;
    pedC4XR[iChx4]=new MwpcPedestal(binNum,lowerValue,upperValue,name.str().c_str());
    name.str("");
  }

  return 0;

}

Long_t Det_MWPC::startup_ped()
{
  getBranchObject("RawMwpcInfo",(TObject **) &treeRaw);
  gStyle->SetOptStat(0000);
  gStyle->SetOptFit(1111);
  return 0;
};


Long_t Det_MWPC::process_ped()
{
  for(int iChy=0;iChy<16;iChy++){
    if(treeRaw->ADC_C2Y_L[iChy]>0)
      pedC2YL[iChy]->addData(treeRaw->ADC_C2Y_L[iChy]);
    if(treeRaw->ADC_C2Y_R[iChy]>0)
      pedC2YR[iChy]->addData(treeRaw->ADC_C2Y_R[iChy]);
    if(treeRaw->ADC_C3Y_L[iChy]>0)
      pedC3YL[iChy]->addData(treeRaw->ADC_C3Y_L[iChy]);
    if(treeRaw->ADC_C3Y_R[iChy]>0)
      pedC3YR[iChy]->addData(treeRaw->ADC_C3Y_R[iChy]);
    if(treeRaw->ADC_C4Y_L[iChy]>0)
      pedC4YL[iChy]->addData(treeRaw->ADC_C4Y_L[iChy]);
    if(treeRaw->ADC_C4Y_R[iChy]>0)
      pedC4YR[iChy]->addData(treeRaw->ADC_C4Y_R[iChy]);
  }
  for(int iChx2=0;iChx2<56;iChx2++){
    if(treeRaw->ADC_C2X_L[iChx2]>0)
      pedC2XL[iChx2]->addData(treeRaw->ADC_C2X_L[iChx2]); // Strip transformation
    if(treeRaw->ADC_C2X_R[iChx2]>0)
      pedC2XR[iChx2]->addData(treeRaw->ADC_C2X_R[iChx2]); // Strip transformation
  }
  for(int iChx3=0;iChx3<64;iChx3++){
    if(treeRaw->ADC_C3X_L[iChx3]>0)
      pedC3XL[iChx3]->addData(treeRaw->ADC_C3X_L[iChx3]);
    if(treeRaw->ADC_C3X_R[iChx3]>0)
      pedC3XR[iChx3]->addData(treeRaw->ADC_C3X_R[iChx3]);
  }
  for(int iChx4=0;iChx4<72;iChx4++){
    if(treeRaw->ADC_C4X_L[iChx4]>0)
      pedC4XL[iChx4]->addData(treeRaw->ADC_C4X_L[iChx4]);
    if(treeRaw->ADC_C4X_R[iChx4]>0)
      pedC4XR[iChx4]->addData(treeRaw->ADC_C4X_R[iChx4]);
  }
  return 0;
};


Long_t Det_MWPC::done_ped(){
  double pddC2YL[16],pddC2YR[16],pddC3YL[16],pddC3YR[16],pddC4YL[16],pddC4YR[16];
  double pddC2XL[56],pddC2XR[56],pddC3XL[64],pddC3XR[64],pddC4XL[72],pddC4XR[72];
  double errC2YL[16],errC2YR[16],errC3YL[16],errC3YR[16],errC4YL[16],errC4YR[16]; 
  double errC2XL[56],errC2XR[56],errC3XL[64],errC3XR[64],errC4XL[72],errC4XR[72];

  double sigmaC2YL[16],sigmaC2YR[16],sigmaC3YL[16],sigmaC3YR[16],sigmaC4YL[16],sigmaC4YR[16];
  double sigmaC2XL[56],sigmaC2XR[56],sigmaC3XL[64],sigmaC3XR[64],sigmaC4XL[72],sigmaC4XR[72];
  double sigmaErrC2YL[16],sigmaErrC2YR[16],sigmaErrC3YL[16],sigmaErrC3YR[16],sigmaErrC4YL[16],sigmaErrC4YR[16];
  double sigmaErrC2XL[56],sigmaErrC2XR[56],sigmaErrC3XL[64],sigmaErrC3XR[64],sigmaErrC4XL[72],sigmaErrC4XR[72];

  for(int iChy=0;iChy<16;iChy++){
    pddC2YL[iChy]=pedC2YL[iChy]->getPedestal();
    errC2YL[iChy]=pedC2YL[iChy]->getError();
    sigmaC2YL[iChy]=pedC2YL[iChy]->getSigma();
    sigmaErrC2YL[iChy]=pedC2YL[iChy]->getSigmaError();
    h1PedestalC2YL->SetBinContent(1+iChy,pddC2YL[iChy]);
    h1PedestalC2YL->SetBinError(1+iChy,sigmaC2YL[iChy]);
    pedC2YL[iChy]->savePlot();
    pddC3YL[iChy]=pedC3YL[iChy]->getPedestal();
    errC3YL[iChy]=pedC3YL[iChy]->getError();
    sigmaC3YL[iChy]=pedC3YL[iChy]->getSigma();
    sigmaErrC3YL[iChy]=pedC3YL[iChy]->getSigmaError();
    h1PedestalC3YL->SetBinContent(1+iChy,pddC3YL[iChy]);
    h1PedestalC3YL->SetBinError(1+iChy,sigmaC3YL[iChy]);
    pedC3YL[iChy]->savePlot();
    pddC4YL[iChy]=pedC4YL[iChy]->getPedestal();
    errC4YL[iChy]=pedC4YL[iChy]->getError();
    sigmaC4YL[iChy]=pedC4YL[iChy]->getSigma();
    sigmaErrC4YL[iChy]=pedC4YL[iChy]->getSigmaError();
    h1PedestalC4YL->SetBinContent(1+iChy,pddC4YL[iChy]);
    h1PedestalC4YL->SetBinError(1+iChy,sigmaC4YL[iChy]);
    pedC4YL[iChy]->savePlot();
    pddC2YR[iChy]=pedC2YR[iChy]->getPedestal();
    errC2YR[iChy]=pedC2YR[iChy]->getError();
    sigmaC2YR[iChy]=pedC2YR[iChy]->getSigma();
    sigmaErrC2YR[iChy]=pedC2YR[iChy]->getSigmaError();
    h1PedestalC2YR->SetBinContent(1+iChy,pddC2YR[iChy]);
    h1PedestalC2YR->SetBinError(1+iChy,sigmaC2YR[iChy]);
    pedC2YR[iChy]->savePlot();
    pddC3YR[iChy]=pedC3YR[iChy]->getPedestal();
    errC3YR[iChy]=pedC3YR[iChy]->getError();
    sigmaC3YR[iChy]=pedC3YR[iChy]->getSigma();
    sigmaErrC3YR[iChy]=pedC3YR[iChy]->getSigmaError();
    h1PedestalC3YR->SetBinContent(1+iChy,pddC3YR[iChy]);
    h1PedestalC3YR->SetBinError(1+iChy,sigmaC3YR[iChy]);
    pedC3YR[iChy]->savePlot();
    pddC4YR[iChy]=pedC4YR[iChy]->getPedestal();
    errC4YR[iChy]=pedC4YR[iChy]->getError();
    sigmaC4YR[iChy]=pedC4YR[iChy]->getSigma();
    sigmaErrC4YR[iChy]=pedC4YR[iChy]->getSigmaError();
    h1PedestalC4YR->SetBinContent(1+iChy,pddC4YR[iChy]);
    h1PedestalC4YR->SetBinError(1+iChy,sigmaC4YR[iChy]);
    pedC4YR[iChy]->savePlot();
  }

  for(int iChx2=0;iChx2<56;iChx2++){
    pddC2XL[iChx2]=pedC2XL[iChx2]->getPedestal();
    errC2XL[iChx2]=pedC2XL[iChx2]->getError();
    sigmaC2XL[iChx2]=pedC2XL[iChx2]->getSigma();
    sigmaErrC2XL[iChx2]=pedC2XL[iChx2]->getSigmaError();
    h1PedestalC2XL->SetBinContent(1+iChx2,pddC2XL[iChx2]);
    h1PedestalC2XL->SetBinError(1+iChx2,sigmaC2XL[iChx2]);
    pedC2XL[iChx2]->savePlot();
    pddC2XR[iChx2]=pedC2XR[iChx2]->getPedestal();
    errC2XR[iChx2]=pedC2XR[iChx2]->getError();
    sigmaC2XR[iChx2]=pedC2XR[iChx2]->getSigma();
    sigmaErrC2XR[iChx2]=pedC2XR[iChx2]->getSigmaError();
    h1PedestalC2XR->SetBinContent(1+iChx2,pddC2XR[iChx2]);
    h1PedestalC2XR->SetBinError(1+iChx2,sigmaC2XR[iChx2]);
    pedC2XR[iChx2]->savePlot();
  }

  for(int iChx3=0;iChx3<64;iChx3++){
    pddC3XL[iChx3]=pedC3XL[iChx3]->getPedestal();
    errC3XL[iChx3]=pedC3XL[iChx3]->getError();
    sigmaC3XL[iChx3]=pedC3XL[iChx3]->getSigma();
    sigmaErrC3XL[iChx3]=pedC3XL[iChx3]->getSigmaError();
    h1PedestalC3XL->SetBinContent(1+iChx3,pddC3XL[iChx3]);
    h1PedestalC3XL->SetBinError(1+iChx3,sigmaC3XL[iChx3]);
    pedC3XL[iChx3]->savePlot();
    pddC3XR[iChx3]=pedC3XR[iChx3]->getPedestal();
    errC3XR[iChx3]=pedC3XR[iChx3]->getError();
    sigmaC3XR[iChx3]=pedC3XR[iChx3]->getSigma();
    sigmaErrC3XR[iChx3]=pedC3XR[iChx3]->getSigmaError();
    h1PedestalC3XR->SetBinContent(1+iChx3,pddC3XR[iChx3]);
    h1PedestalC3XR->SetBinError(1+iChx3,sigmaC3XR[iChx3]);
    pedC3XR[iChx3]->savePlot();
  }

  for(int iChx4=0;iChx4<72;iChx4++){
    pddC4XL[iChx4]=pedC4XL[iChx4]->getPedestal();
    errC4XL[iChx4]=pedC4XL[iChx4]->getError();
    sigmaC4XL[iChx4]=pedC4XL[iChx4]->getSigma();
    sigmaErrC4XL[iChx4]=pedC4XL[iChx4]->getSigmaError();
    h1PedestalC4XL->SetBinContent(1+iChx4,pddC4XL[iChx4]);
    h1PedestalC4XL->SetBinError(1+iChx4,sigmaC4XL[iChx4]);
    pedC4XL[iChx4]->savePlot();
    pddC4XR[iChx4]=pedC4XR[iChx4]->getPedestal();
    errC4XR[iChx4]=pedC4XR[iChx4]->getError();
    sigmaC4XR[iChx4]=pedC4XR[iChx4]->getSigma();
    sigmaErrC4XR[iChx4]=pedC4XR[iChx4]->getSigmaError();
    h1PedestalC4XR->SetBinContent(1+iChx4,pddC4XR[iChx4]);
    h1PedestalC4XR->SetBinError(1+iChx4,sigmaC4XR[iChx4]);
    pedC4XR[iChx4]->savePlot();
  }

  TFile *file=new TFile("pedestal.root","recreate");
  TTree *tree=new TTree("pedestal","");
  tree->Branch("pddC2XL",pddC2XL,"pddC2XL[56]/D");
  tree->Branch("pddC2XR",pddC2XR,"pddC2XR[56]/D");
  tree->Branch("pddC3XL",pddC3XL,"pddC3XL[64]/D");
  tree->Branch("pddC3XR",pddC3XR,"pddC3XR[64]/D");
  tree->Branch("pddC4XL",pddC4XL,"pddC4XL[72]/D");
  tree->Branch("pddC4XR",pddC4XR,"pddC4XR[72]/D");

  tree->Branch("pddC2YL",pddC2YL,"pddC2YL[16]/D");
  tree->Branch("pddC2YR",pddC2YR,"pddC2YR[16]/D");
  tree->Branch("pddC3YL",pddC3YL,"pddC3YL[16]/D");
  tree->Branch("pddC3YR",pddC3YR,"pddC3YR[16]/D");
  tree->Branch("pddC4YL",pddC4YL,"pddC4YL[16]/D");
  tree->Branch("pddC4YR",pddC4YR,"pddC4YR[16]/D");

  tree->Branch("errC2XL",errC2XL,"errC2XL[56]/D");
  tree->Branch("errC2XR",errC2XR,"errC2XR[56]/D");
  tree->Branch("errC3XL",errC3XL,"errC3XL[64]/D");
  tree->Branch("errC3XR",errC3XR,"errC3XR[64]/D");
  tree->Branch("errC4XL",errC4XL,"errC4XL[72]/D");
  tree->Branch("errC4XR",errC4XR,"errC4XR[72]/D");

  tree->Branch("errC2YL",errC2YL,"errC2YL[16]/D");
  tree->Branch("errC2YR",errC2YR,"errC2YR[16]/D");
  tree->Branch("errC3YL",errC3YL,"errC3YL[16]/D");
  tree->Branch("errC3YR",errC3YR,"errC3YR[16]/D");
  tree->Branch("errC4YL",errC4YL,"errC4YL[16]/D");
  tree->Branch("errC4YR",errC4YR,"errC4YR[16]/D");

  tree->Branch("sigmaC2XL",sigmaC2XL,"sigmaC2XL[56]/D");
  tree->Branch("sigmaC2XR",sigmaC2XR,"sigmaC2XR[56]/D");
  tree->Branch("sigmaC3XL",sigmaC3XL,"sigmaC3XL[64]/D");
  tree->Branch("sigmaC3XR",sigmaC3XR,"sigmaC3XR[64]/D");
  tree->Branch("sigmaC4XL",sigmaC4XL,"sigmaC4XL[72]/D");
  tree->Branch("sigmaC4XR",sigmaC4XR,"sigmaC4XR[72]/D");

  tree->Branch("sigmaC2YL",sigmaC2YL,"sigmaC2YL[16]/D");
  tree->Branch("sigmaC2YR",sigmaC2YR,"sigmaC2YR[16]/D");
  tree->Branch("sigmaC3YL",sigmaC3YL,"sigmaC3YL[16]/D");
  tree->Branch("sigmaC3YR",sigmaC3YR,"sigmaC3YR[16]/D");
  tree->Branch("sigmaC4YL",sigmaC4YL,"sigmaC4YL[16]/D");
  tree->Branch("sigmaC4YR",sigmaC4YR,"sigmaC4YR[16]/D");

  tree->Branch("sigmaErrC2XL",sigmaErrC2XL,"sigmaErrC2XL[56]/D");
  tree->Branch("sigmaErrC2XR",sigmaErrC2XR,"sigmaErrC2XR[56]/D");
  tree->Branch("sigmaErrC3XL",sigmaErrC3XL,"sigmaErrC3XL[64]/D");
  tree->Branch("sigmaErrC3XR",sigmaErrC3XR,"sigmaErrC3XR[64]/D");
  tree->Branch("sigmaErrC4XL",sigmaErrC4XL,"sigmaErrC4XL[72]/D");
  tree->Branch("sigmaErrC4XR",sigmaErrC4XR,"sigmaErrC4XR[72]/D");

  tree->Branch("sigmaErrC2YL",sigmaErrC2YL,"sigmaErrC2YL[16]/D");
  tree->Branch("sigmaErrC2YR",sigmaErrC2YR,"sigmaErrC2YR[16]/D");
  tree->Branch("sigmaErrC3YL",sigmaErrC3YL,"sigmaErrC3YL[16]/D");
  tree->Branch("sigmaErrC3YR",sigmaErrC3YR,"sigmaErrC3YR[16]/D");
  tree->Branch("sigmaErrC4YL",sigmaErrC4YL,"sigmaErrC4YL[16]/D");
  tree->Branch("sigmaErrC4YR",sigmaErrC4YR,"sigmaErrC4YR[16]/D");

  tree->Fill();
  tree->Write();

  TCanvas *canv=new TCanvas("mwpcPedestalVsSt","mwpcPedestalvsSt",1200,900);
  canv->Divide(4,3);
  canv->cd(1);
  h1PedestalC2XL->SetMinimum(0);
  h1PedestalC2XL->SetMaximum(200);
  h1PedestalC2XL->Draw();
  canv->cd(2);
  h1PedestalC2XR->SetMinimum(0);
  h1PedestalC2XR->SetMaximum(200);
  h1PedestalC2XR->Draw();
  canv->cd(3);
  h1PedestalC2YL->SetMinimum(0);
  h1PedestalC2YL->SetMaximum(200);
  h1PedestalC2YL->Draw();
  canv->cd(4);
  h1PedestalC2YR->SetMinimum(0);
  h1PedestalC2YR->SetMaximum(200);
  h1PedestalC2YR->Draw();
  canv->cd(5);
  h1PedestalC3XL->SetMinimum(0);
  h1PedestalC3XL->SetMaximum(200);
  h1PedestalC3XL->Draw();
  canv->cd(6);
  h1PedestalC3XR->SetMinimum(0);
  h1PedestalC3XR->SetMaximum(200);
  h1PedestalC3XR->Draw();
  canv->cd(7);
  h1PedestalC3YL->SetMinimum(0);
  h1PedestalC3YL->SetMaximum(200);
  h1PedestalC3YL->Draw();
  canv->cd(8);
  h1PedestalC3YR->SetMinimum(0);
  h1PedestalC3YR->SetMaximum(200);
  h1PedestalC3YR->Draw();
  canv->cd(9);
  h1PedestalC4XL->SetMinimum(0);
  h1PedestalC4XL->SetMaximum(200);
  h1PedestalC4XL->Draw();
  canv->cd(10);
  h1PedestalC4XR->SetMinimum(0);
  h1PedestalC4XR->SetMaximum(200);
  h1PedestalC4XR->Draw();
  canv->cd(11);
  h1PedestalC4YL->SetMinimum(0);
  h1PedestalC4YL->SetMaximum(200);
  h1PedestalC4YL->Draw();
  canv->cd(12);
  h1PedestalC4YR->SetMinimum(0);
  h1PedestalC4YR->SetMaximum(200);
  h1PedestalC4YR->Draw();
  canv->SaveAs("mwpcPedestalVsSt.pdf");

  ofstream output("pedestal.h");

  output<<"double pddC2YL[16]={";
  for(int i=0;i<16;i++){
    if(i<15)
      output<<pddC2YL[i]<<", ";
    else
      output<<pddC2YL[i]<<"};"<<endl;
   }

  output<<"double sigmaC2YL[16]={";
  for(int i=0;i<16;i++){
    if(i<15)
      output<<sigmaC2YL[i]<<", ";
    else
      output<<sigmaC2YL[i]<<"};"<<endl;
  }

  output<<"double pddC2YR[16]={";
  for(int i=0;i<16;i++){
    if(i<15)
      output<<pddC2YR[i]<<", ";
    else
      output<<pddC2YR[i]<<"};"<<endl;
  }

  output<<"double sigmaC2YR[16]={";
  for(int i=0;i<16;i++){
    if(i<15)
      output<<sigmaC2YR[i]<<", ";
    else
      output<<sigmaC2YR[i]<<"};"<<endl;
  }

  output<<"double pddC3YL[16]={";
  for(int i=0;i<16;i++){
    if(i<15)
      output<<pddC3YL[i]<<", ";
    else
      output<<pddC3YL[i]<<"};"<<endl;
  }

  output<<"double sigmaC3YL[16]={";
  for(int i=0;i<16;i++){
    if(i<15)
      output<<sigmaC3YL[i]<<", ";
    else
      output<<sigmaC3YL[i]<<"};"<<endl;
  }

  output<<"double pddC3YR[16]={";
  for(int i=0;i<16;i++){
    if(i<15)
      output<<pddC3YR[i]<<", ";
    else
      output<<pddC3YR[i]<<"};"<<endl;
  }

  output<<"double sigmaC3YR[16]={";
  for(int i=0;i<16;i++){
    if(i<15)
      output<<sigmaC3YR[i]<<", ";
    else
      output<<sigmaC3YR[i]<<"};"<<endl;
  }

  output<<"double pddC4YL[16]={";
  for(int i=0;i<16;i++){
    if(i<15)
      output<<pddC4YL[i]<<", ";
    else
      output<<pddC4YL[i]<<"};"<<endl;
  }

  output<<"double sigmaC4YL[16]={";
  for(int i=0;i<16;i++){
    if(i<15)
      output<<sigmaC4YL[i]<<", ";
    else
      output<<sigmaC4YL[i]<<"};"<<endl;
  }

  output<<"double pddC4YR[16]={";
  for(int i=0;i<16;i++){
    if(i<15)
      output<<pddC4YR[i]<<", ";
    else
      output<<pddC4YR[i]<<"};"<<endl;
  }

  output<<"double sigmaC4YR[16]={";
  for(int i=0;i<16;i++){
    if(i<15)
      output<<sigmaC4YR[i]<<", ";
    else
      output<<sigmaC4YR[i]<<"};"<<endl;
  }

  output<<"double pddC2XL[56]={";
  for(int i=0;i<56;i++){
    if(i<55)
      output<<pddC2XL[i]<<", ";
    else
      output<<pddC2XL[i]<<"};"<<endl;
  }

  output<<"double sigmaC2XL[56]={";
  for(int i=0;i<56;i++){
    if(i<55)
      output<<sigmaC2XL[i]<<", ";
    else
      output<<sigmaC2XL[i]<<"};"<<endl;
  }

  output<<"double pddC2XR[56]={";
  for(int i=0;i<56;i++){
    if(i<55)
      output<<pddC2XR[i]<<", ";
    else
      output<<pddC2XR[i]<<"};"<<endl;
  }

  output<<"double sigmaC2XR[56]={";
  for(int i=0;i<56;i++){
    if(i<55)
      output<<sigmaC2XR[i]<<", ";
    else
      output<<sigmaC2XR[i]<<"};"<<endl;
  }

  output<<"double pddC3XL[64]={";
  for(int i=0;i<64;i++){
    if(i<63)
      output<<pddC3XL[i]<<", ";
    else
      output<<pddC3XL[i]<<"};"<<endl;
  }

  output<<"double sigmaC3XL[64]={";
  for(int i=0;i<64;i++){
    if(i<63)
      output<<sigmaC3XL[i]<<", ";
    else
      output<<sigmaC3XL[i]<<"};"<<endl;
  }

  output<<"double pddC3XR[64]={";
  for(int i=0;i<64;i++){
    if(i<63)
      output<<pddC3XR[i]<<", ";
    else
      output<<pddC3XR[i]<<"};"<<endl;
  }

  output<<"double sigmaC3XR[64]={";
  for(int i=0;i<64;i++){
    if(i<63)
      output<<sigmaC3XR[i]<<", ";
    else
      output<<sigmaC3XR[i]<<"};"<<endl;
  }

  output<<"double pddC4XL[72]={";
  for(int i=0;i<72;i++){
    if(i<71)
      output<<pddC4XL[i]<<", ";
    else
      output<<pddC4XL[i]<<"};"<<endl;
  }

  output<<"double sigmaC4XL[72]={";
  for(int i=0;i<72;i++){
    if(i<71)
    output<<sigmaC4XL[i]<<", ";
  else
    output<<sigmaC4XL[i]<<"};"<<endl;
  }

  output<<"double pddC4XR[72]={";
  for(int i=0;i<72;i++){
    if(i<71)
      output<<pddC4XR[i]<<", ";
    else
      output<<pddC4XR[i]<<"};"<<endl;
  }

  output<<"double sigmaC4XR[72]={";
  for(int i=0;i<72;i++){
    if(i<71)
      output<<sigmaC4XR[i]<<", ";
    else
      output<<sigmaC4XR[i]<<"};"<<endl;
  }

  output.close();

  return 0;
};
