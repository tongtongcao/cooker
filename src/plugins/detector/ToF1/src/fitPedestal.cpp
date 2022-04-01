#include <Det_ToF1.h>

#include<cmath>

#include <fstream>
#include <iostream>
#include <TStyle.h>
using namespace std;



Long_t Det_ToF1::histos_ped()
{
  h1Pedestal=new TH1D("Pedestal","Pedestal for TOF1;TOF1 Paddle [U1-12,D1-12]",24,0.5,24.5);
  for(int iTof=0;iTof<12;iTof++){
    string name("Tof1U_");
    name+=iTof/10+'0';
    name+=iTof%10+'0';
    adcTof1U[iTof]=new TrekAdcPedestal(100,0,300,name);
  }
  for(int iTof=0;iTof<12;iTof++){
    string name("Tof1D_");
    name+=iTof/10+'0';
    name+=iTof%10+'0';
    adcTof1D[iTof]=new TrekAdcPedestal(100,0,300,name);
  }
  return 0;

}

Long_t Det_ToF1::startup_ped()
{
  getBranchObject("RawTof1Info",(TObject **) &treeRaw);
  gStyle->SetOptStat(0);
  return 0;
};


Long_t Det_ToF1::process_ped()
{
  for(int iTof=0;iTof<12;iTof++){
    if(treeRaw->ADC_TOF1U[iTof]>=0)
      adcTof1U[iTof]->addData(treeRaw->ADC_TOF1U[iTof]);
    if(treeRaw->ADC_TOF1D[iTof]>=0)
      adcTof1D[iTof]->addData(treeRaw->ADC_TOF1D[iTof]);
  }
  //  cout<<treeBeam->TDC_Trig[0][0]<<endl;
  return 0;
};


Long_t Det_ToF1::done_ped()
{

  double pedestal[24];
  for(int iTof=0;iTof<12;iTof++){
    cout<<iTof+1<<"&";
    double ped1U=adcTof1U[iTof]->getPedestal();
    double err1U=adcTof1U[iTof]->getError();
    h1Pedestal->SetBinContent(1+iTof,ped1U);
    h1Pedestal->SetBinError(1+iTof,err1U);
    pedestal[iTof]=ped1U;
    cout<<ped1U<<"&";

    double ped1D=adcTof1D[iTof]->getPedestal();
    double err1D=adcTof1D[iTof]->getError();
    h1Pedestal->SetBinContent(13+iTof,ped1D);
    h1Pedestal->SetBinError(13+iTof,err1D);
    cout<<ped1D<<"\\"<<endl;
    pedestal[iTof+12]=ped1D;
  }
  cout<<"vector<double> listPedestal(24,0);"<<endl;
  for(int iTof=0;iTof<12;iTof++){
    cout<<"listPedestal["<<iTof<<"]="<<pedestal[iTof]<<";"<<endl;
    cout<<"listPedestal["<<iTof+12<<"]="<<pedestal[iTof+12]<<";"<<endl;
  }
  return 0;
};
