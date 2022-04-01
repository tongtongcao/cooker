#include <Det_Ac.h>

#include<cmath>

#include <fstream>
#include <iostream>
#include <TStyle.h>
using namespace std;

Long_t Det_Ac::hist_calib(){
  h1Pedestal=new TH1D("Peak","Calibration Peak for AC;AC Paddle [U1-12,D1-12]",24,0.5,24.5);
  for(int iAc=0;iAc<12;iAc++){
    string name("AcU_");
    name+=(iAc+1)/10+'0';
    name+=(iAc+1)%10+'0';
    pedADC_AcU[iAc]=new peakSelect(150,0,4000,name);
  }
  for(int iAc=0;iAc<12;iAc++){
    string name("AcD_");
    name+=(iAc+1)/10+'0';
    name+=(iAc+1)%10+'0';
    pedADC_AcD[iAc]=new peakSelect(150,0,4000,name);
  }
  return 0;

}

Long_t Det_Ac::start_calib(){
  getBranchObject("RawAcInfo",(TObject **) &treeRaw);
  gStyle->SetOptStat(0);
  return 0;
};

Long_t Det_Ac::calib_peak(){
  for(int iAc = 0; iAc < 12; iAc++){
    if(treeRaw->ADC_ACU[iAc] > 0)
      pedADC_AcU[iAc]->addData(treeRaw->ADC_ACU[iAc]);
    if(treeRaw->ADC_ACD[iAc] > 0)
      pedADC_AcD[iAc]->addData(treeRaw->ADC_ACD[iAc]);
  }
  return 0;
};

Long_t Det_Ac::done_calib(){
  double pedestal[24];
  for(int iAc = 0; iAc < 12; iAc++){
    cout<<iAc+1<<"&";
    double ped1U=pedADC_AcU[iAc]->getPedestal();
    double err1U=pedADC_AcU[iAc]->getError();
    double sigmU=pedADC_AcU[iAc]->getSigma();
    double errSU=pedADC_AcU[iAc]->getSigmaError();
    h1Pedestal->SetBinContent(1+iAc,ped1U);
    h1Pedestal->SetBinError(1+iAc,err1U);
    h1Pedestal->SetBinError(1+iAc,sigmU);
    h1Pedestal->SetBinError(1+iAc,errSU);
    pedestal[iAc]=ped1U;
    cout<<ped1U<<"&";

    double ped1D=pedADC_AcD[iAc]->getPedestal();
    double err1D=pedADC_AcD[iAc]->getError();
    double sigmD=pedADC_AcU[iAc]->getSigma();
    double errSD=pedADC_AcU[iAc]->getSigmaError();
    h1Pedestal->SetBinContent(13+iAc,ped1D);
    h1Pedestal->SetBinError(13+iAc,err1D);
    h1Pedestal->SetBinError(1+iAc,sigmD);
    h1Pedestal->SetBinError(1+iAc,errSD);
    cout<<ped1D<<"\\"<<endl;
    pedestal[iAc+12]=ped1D;
  }
  cout<<"vector<double> listPedestal(24,0);"<<endl;
  for(int iAc=0;iAc<12;iAc++){
    cout<<"CalibPeak["<<iAc<<"]="<<pedestal[iAc]<<";"<<endl;
    cout<<"CalibPeak["<<iAc+12<<"]="<<pedestal[iAc+12]<<";"<<endl;
  }
  cout<<"  " <<endl;
  for(int iAc=0;iAc<12;iAc++){
    cout<<pedestal[iAc]<<"    "<<pedestal[iAc+12]<<endl;
  }

  return 0;
};
