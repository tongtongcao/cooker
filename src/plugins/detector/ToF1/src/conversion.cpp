#include <Det_ToF1.h>

#include<cmath>

#include <fstream>
#include <iostream>
#include <TStyle.h>
using namespace std;



Long_t Det_ToF1::histos_conversion()
{
  h1Multi=new TH1D("Mulplicity_TOF1","Number of gaps having both U and D ADCs, TDCs",13,-0.5,12.5);
  h1TDC=new TH1D("TDC_TOF1","tdc",4097,-0.5,4096.5);
  h1ADC=new TH1D("ADC_TOF1","adc",4097,-0.5,4096.5);
  return 0;

}

Long_t Det_ToF1::startup_conversion()
{
  getBranchObject("RawTof1Info",(TObject **) &treeRaw);
  treeCali=new CRTCaliTof1();

  // Make new branch on output tree                                                            
  makeBranch("hitTof1",(TObject **)&treeCali);

  gStyle->SetOptStat(0);

  return 0;
};


Long_t Det_ToF1::process_conversion()
{

  static bool first=true;
  static vector<double> listPedestal(24,0);
  if(first){
    listPedestal[0]=45.0186;
    listPedestal[12]=60.0226;
    listPedestal[1]=67.484;
    listPedestal[13]=48.4007;
    listPedestal[2]=56.1659;
    listPedestal[14]=79.1373;
    listPedestal[3]=37.9318;
    listPedestal[15]=44.1452;
    listPedestal[4]=59.794;
    listPedestal[16]=62.797;
    listPedestal[5]=60.4463;
    listPedestal[17]=72.4111;
    listPedestal[6]=52.0727;
    listPedestal[18]=82.3057;
    listPedestal[7]=66.2435;
    listPedestal[19]=65.7156;
    listPedestal[8]=63.2392;
    listPedestal[20]=69.8056;
    listPedestal[9]=51.6298;
    listPedestal[21]=64.618;
    listPedestal[10]=67.5943;
    listPedestal[22]=56.4271;
    listPedestal[11]=49.2037;
    listPedestal[23]=63.7336;
    first=false;
  }
  treeCali->run=treeRaw->run;
  treeCali->event=treeRaw->event;
  treeCali->gap.clear();
  treeCali->adcU.clear();
  treeCali->adcD.clear();
  treeCali->timeU.clear();
  treeCali->timeD.clear();
  UInt_t nGap=0;
  for(int iTof=0;iTof<12;iTof++){
    if(treeRaw->TDC_TOF1U[iTof]>0 || treeRaw->TDC_TOF1D[iTof]>0){

      if(treeRaw->ADC_TOF1U[iTof]<0.5) treeRaw->ADC_TOF1U[iTof]=4096;
      
      if(treeRaw->ADC_TOF1D[iTof]<0.5) treeRaw->ADC_TOF1D[iTof]=4096;
      Float_t myAdcU=treeRaw->ADC_TOF1U[iTof]-listPedestal[iTof];
      Float_t myAdcD=treeRaw->ADC_TOF1D[iTof]-listPedestal[iTof+12];
      if(myAdcU<10 && myAdcD<10 && treeRaw->TDC_TOF1U[iTof]<0 && treeRaw->TDC_TOF1D[iTof]<0) continue;
      h1TDC->Fill(treeRaw->TDC_TOF1U[iTof]);
      h1ADC->Fill(myAdcU);
      h1ADC->Fill(myAdcD);
      h1TDC->Fill(treeRaw->TDC_TOF1D[iTof]);
      nGap++;
      treeCali->gap.push_back(iTof+1);
      treeCali->adcU.push_back(myAdcU);
      treeCali->adcD.push_back(myAdcD);
      treeCali->timeU.push_back(treeRaw->TDC_TOF1U[iTof]*0.025);
      treeCali->timeD.push_back(treeRaw->TDC_TOF1D[iTof]*0.025);
      
    }
  }
  h1Multi->Fill(nGap); 
  //  cout<<treeBeam->TDC_Trig[0][0]<<endl;
  return 0;
};


Long_t Det_ToF1::done_conversion()
{

  return 0;
};
