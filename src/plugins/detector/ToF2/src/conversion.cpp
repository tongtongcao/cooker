#include <Det_ToF2.h>

#include<cmath>

#include <fstream>
#include <iostream>
#include <TStyle.h>
#include <TCanvas.h>
#include <TLine.h>
#include <TPad.h>
using namespace std;



Long_t Det_ToF2::histos_conversion()
{
  h1Multi=new TH1D("Mulplicity_TOF2","Number of valid gaps",13,-0.5,12.5);
  h1TDC=new TH1D("TDC_TOF2","tdc",100,0.,100);
  h1ADC=new TH1D("ADC_TOF2","adc",4097,-0.5,4096.5);

  return 0;

}

Long_t Det_ToF2::startup_conversion()
{
  getBranchObject("RawTof2Info",(TObject **) &treeRaw);
  gStyle->SetOptStat(0);
  treeCali=new CRTCaliTof2();
  makeBranch("hitTof2",(TObject **)&treeCali);

  return 0;
};


Long_t Det_ToF2::process_conversion()
{
  static bool first=true;
  static vector<double> listPedestal(48,0);
  if(first){
    listPedestal[0]=60.1263;
    listPedestal[1]=33.8924;
    listPedestal[2]=148.012;
    listPedestal[3]=158.347;
    listPedestal[4]=186.154;
    listPedestal[5]=112.004;
    listPedestal[6]=176.498;
    listPedestal[7]=139.647;
    listPedestal[8]=38.4092;
    listPedestal[9]=39.1898;
    listPedestal[10]=77.4019;
    listPedestal[11]=48.8611;
    listPedestal[12]=50.8747;
    listPedestal[13]=56.8653;
    listPedestal[14]=58.9634;
    listPedestal[15]=71.2165;
    listPedestal[16]=60.5193;
    listPedestal[17]=47.1619;
    listPedestal[18]=40.8477;
    listPedestal[19]=70.9978;
    listPedestal[20]=62.4487;
    listPedestal[21]=41.1895;
    listPedestal[22]=49.0886;
    listPedestal[23]=41.6453;
    listPedestal[24]=53.0709;
    listPedestal[25]=44.9932;
    listPedestal[26]=55.0854;
    listPedestal[27]=44.9972;
    listPedestal[28]=62.1591;
    listPedestal[29]=54.9487;
    listPedestal[30]=75.0317;
    listPedestal[31]=67.0182;
    listPedestal[32]=54.7372;
    listPedestal[33]=71.1766;
    listPedestal[34]=48.9068;
    listPedestal[35]=59.0462;
    listPedestal[36]=59.0347;
    listPedestal[37]=47.5102;
    listPedestal[38]=56.9434;
    listPedestal[39]=47.29;
    listPedestal[40]=67.0051;
    listPedestal[41]=63.0649;
    listPedestal[42]=84.3897;
    listPedestal[43]=49.0058;
    listPedestal[44]=69.017;
    listPedestal[45]=49.5384;
    listPedestal[46]=25.2035;
    listPedestal[47]=50.3541;
    first=false;
  }
  treeCali->run=treeRaw->run;
  treeCali->event=treeRaw->event;
  treeCali->gap.clear();
  treeCali->adcAI.clear();
  treeCali->adcAO.clear();
  treeCali->adcBI.clear();
  treeCali->adcBO.clear();
  treeCali->timeAI.clear();
  treeCali->timeAO.clear();
  treeCali->timeBI.clear();
  treeCali->timeBO.clear();
  UInt_t nGap=0;
  for(int iTof=0;iTof<12;iTof++){
    if( (treeRaw->TDC_TOF2AI[iTof]>0 && treeRaw->TDC_TOF2AO[iTof]>0)||(treeRaw->TDC_TOF2BI[iTof]>0 && treeRaw->TDC_TOF2BO[iTof]>0)){
      if(treeRaw->ADC_TOF2AI[iTof]<0.5) treeRaw->ADC_TOF2AI[iTof]=4096;
      if(treeRaw->ADC_TOF2AO[iTof]<0.5) treeRaw->ADC_TOF2AO[iTof]=4096;
      if(treeRaw->ADC_TOF2BI[iTof]<0.5) treeRaw->ADC_TOF2BI[iTof]=4096;
      if(treeRaw->ADC_TOF2BO[iTof]<0.5) treeRaw->ADC_TOF2BO[iTof]=4096;
      Float_t myAdcAO=treeRaw->ADC_TOF2AO[iTof]-listPedestal[iTof]; 
      Float_t myAdcAI=treeRaw->ADC_TOF2AI[iTof]-listPedestal[iTof+12]; 
      Float_t myAdcBO=treeRaw->ADC_TOF2BO[iTof]-listPedestal[iTof+24]; 
      Float_t myAdcBI=treeRaw->ADC_TOF2BI[iTof]-listPedestal[iTof+36]; 

      h1ADC->Fill(myAdcAO);
      h1ADC->Fill(myAdcAI);
      h1ADC->Fill(myAdcBO);
      h1ADC->Fill(myAdcBI);

      h1TDC->Fill(treeRaw->TDC_TOF2AI[iTof]*0.025);
      h1TDC->Fill(treeRaw->TDC_TOF2AO[iTof]*0.025);
      h1TDC->Fill(treeRaw->TDC_TOF2BI[iTof]*0.025);
      h1TDC->Fill(treeRaw->TDC_TOF2BO[iTof]*0.025);
      nGap++;
      treeCali->gap.push_back(iTof+1);
      treeCali->adcAI.push_back(myAdcAI);
      treeCali->adcAO.push_back(myAdcAO);
      treeCali->adcBI.push_back(myAdcBI);
      treeCali->adcBO.push_back(myAdcBO);

      treeCali->timeAI.push_back(treeRaw->TDC_TOF2AI[iTof]*0.025);
      treeCali->timeAO.push_back(treeRaw->TDC_TOF2AO[iTof]*0.025);
      treeCali->timeBI.push_back(treeRaw->TDC_TOF2BI[iTof]*0.025);
      treeCali->timeBO.push_back(treeRaw->TDC_TOF2BO[iTof]*0.025);

      
    }

  }
  h1Multi->Fill(nGap); 
  
  return 0;
};


Long_t Det_ToF2::done_conversion()
{
  return 0;
};
