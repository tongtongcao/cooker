#include <Det_ToF2.h>

#include<cmath>

#include <fstream>
#include <iostream>
#include <TStyle.h>
#include <TCanvas.h>
#include <TLine.h>
#include <TPad.h>
using namespace std;



Long_t Det_ToF2::histos_checkPed()
{
  for(int iTof=0;iTof<12;iTof++){
    string name("Tof2AO_");
    name+=iTof/10+'0';
    name+=iTof%10+'0';
    h1Adc[iTof]=new TH1D(name.c_str(),"ADC Spectrum",2000,0,5000);
  }
  for(int iTof=0;iTof<12;iTof++){
    string name("Tof2AI_");
    name+=iTof/10+'0';
    name+=iTof%10+'0';
    h1Adc[iTof+12]=new TH1D(name.c_str(),"ADC Spectrum",2000,0,5000);
  }
  for(int iTof=0;iTof<12;iTof++){
    string name("Tof2BO_");
    name+=iTof/10+'0';
    name+=iTof%10+'0';
    h1Adc[iTof+24]=new TH1D(name.c_str(),"ADC Spectrum",2000,0,5000);
  }
  for(int iTof=0;iTof<12;iTof++){
    string name("Tof2BI_");
    name+=iTof/10+'0';
    name+=iTof%10+'0';
    h1Adc[iTof+36]=new TH1D(name.c_str(),"ADC Spectrum",2000,0,5000);
  }

  return 0;

}

Long_t Det_ToF2::startup_checkPed()
{
  getBranchObject("RawTof2Info",(TObject **) &treeRaw);
  gStyle->SetOptStat(0);
  return 0;
};


Long_t Det_ToF2::process_checkPed()
{
  for(int iTof=0;iTof<12;iTof++){
    if(treeRaw->ADC_TOF2AO[iTof]>0)
      h1Adc[iTof]->Fill(treeRaw->ADC_TOF2AO[iTof]);
    if(treeRaw->ADC_TOF2AI[iTof]>0)
      h1Adc[iTof+12]->Fill(treeRaw->ADC_TOF2AI[iTof]);
    if(treeRaw->ADC_TOF2BO[iTof]>0)
      h1Adc[iTof]->Fill(treeRaw->ADC_TOF2BO[iTof]);
    if(treeRaw->ADC_TOF2BI[iTof]>0)
      h1Adc[iTof+12]->Fill(treeRaw->ADC_TOF2BI[iTof]);
  }

  return 0;
};


Long_t Det_ToF2::done_checkPed()
{
  vector<double> listPedestal(48,0);
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
  /*
  TCanvas* c1=new TCanvas("TOF1U","TOF1U",3508,2480);
  TCanvas* c2=new TCanvas("TOF1D","TOF1D",3508,2480);
  c1->Divide(3,4);
  c2->Divide(3,4);
  for(int iTof=0;iTof<12;iTof++){
    c1->cd(iTof+1);
    h1Adc[iTof]->Draw();
    double ymax=h1Adc[iTof]->GetMaximum();
    TLine* line=new TLine(listPedestal[iTof],0,listPedestal[iTof],ymax);
    line->SetLineColor(kRed);
    line->Draw();
    gPad->SetLogx();
    gPad->SetLogy();
    c2->cd(iTof+1);
    h1Adc[iTof+12]->Draw();
    ymax=h1Adc[iTof+12]->GetMaximum();
    line=new TLine(listPedestal[iTof+12],0,listPedestal[iTof+12],ymax);
    line->SetLineColor(kRed);
    line->Draw();
    gPad->SetLogx();
    gPad->SetLogy();
    
  }
  c1->Write();
  c2->Write();
  c1->SaveAs("checkPed1U.png");
  c2->SaveAs("checkPed1D.png");
  */
  return 0;
};
