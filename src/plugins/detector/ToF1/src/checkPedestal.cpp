#include <Det_ToF1.h>

#include<cmath>

#include <fstream>
#include <iostream>
#include <TStyle.h>
#include <TCanvas.h>
#include <TLine.h>
#include <TPad.h>
using namespace std;



Long_t Det_ToF1::histos_checkPed()
{
  for(int iTof=0;iTof<12;iTof++){
    string name("Tof1U_");
    name+=iTof/10+'0';
    name+=iTof%10+'0';
    h1Adc[iTof]=new TH1D(name.c_str(),"ADC Spectrum",2000,0,5000);
  }
  for(int iTof=0;iTof<12;iTof++){
    string name("Tof1D_");
    name+=iTof/10+'0';
    name+=iTof%10+'0';
    h1Adc[iTof+12]=new TH1D(name.c_str(),"ADC Spectrum",2000,0,5000);
  }
  return 0;

}

Long_t Det_ToF1::startup_checkPed()
{
  getBranchObject("RawTof1Info",(TObject **) &treeRaw);
  gStyle->SetOptStat(0);
  return 0;
};


Long_t Det_ToF1::process_checkPed()
{
  for(int iTof=0;iTof<12;iTof++){
    if(treeRaw->ADC_TOF1U[iTof]>=0)
      h1Adc[iTof]->Fill(treeRaw->ADC_TOF1U[iTof]);
    if(treeRaw->ADC_TOF1D[iTof]>=0)
      h1Adc[iTof+12]->Fill(treeRaw->ADC_TOF1D[iTof]);
  }

  return 0;
};


Long_t Det_ToF1::done_checkPed()
{
  vector<double> listPedestal(24,0);
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

  return 0;
};
