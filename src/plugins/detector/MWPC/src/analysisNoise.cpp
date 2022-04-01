#include <Det_MWPC.h>
#include <cmath>
#include <vector>
#include <fstream>
#include <iostream>
#include <TStyle.h>
#include <sstream>
#include "TLine.h"
#include <string>
#include "TH1D.h"
#include "TCanvas.h"
#include "TText.h"
#include "TLegend.h"
#include "TLine.h"
#include "TLatex.h"
#include "TTree.h"
#include "TF1.h"
using namespace std;

static TH2D *histAdcVsStNoiseThTCCut[12];

Long_t Det_MWPC::histos_analysisNoise(){
  ostringstream sstr_name,sstr_title;
  for(int i=0;i<12;i++){
    sstr_name<<"Adc Vs Strip for noise of "<<"MwpcC"<<i/4+2<<(char)((i%4)/2+88)<<(char)((i%4%2)*6+76)<<" for Threshold and Total-Charge Cuts";
    sstr_title<<sstr_name.str()<<";Strip Number; ADC Value";
    histAdcVsStNoiseThTCCut[i]=new TH2D(sstr_name.str().c_str(),sstr_title.str().c_str(),((i%4)/2*(-5)+(i/4)*(!((i%4)/2))+7)*8,-0.5,((i%4)/2*(-5)+(i/4)*(!((i%4)/2))+7)*8-0.5,150,0,150);
    sstr_name.str("");
    sstr_title.str("");
  }

  return 0;
}

Long_t Det_MWPC::startup_analysisNoise(){
  getBranchObject("CaliMwpcInfo",(TObject **) &treeCali);
  gStyle->SetOptStat(0000);
  gStyle->SetOptFit(1111);

  return 0;
};

static int flag1,flag2;
Long_t Det_MWPC::process_analysisNoise(){
  for(UInt_t i=0;i<12;i++){
    for(UInt_t j=0;j<treeCali->adcValuePedSub[i].size();j++){
      flag2=0;
      flag1=0;
      for(UInt_t m=0;m<treeCali->thTCCutClusterStNum[i].size();m++){
        for(UInt_t n=0;n<treeCali->thTCCutClusterStNum[i][m].size();n++){
          if(j == (UInt_t)treeCali->thTCCutClusterStNum[i][m][n]){
            flag1=1;
            break;
          }
        } 
        if(flag1==1){
          flag2=1;
          break;
        }
      }
      if(flag2==0) histAdcVsStNoiseThTCCut[i]->Fill(j,treeCali->adcValuePedSub[i][j]);
    }

  }    



  return 0;
}

Long_t Det_MWPC::done_analysisNoise(){
  TCanvas *canvAdcVsStNoiseThTCCut=new TCanvas("Adc Vs Strip for noise for Threshold and Total-Charge Cuts","Adc Vs Strip for noise for Threshold and Total-Charge Cuts",2000,1500);
  canvAdcVsStNoiseThTCCut->Divide(4,3);
  for(int i=0;i<12;i++){
    canvAdcVsStNoiseThTCCut->cd(i+1);
    gPad->SetLogz();
    histAdcVsStNoiseThTCCut[i]->Draw("colz");
  }
  canvAdcVsStNoiseThTCCut->SaveAs("Adc Vs Strip for noise for Threshold and Total-Charge Cuts.pdf");

  return 0;
};

