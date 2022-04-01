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
static TH2D *histAdcVsStNoCut[12];
static TH2D *histAdcVsStThTCCut[12];
static TH2D *histAdcVsStThTCFFCut[12];
static TH2D *histAdcVsStThTCFFPairCut[12];
static TH2D *histAdcVsStThTCFFPairSFCut[12];
static TH1D *histPosNoCut[12];
static TH1D *histPosThTCCut[12];
static TH1D *histPosThTCFFCut[12];
static TH1D *histPosThTCFFPairCut[12];
static TH1D *histPosThTCFFPairSFCut[12];
Long_t Det_MWPC::histos_analysisExtractedCluster2(){
  ostringstream sstr_name,sstr_title;
  for(int i=0;i<12;i++){
    sstr_name<<"Adc Vs Strip for all clusters of "<<"MwpcC"<<i/4+2<<(char)((i%4)/2+88)<<(char)((i%4%2)*6+76)<<" for No Cut";
    sstr_title<<sstr_name.str()<<";Strip Number; ADC Value";
    histAdcVsStNoCut[i]=new TH2D(sstr_name.str().c_str(),sstr_title.str().c_str(),((i%4)/2*(-5)+(i/4)*(!((i%4)/2))+7)*8,-0.5,((i%4)/2*(-5)+(i/4)*(!((i%4)/2))+7)*8-0.5,150,0,5000);
    sstr_name.str("");
    sstr_title.str("");

    sstr_name<<"Adc Vs Strip for all clusters of "<<"MwpcC"<<i/4+2<<(char)((i%4)/2+88)<<(char)((i%4%2)*6+76)<<" for Threshold and Total-Charge Cuts";
    sstr_title<<sstr_name.str()<<";Strip Number; ADC Value";
    histAdcVsStThTCCut[i]=new TH2D(sstr_name.str().c_str(),sstr_title.str().c_str(),((i%4)/2*(-5)+(i/4)*(!((i%4)/2))+7)*8,-0.5,((i%4)/2*(-5)+(i/4)*(!((i%4)/2))+7)*8-0.5,150,0,5000);
    sstr_name.str("");
    sstr_title.str("");  

    sstr_name<<"Adc Vs Strip for all clusters of "<<"MwpcC"<<i/4+2<<(char)((i%4)/2+88)<<(char)((i%4%2)*6+76)<<" for Threshold, Total-Charge and First-Flag Cuts";
    sstr_title<<sstr_name.str()<<";Strip Number; ADC Value";
    histAdcVsStThTCFFCut[i]=new TH2D(sstr_name.str().c_str(),sstr_title.str().c_str(),((i%4)/2*(-5)+(i/4)*(!((i%4)/2))+7)*8,-0.5,((i%4)/2*(-5)+(i/4)*(!((i%4)/2))+7)*8-0.5,150,0,5000);
    sstr_name.str("");
    sstr_title.str("");

    sstr_name<<"Adc Vs Strip for all pairs of "<<"MwpcC"<<i/4+2<<(char)((i%4)/2+88)<<(char)((i%4%2)*6+76)<<" for Threshold, Total-Charge, First-Flag and Pair Cuts";
    sstr_title<<sstr_name.str()<<";Strip Number; ADC Value";
    histAdcVsStThTCFFPairCut[i]=new TH2D(sstr_name.str().c_str(),sstr_title.str().c_str(),((i%4)/2*(-5)+(i/4)*(!((i%4)/2))+7)*8,-0.5,((i%4)/2*(-5)+(i/4)*(!((i%4)/2))+7)*8-0.5,150,0,5000);
    sstr_name.str("");
    sstr_title.str("");

    sstr_name<<"Adc Vs Strip for all pairs of "<<"MwpcC"<<i/4+2<<(char)((i%4)/2+88)<<(char)((i%4%2)*6+76)<<" for Threshold, Total-Charge, First-Flag, Pair and Second-Flag Cuts";
    sstr_title<<sstr_name.str()<<";Strip Number; ADC Value";
    histAdcVsStThTCFFPairSFCut[i]=new TH2D(sstr_name.str().c_str(),sstr_title.str().c_str(),((i%4)/2*(-5)+(i/4)*(!((i%4)/2))+7)*8,-0.5,((i%4)/2*(-5)+(i/4)*(!((i%4)/2))+7)*8-0.5,150,0,5000);
    sstr_name.str("");
    sstr_title.str("");

    sstr_name<<"Position for all clusters of "<<"MwpcC"<<i/4+2<<(char)((i%4)/2+88)<<(char)((i%4%2)*6+76)<<" for No Cut";
    sstr_title<<sstr_name.str()<<";Position; Counts";
    histPosNoCut[i]=new TH1D(sstr_name.str().c_str(),sstr_title.str().c_str(),((i%4)/2*(-5)+(i/4)*(!((i%4)/2))+7)*8,0,((i%4)/2*(-5)+(i/4)*(!((i%4)/2))+7)*8);
    sstr_name.str("");
    sstr_title.str("");

    sstr_name<<"Position for all clusters of "<<"MwpcC"<<i/4+2<<(char)((i%4)/2+88)<<(char)((i%4%2)*6+76)<<" for Threshold and Total-Charge Cuts";
    sstr_title<<sstr_name.str()<<";Position; Counts";
    histPosThTCCut[i]=new TH1D(sstr_name.str().c_str(),sstr_title.str().c_str(),((i%4)/2*(-5)+(i/4)*(!((i%4)/2))+7)*8,0,((i%4)/2*(-5)+(i/4)*(!((i%4)/2))+7)*8);
    sstr_name.str("");
    sstr_title.str("");

    sstr_name<<"Position for all clusters of "<<"MwpcC"<<i/4+2<<(char)((i%4)/2+88)<<(char)((i%4%2)*6+76)<<" for Threshold, Total-Charge and First-Flag Cuts";
    sstr_title<<sstr_name.str()<<";Position; Counts";
    histPosThTCFFCut[i]=new TH1D(sstr_name.str().c_str(),sstr_title.str().c_str(),((i%4)/2*(-5)+(i/4)*(!((i%4)/2))+7)*8,0,((i%4)/2*(-5)+(i/4)*(!((i%4)/2))+7)*8);
    sstr_name.str("");
    sstr_title.str("");

    sstr_name<<"Position for all pairs of "<<"MwpcC"<<i/4+2<<(char)((i%4)/2+88)<<(char)((i%4%2)*6+76)<<" for Threshold, Total-Charge, First-Flag and Pair Cuts";
    sstr_title<<sstr_name.str()<<";Position; Counts";
    histPosThTCFFPairCut[i]=new TH1D(sstr_name.str().c_str(),sstr_title.str().c_str(),((i%4)/2*(-5)+(i/4)*(!((i%4)/2))+7)*8,0,((i%4)/2*(-5)+(i/4)*(!((i%4)/2))+7)*8);
    sstr_name.str("");
    sstr_title.str("");

    sstr_name<<"Position for all pairs of "<<"MwpcC"<<i/4+2<<(char)((i%4)/2+88)<<(char)((i%4%2)*6+76)<<" for Threshold, Total-Charge, First-Flag, Pair and Second-Flag Cuts";
    sstr_title<<sstr_name.str()<<";Position; Counts";
    histPosThTCFFPairSFCut[i]=new TH1D(sstr_name.str().c_str(),sstr_title.str().c_str(),((i%4)/2*(-5)+(i/4)*(!((i%4)/2))+7)*8,0,((i%4)/2*(-5)+(i/4)*(!((i%4)/2))+7)*8);
    sstr_name.str("");
    sstr_title.str("");
  }
  return 0;
}

Long_t Det_MWPC::startup_analysisExtractedCluster2(){
  getBranchObject("CaliMwpcInfo",(TObject **) &treeCali);
  gStyle->SetOptStat(0000);
  gStyle->SetOptFit(1111);

  return 0;
};


Long_t Det_MWPC::process_analysisExtractedCluster2(){
  for(int i=0;i<12;i++){
    for(UInt_t j=0;j<treeCali->noCutClusterStNum[i].size();j++){
      histPosNoCut[i]->Fill(treeCali->noCutClusterPosCentroid[i][j]);
      for(UInt_t k=0;k<treeCali->noCutClusterStNum[i][j].size();k++){
        histAdcVsStNoCut[i]->Fill(treeCali->noCutClusterStNum[i][j][k],treeCali->noCutClusterStCharge[i][j][k]);
      }
    }
    for(UInt_t j=0;j<treeCali->thTCCutClusterStNum[i].size();j++){
      histPosThTCCut[i]->Fill(treeCali->thTCCutClusterPosCentroid[i][j]);
      for(UInt_t k=0;k<treeCali->thTCCutClusterStNum[i][j].size();k++){
        histAdcVsStThTCCut[i]->Fill(treeCali->thTCCutClusterStNum[i][j][k],treeCali->thTCCutClusterStCharge[i][j][k]);
      }
    }
    if(treeCali->thTCCutLeftFlagFirst==1 || treeCali->thTCCutRightFlagFirst==1){
      for(UInt_t j=0;j<treeCali->thTCCutClusterStNum[i].size();j++){
        histPosThTCFFCut[i]->Fill(treeCali->thTCCutClusterPosCentroid[i][j]);
        for(UInt_t k=0;k<treeCali->thTCCutClusterStNum[i][j].size();k++){
          histAdcVsStThTCFFCut[i]->Fill(treeCali->thTCCutClusterStNum[i][j][k],treeCali->thTCCutClusterStCharge[i][j][k]);
        }
      }
    }
    if(treeCali->thTCPairCutLeftFlagFirst==1 || treeCali->thTCPairCutRightFlagFirst==1){
      for(UInt_t j=0;j<treeCali->thTCPairCutPairClusterStNum[i].size();j++){
        histPosThTCFFPairCut[i]->Fill(treeCali->thTCPairCutPairPosCentroid[i][j]);
        for(UInt_t k=0;k<treeCali->thTCPairCutPairClusterStNum[i][j].size();k++){
          histAdcVsStThTCFFPairCut[i]->Fill(treeCali->thTCPairCutPairClusterStNum[i][j][k],treeCali->thTCPairCutPairClusterStCharge[i][j][k]);
        }
      }
    }
    if((treeCali->thTCPairCutLeftFlagFirst==1 || treeCali->thTCPairCutRightFlagFirst==1) && (treeCali->thTCPairCutLeftFlagSecond==1 || treeCali->thTCPairCutRightFlagSecond==1)){
      for(UInt_t j=0;j<treeCali->thTCPairCutPairClusterStNum[i].size();j++){
        histPosThTCFFPairSFCut[i]->Fill(treeCali->thTCPairCutPairPosCentroid[i][j]);
        for(UInt_t k=0;k<treeCali->thTCPairCutPairClusterStNum[i][j].size();k++){
          histAdcVsStThTCFFPairSFCut[i]->Fill(treeCali->thTCPairCutPairClusterStNum[i][j][k],treeCali->thTCPairCutPairClusterStCharge[i][j][k]);
        }
      }
    }
  }
  return 0;
}

Long_t Det_MWPC::done_analysisExtractedCluster2(){
  TCanvas *canvAdcVsStNoCut=new TCanvas("Adc Vs Strip for all clusters for No Cut","Adc Vs Strip for all clusters for No Cut",2000,1500);
  canvAdcVsStNoCut->Divide(4,3);
  for(int i=0;i<12;i++){
    canvAdcVsStNoCut->cd(i+1);
    gPad->SetLogz();
    histAdcVsStNoCut[i]->Draw("colz");
  }
  canvAdcVsStNoCut->SaveAs("Adc Vs Strip for all clusters for No Cut.pdf");


  TCanvas *canvAdcVsStThTCCut=new TCanvas("Adc Vs Strip for all clusters for Threshold and Total-Charge Cuts","Adc Vs Strip for all clusters for Threshold and Total-Charge Cuts",2000,1500);
  canvAdcVsStThTCCut->Divide(4,3);
  for(int i=0;i<12;i++){
    canvAdcVsStThTCCut->cd(i+1);
    gPad->SetLogz();
    histAdcVsStThTCCut[i]->Draw("colz");
  }
  canvAdcVsStThTCCut->SaveAs("Adc Vs Strip for all clusters for Threshold and Total-Charge Cuts.pdf");

  TCanvas *canvAdcVsStThTCFFCut=new TCanvas("Adc Vs Strip for all clusters for Threshold, Total-Charge and First-Flag Cuts","Adc Vs Strip for all clusters for Threshold, Total-Charge and First-Flag Cuts",2000,1500);
  canvAdcVsStThTCFFCut->Divide(4,3);
  for(int i=0;i<12;i++){
    canvAdcVsStThTCFFCut->cd(i+1);
    gPad->SetLogz();
    histAdcVsStThTCFFCut[i]->Draw("colz");
  }
  canvAdcVsStThTCFFCut->SaveAs("Adc Vs Strip for all clusters for Threshold, Total-Charge and First-Flag Cuts.pdf");

  TCanvas *canvAdcVsStThTCFFPairCut=new TCanvas("Adc Vs Strip for all pairs for Threshold, Total-Charge, First-Flag and Pair Cuts","Adc Vs Strip for all pairs for Threshold, Total-Charge, First-Flag and Pair Cuts",2000,1500);
  canvAdcVsStThTCFFPairCut->Divide(4,3);
  for(int i=0;i<12;i++){
    canvAdcVsStThTCFFPairCut->cd(i+1);
    gPad->SetLogz();
    histAdcVsStThTCFFPairCut[i]->Draw("colz");
  }
  canvAdcVsStThTCFFPairCut->SaveAs("Adc Vs Strip for all pairs for Threshold, Total-Charge, First-Flag and Pair Cuts.pdf");

  TCanvas *canvAdcVsStThTCFFPairSFCut=new TCanvas("Adc Vs Strip for all pairs for Threshold, Total-Charge, First-Flag, Pair and Second-Flag Cuts","Adc Vs Strip for all pairs for Threshold, Total-Charge, First-Flag, Pair and Second-Flag Cuts",2000,1500);
  canvAdcVsStThTCFFPairSFCut->Divide(4,3);
  for(int i=0;i<12;i++){
    canvAdcVsStThTCFFPairSFCut->cd(i+1);
    gPad->SetLogz();
    histAdcVsStThTCFFPairSFCut[i]->Draw("colz");
  }
  canvAdcVsStThTCFFPairSFCut->SaveAs("Adc Vs Strip for all pairs for Threshold, Total-Charge, First-Flag, Pair and Second-Flag Cuts.pdf");

  TCanvas *canvPosNoCut=new TCanvas("Position for all clusters for No Cut","Position for all clusters for No Cut",2000,900);
  canvPosNoCut->Divide(4,3);
  for(int i=0;i<12;i++){
    canvPosNoCut->cd(i+1);
    gPad->SetLogz();
    histPosNoCut[i]->Draw("colz");
  }
  canvPosNoCut->SaveAs("Position for all clusters for No Cut.pdf");

  TCanvas *canvPosThTCCut=new TCanvas("Position for all clusters for Threshold and Total-Charge Cuts","Position for all clusters for Threshold and Total-Charge Cuts",2000,900);
  canvPosThTCCut->Divide(4,3);
  for(int i=0;i<12;i++){
    canvPosThTCCut->cd(i+1);
    gPad->SetLogz();
    histPosThTCCut[i]->Draw("colz");
  }
  canvPosThTCCut->SaveAs("Position for all clusters for Threshold and Total-Charge Cuts.pdf");

  TCanvas *canvPosThTCFFCut=new TCanvas("Position for all clusters for Threshold, Total-Charge and First-Flag Cuts","Position for all clusters for Threshold, Total-Charge and First-Flag Cuts",2000,900);
  canvPosThTCFFCut->Divide(4,3);
  for(int i=0;i<12;i++){
    canvPosThTCFFCut->cd(i+1);
    gPad->SetLogz();
    histPosThTCFFCut[i]->Draw("colz");
  }
  canvPosThTCFFCut->SaveAs("Position for all clusters for Threshold, Total-Charge and First-Flag Cuts.pdf");

  TCanvas *canvPosThTCFFPairCut=new TCanvas("Position for all pairs for Threshold, Total-Charge, First-Flag and Pair Cuts","Position for all pairs for Threshold, Total-Charge, First-Flag and Pair Cuts",2000,900);
  canvPosThTCFFPairCut->Divide(4,3);
  for(int i=0;i<12;i++){
    canvPosThTCFFPairCut->cd(i+1);
    gPad->SetLogz();
    histPosThTCFFPairCut[i]->Draw("colz");
  }
  canvPosThTCFFPairCut->SaveAs("Position for all pairs for Threshold, Total-Charge, First-Flag and Pair Cuts.pdf");

  TCanvas *canvPosThTCFFPairSFCut=new TCanvas("Position for all pairs for Threshold, Total-Charge, First-Flag, Pair and Second-Flag Cuts","Position for all pairs for Threshold, Total-Charge, First-Flag, Pair and Second-Flag Cuts",2000,900);
  canvPosThTCFFPairSFCut->Divide(4,3);
  for(int i=0;i<12;i++){
    canvPosThTCFFPairSFCut->cd(i+1);
    gPad->SetLogz();
    histPosThTCFFPairSFCut[i]->Draw("colz");
  }
  canvPosThTCFFPairSFCut->SaveAs("Position for all pairs for Threshold, Total-Charge, First-Flag, Pair and Second-Flag Cuts.pdf");

  return 0;
};

