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

static string inputFilePath(getenv("mwpc_inputFiles"));
static string inputFile=inputFilePath+"/gapNumTable/run3990-4059";
static ifstream input(inputFile.c_str());
static int runNum, eventNum, gapNumTOF1, gapNumTOF2;

TH2D *histDiffVsSum[6][16];
TH2D *histDiffVsSumSt0313[6], *histDiffVsSumSt2429[6];
Long_t Det_MWPC::histos_handleDeadStrip(){
  ostringstream sstr_name,sstr_title;
  for(int i=0;i<6;i++){
    if(i!=1){
      for(int j=0;j<16;j++){
        sstr_name<<"Diff vs sum of XY total-charge for X1Y1 for strip "<<j+12<<" of gap "<<i+1<<" of C4XR before gain calibration";
        sstr_title<<"Strip "<<j+12<<";TC_{X}+TC_{Y};TC_{X}-TC_{Y}";
        histDiffVsSum[i][j]=new TH2D(sstr_name.str().c_str(),sstr_title.str().c_str(),200,0,5000,200,-1000,3000);
        sstr_name.str("");
        sstr_title.str("");
      }

      sstr_name<<"Diff vs sum of XY total-charge for X1Y1 for strip 03 - 13 of gap "<<i+1<<" of C4XR before gain calibration";
      sstr_title<<"Gap "<<i+1<<";TC_{X}+TC_{Y};TC_{X}-TC_{Y}";
      histDiffVsSumSt0313[i]=new TH2D(sstr_name.str().c_str(),sstr_title.str().c_str(),200,0,5000,200,-1000,3000);
      sstr_name.str("");
      sstr_title.str("");

      sstr_name<<"Diff vs sum of XY total-charge for X1Y1 for strip 24 - 29 of gap "<<i+1<<" of C4XR before gain calibration";
      sstr_title<<"Gap "<<i+1<<";TC_{X}+TC_{Y};TC_{X}-TC_{Y}";
      histDiffVsSumSt2429[i]=new TH2D(sstr_name.str().c_str(),sstr_title.str().c_str(),200,0,5000,200,-1000,3000);
      sstr_name.str("");
      sstr_title.str("");
    }
    else{
      for(int j=0;j<16;j++){
        sstr_name<<"Diff vs sum of XY total-charge for X1Y1 for strip "<<j+12<<" of gap "<<i+1<<" of C4XR before gain calibration";
        sstr_title<<"Strip "<<j+12<<";TC_{X}+TC_{Y};TC_{X}-TC_{Y}";
        histDiffVsSum[i][j]=new TH2D(sstr_name.str().c_str(),sstr_title.str().c_str(),200,0,8000,200,-1000,5000);
        sstr_name.str("");
        sstr_title.str("");
      }

      sstr_name<<"Diff vs sum of XY total-charge for X1Y1 for strip 03 - 13 of gap "<<i+1<<" of C4XR before gain calibration";
      sstr_title<<"Gap "<<i+1<<";TC_{X}+TC_{Y};TC_{X}-TC_{Y}";
      histDiffVsSumSt0313[i]=new TH2D(sstr_name.str().c_str(),sstr_title.str().c_str(),200,0,8000,200,-1000,5000);
      sstr_name.str("");
      sstr_title.str("");

      sstr_name<<"Diff vs sum of XY total-charge for X1Y1 for strip 24 - 29 of gap "<<i+1<<" of C4XR before gain calibration";
      sstr_title<<"Gap "<<i+1<<";TC_{X}+TC_{Y};TC_{X}-TC_{Y}";
      histDiffVsSumSt2429[i]=new TH2D(sstr_name.str().c_str(),sstr_title.str().c_str(),200,0,8000,200,-1000,5000);
      sstr_name.str("");
      sstr_title.str("");
    }
  }
  
  return 0;
}

Long_t Det_MWPC::startup_handleDeadStrip(){
  getBranchObject("CaliMwpcInfo",(TObject **) &treeCali);
  gStyle->SetOptStat(0000);
  gStyle->SetOptFit(1111);

  return 0;
};

Long_t Det_MWPC::process_handleDeadStrip(){
  input>>runNum>>eventNum>>gapNumTOF1>>gapNumTOF2;
  if((UInt_t)runNum!=treeCali->run || (UInt_t)eventNum!=treeCali->event){
    cout<<"Error: Run number or event number are not consistent between two files!"<<endl;
    cout<<"For MWPC file: "<<"Run number - "<<treeCali->run<<"and  Event number - "<<treeCali->event<<endl;
    cout<<"For TRIUMF file: "<<"Run number - "<<runNum<<"and  Event number - "<<eventNum<<endl;
  }

  MwpcEvent *event=new MwpcEvent(72,1);
  if(treeCali->thTCCutNumCluster[9]==1 && treeCali->thTCCutNumCluster[11]==1 && gapNumTOF2!=-1){
    for(int i=0;i<6;i++){
      if(gapNumTOF2==i+1){
        for(int j=0;j<16;j++){
          if(event->getEachClusterStNumWithHighestADC(treeCali->thTCCutClusterStNum[9][0],treeCali->thTCCutClusterStCharge[9][0])==12+j){
            histDiffVsSum[i][j]->Fill(treeCali->thTCCutClusterTC[9][0]+treeCali->thTCCutClusterTC[11][0],treeCali->thTCCutClusterTC[9][0]-treeCali->thTCCutClusterTC[11][0]);
          }
        if(event->getEachClusterStNumWithHighestADC(treeCali->thTCCutClusterStNum[9][0],treeCali->thTCCutClusterStCharge[9][0])>=3 && event->getEachClusterStNumWithHighestADC(treeCali->thTCCutClusterStNum[9][0],treeCali->thTCCutClusterStCharge[9][0])<=13) histDiffVsSumSt0313[i]->Fill(treeCali->thTCCutClusterTC[9][0]+treeCali->thTCCutClusterTC[11][0],treeCali->thTCCutClusterTC[9][0]-treeCali->thTCCutClusterTC[11][0]);
        if(event->getEachClusterStNumWithHighestADC(treeCali->thTCCutClusterStNum[9][0],treeCali->thTCCutClusterStCharge[9][0])>=24 && event->getEachClusterStNumWithHighestADC(treeCali->thTCCutClusterStNum[9][0],treeCali->thTCCutClusterStCharge[9][0])<=29) histDiffVsSumSt2429[i]->Fill(treeCali->thTCCutClusterTC[9][0]+treeCali->thTCCutClusterTC[11][0],treeCali->thTCCutClusterTC[9][0]-treeCali->thTCCutClusterTC[11][0]);
        }
      }
    }
  }

  return 0;
}

Long_t Det_MWPC::done_handleDeadStrip(){
  TFile *histroot=new TFile("hist.root","recreate");
  ostringstream sstr;
  TCanvas *canvDiffVsSum[6];
  for(int i=0;i<6;i++){
    sstr<<"Diff vs sum of XY total-charge for X1Y1 of gap "<<i+1<<" of C4YR";
    canvDiffVsSum[i]=new TCanvas(sstr.str().c_str(),sstr.str().c_str(),1200,1200);
    sstr.str("");
    canvDiffVsSum[i]->Divide(4,4);
    for(int j=0;j<16;j++){
      canvDiffVsSum[i]->cd(j+1);
      histDiffVsSum[i][j]->Draw("colz");
      histDiffVsSum[i][j]->Write();
    }
    sstr<<"Diff vs sum of XY total-charge for X1Y1 of gap "<<i+1<<" of C4YR.pdf";
    canvDiffVsSum[i]->SaveAs(sstr.str().c_str());
    sstr.str("");
  }

  TCanvas *canvDiffVsSumSt0313=new TCanvas("Diff vs sum of XY total-charge for X1Y1 of strip 03 - 13","Diff vs sum of XY total-charge for X1Y1 of strip 03 - 13",800,1200);
  canvDiffVsSumSt0313->Divide(2,3);
  for(int i=0;i<6;i++){
    canvDiffVsSumSt0313->cd(i+1);
    gPad->SetLogz();
    histDiffVsSumSt0313[i]->Draw("colz");
    histDiffVsSumSt0313[i]->Write();
  }
  canvDiffVsSumSt0313->SaveAs("Diff vs sum of XY total-charge for X1Y1 of strip 03 - 13.pdf");

  TCanvas *canvDiffVsSumSt2429=new TCanvas("Diff vs sum of XY total-charge for X1Y1 of strip 24 - 29","Diff vs sum of XY total-charge for X1Y1 of strip 24 - 29",800,1200);
  canvDiffVsSumSt2429->Divide(2,3);
  for(int i=0;i<6;i++){
    canvDiffVsSumSt2429->cd(i+1);
    gPad->SetLogz();
    histDiffVsSumSt2429[i]->Draw("colz");
    histDiffVsSumSt2429[i]->Write();
  }
  canvDiffVsSumSt2429->SaveAs("Diff vs sum of XY total-charge for X1Y1 of strip 24 - 29.pdf");

  histroot->Write();

  return 0;
};

