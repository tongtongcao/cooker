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
static string inputFileC3XR=inputFilePath+"/gain-weight-C3XR";
static string inputFileC3YR=inputFilePath+"/gain-weight-C3YR";
static ifstream inputC3XR(inputFileC3XR.c_str());
static ifstream inputC3YR(inputFileC3YR.c_str());
static double gainWeightC3XR[6][64], gainWeightC3YR[6][16];
static string inputFileNormC3R=inputFilePath+"/normalized-factor-C3R";
static ifstream inputNormC3R(inputFileNormC3R.c_str());
static double normC3R[6];

TH1D *histRatioDoubY[6][16];
TH1D *histRatioTripY1[6][16];

Long_t Det_MWPC::histos_posPuzzle(){
  ostringstream sstr_name,sstr_title;
  for(int i=0;i<6;i++){
    for(int j=0;j<16;j++){
      sstr_name<<"Ratio of low and high ADC for Doublet of X1Y1 for strip "<<j<<" of gap "<<i+1<<" of C3XR after gain calibration";
      sstr_title<<"Strip "<<j<<";Ratio;Counts";
      histRatioDoubY[i][j]=new TH1D(sstr_name.str().c_str(),sstr_title.str().c_str(),100,0,1);
      sstr_name.str("");
      sstr_title.str("");
    }
  }   

  return 0;
}

Long_t Det_MWPC::startup_posPuzzle(){
  getBranchObject("CaliMwpcInfo",(TObject **) &treeCali);
  gStyle->SetOptStat(0000);
  gStyle->SetOptFit(1111);
  for(int i=0;i<6;i++){
    for(int j=0;j<64;j++){
      inputC3XR>>gainWeightC3XR[i][j];
    }
  }
  for(int i=0;i<6;i++){
    for(int j=0;j<16;j++){
      inputC3YR>>gainWeightC3YR[i][j];
    }
  }
  for(int i=0;i<6;i++){
    inputNormC3R>>normC3R[i];
  }

  return 0;
};

Long_t Det_MWPC::process_posPuzzle(){
  input>>runNum>>eventNum>>gapNumTOF1>>gapNumTOF2;
  if((UInt_t)runNum!=treeCali->run || (UInt_t)eventNum!=treeCali->event){
    cout<<"Error: Run number or event number are not consistent between two files!"<<endl;
    cout<<"For MWPC file: "<<"Run number - "<<treeCali->run<<"and  Event number - "<<treeCali->event<<endl;
    cout<<"For TRIUMF file: "<<"Run number - "<<runNum<<"and  Event number - "<<eventNum<<endl;
  }

  vector<double> stCharge;
  if(treeCali->thTCCutNumCluster[5]==1 && treeCali->thTCCutNumCluster[7]==1 && gapNumTOF2!=-1){
    for(int iGap=0;iGap<6;iGap++){
      if(gapNumTOF2==iGap+1){
        for(int i=0;i<16;i++){
          if(treeCali->thTCCutClusterPosCentroid[7][0]>=i && treeCali->thTCCutClusterPosCentroid[7][0]<i+1){
            if(treeCali->thTCCutClusterStNum[7][0].size()==2){
              for(int j=0;j<treeCali->thTCCutClusterStCharge[7][0].size();j++){
                stCharge.push_back(treeCali->thTCCutClusterStCharge[7][0][j]*gainWeightC3YR[iGap][treeCali->thTCCutClusterStNum[7][0][j]]);
              }
              if(stCharge[0]<=stCharge[1]) histRatioDoubY[iGap][i]->Fill(stCharge[0]/stCharge[1]);
              else histRatioDoubY[iGap][i]->Fill(stCharge[1]/stCharge[0]);
              stCharge.clear();
            }
          }
        }
      }
    }
  }


  return 0;
}

Long_t Det_MWPC::done_posPuzzle(){
  ostringstream sstr;
  TCanvas *canvRatioDoubY[6];
  for(int i=0;i<6;i++){
    sstr<<"Ratio of low and high ADC for Doublet of X1Y1 of gap "<<i+1<<" of C3YR";
    canvRatioDoubY[i]=new TCanvas(sstr.str().c_str(),sstr.str().c_str(),2000,1200);
    sstr.str("");
    canvRatioDoubY[i]->Divide(4,4);
    for(int j=0;j<16;j++){
      canvRatioDoubY[i]->cd(j+1);
      histRatioDoubY[i][j]->Draw();
    }
    sstr<<"Ratio of low and high ADC for Doublet of X1Y1 of gap "<<i+1<<" of C3YR.pdf";
    canvRatioDoubY[i]->SaveAs(sstr.str().c_str());
    sstr.str("");
  }

  return 0;
};

