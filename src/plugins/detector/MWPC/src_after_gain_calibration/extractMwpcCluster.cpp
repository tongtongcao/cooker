////// To use cut functions of diff vs. sum of XY-total-charge, a full directory should be given since the functions are saved in a root file. Here, the full directory is /home/caot/trek-e36/cooker/src/plugins/detector/MWPC/src/MwpcDiffVsSumTotalChargeCutFunction.root 
#include <Det_MWPC.h>
#include <cmath>
#include <vector>
#include <fstream>
#include <iostream>
#include <TStyle.h>
#include <sstream>
#include "TLine.h"
#include <string>
#include "TCanvas.h"
#include "TText.h"
#include "TLegend.h"
#include "TLine.h"
#include "TLatex.h"
#include "TTree.h"
#include "TF1.h"
#include "MwpcPedestalValue.h"
#include "MwpcClusterThresholdSecondApproach.h"
#include "MwpcClusterThresholdFromTotalCharge.h"
using namespace std;

static  TFile *fileCutFunc=new TFile("/home/caot/trek-e36/cooker/src/plugins/detector/MWPC/src/MwpcDiffVsSumTotalChargeCutFunction.root");
static  TF1 *cutFuncUp[6],*cutFuncLow[6];
static TH2D *test;
Long_t Det_MWPC::histos_extractCluster(){

  test=new TH2D("test","test",200,0,5000,200,-1000,1500);

  ostringstream sstr;
  for(int i=0;i<6;i++){
    sstr<<"C"<<i/2+2<<(char)((i%2)*6+76)<<"_up";
    cutFuncUp[i]=(TF1*)fileCutFunc->Get(sstr.str().c_str());
    sstr.str("");
    sstr<<"C"<<i/2+2<<(char)((i%2)*6+76)<<"_low";
    cutFuncLow[i]=(TF1*)fileCutFunc->Get(sstr.str().c_str());
    sstr.str("");
  }

  return 0;
}

Long_t Det_MWPC::startup_extractCluster(){
  getBranchObject("RawMwpcInfo",(TObject **) &treeRaw);
  treeCali=0;
  out->Branch("CaliMwpcInfo",&treeCali,160000,0);
  gStyle->SetOptStat(0000);
  gStyle->SetOptFit(1111);
  return 0;
};


Long_t Det_MWPC::process_extractCluster(){
  int* data[12]={treeRaw->ADC_C2X_L,treeRaw->ADC_C2X_R,treeRaw->ADC_C2Y_L,treeRaw->ADC_C2Y_R,treeRaw->ADC_C3X_L,treeRaw->ADC_C3X_R,treeRaw->ADC_C3Y_L,treeRaw->ADC_C3Y_R,treeRaw->ADC_C4X_L,treeRaw->ADC_C4X_R,treeRaw->ADC_C4Y_L,treeRaw->ADC_C4Y_R};
  const double* pedValue[12]={pedValueC2XL,pedValueC2XR,pedValueC2YL,pedValueC2YR,pedValueC3XL,pedValueC3XR,pedValueC3YL,pedValueC3YR,pedValueC4XL,pedValueC4XR,pedValueC4YL,pedValueC4YR};
  const double* pedSigma[12]={pedSigmaC2XL,pedSigmaC2XR,pedSigmaC2YL,pedSigmaC2YR,pedSigmaC3XL,pedSigmaC3XR,pedSigmaC3YL,pedSigmaC3YR,pedSigmaC4XL,pedSigmaC4XR,pedSigmaC4YL,pedSigmaC4YR}; 
  const double* totalChargeCut[12]={totChargeCutC2XL,totChargeCutC2XR,totChargeCutC2YL,totChargeCutC2YR,totChargeCutC3XL,totChargeCutC3XR,totChargeCutC3YL,totChargeCutC3YR,totChargeCutC4XL,totChargeCutC4XR,totChargeCutC4YL,totChargeCutC4YR};
  MwpcEvent *event[12];
  vector<vector<int> > noCutClusterStNum[12],thTCCutClusterStNum[12];
  vector<vector<double> > noCutClusterStCharge[12],thTCCutClusterStCharge[12];
  double *noCutClusterTC[12],*thTCCutClusterTC[12];
  for(int i=0;i<12;i++){
    event[i]=new MwpcEvent(((i%4)/2*(-5)+(i/4)*(!((i%4)/2))+7)*8,!((i/4==0)&&((i%4/2))==0));
    noCutClusterStNum[i]=event[i]->clusterStripNumber(event[i],data[i],pedValue[i],pedSigma[i],0,totalChargeCut[i],0,singMiniSigmaClusterSecondApproach[i],0,doubMiniSigmaXClusterSecondApproach[i],doubMiniTotalChargeXCluster[i]);
    noCutClusterStCharge[i]=event[i]->clusterStripCharge(event[i],data[i],pedValue[i],pedSigma[i],0,totalChargeCut[i],0,singMiniSigmaClusterSecondApproach[i],0,doubMiniSigmaXClusterSecondApproach[i],doubMiniTotalChargeXCluster[i]);
    noCutClusterTC[i]=event[i]->getClusterTotCharge(event[i],data[i],pedValue[i],pedSigma[i],0,totalChargeCut[i],0,singMiniSigmaClusterSecondApproach[i],0,doubMiniSigmaXClusterSecondApproach[i],doubMiniTotalChargeXCluster[i]);
    thTCCutClusterStNum[i]=event[i]->clusterStripNumber(event[i],data[i],pedValue[i],pedSigma[i],1,totalChargeCut[i],1,singMiniSigmaClusterSecondApproach[i],1,doubMiniSigmaXClusterSecondApproach[i],doubMiniTotalChargeXCluster[i]);
    thTCCutClusterStCharge[i]=event[i]->clusterStripCharge(event[i],data[i],pedValue[i],pedSigma[i],1,totalChargeCut[i],1,singMiniSigmaClusterSecondApproach[i],1,doubMiniSigmaXClusterSecondApproach[i],doubMiniTotalChargeXCluster[i]);
    thTCCutClusterTC[i]=event[i]->getClusterTotCharge(event[i],data[i],pedValue[i],pedSigma[i],1,totalChargeCut[i],1,singMiniSigmaClusterSecondApproach[i],1,doubMiniSigmaXClusterSecondApproach[i],doubMiniTotalChargeXCluster[i]);
  }

  treeCali->init(); // Init data memebers
  //// Common variables
  treeCali->run=treeRaw->run;
  treeCali->event=treeRaw->event;
  double* adcValuePedSub;
  for(int i=0;i<12;i++){
    adcValuePedSub=event[i]->getAdcValue(data[i],pedValue[i]);
    for(int j=0;j<((i%4)/2*(-5)+(i/4)*(!((i%4)/2))+7)*8;j++){
      treeCali->adcValueOriginal[i].push_back(data[i][j]);
      treeCali->adcValuePedSub[i].push_back(adcValuePedSub[j]);
    }
  }

  //// Specified variables for no cut
  for(int i=0;i<12;i++){
    treeCali->noCutNumCluster[i]=noCutClusterStNum[i].size();
    treeCali->noCutClusterStNum[i]=noCutClusterStNum[i];
    treeCali->noCutClusterStCharge[i]=noCutClusterStCharge[i];   
    for(UInt_t j=0;j<noCutClusterStNum[i].size();j++){
      treeCali->noCutClusterTC[i].push_back(noCutClusterTC[i][j]);
      treeCali->noCutClusterNumSt[i].push_back(noCutClusterStNum[i].size());
      treeCali->noCutClusterPosCentroid[i].push_back(event[i]->getEachClusterPositionCentroid(treeCali->noCutClusterStNum[i][j],treeCali->noCutClusterStCharge[i][j]));
      treeCali->noCutClusterPosErrCentroid[i].push_back(event[i]->getEachClusterPositionErrorCentroid(treeCali->noCutClusterStCharge[i][j]));
    }
  }

  if(treeCali->noCutNumCluster[0]>0 && treeCali->noCutNumCluster[2]>0 && treeCali->noCutNumCluster[4]>0 && treeCali->noCutNumCluster[6]>0 && treeCali->noCutNumCluster[8]>0 && treeCali->noCutNumCluster[10]>0) treeCali->noCutLeftFlagFirst=1;
  if(treeCali->noCutNumCluster[1]>0 && treeCali->noCutNumCluster[3]>0 && treeCali->noCutNumCluster[5]>0 && treeCali->noCutNumCluster[7]>0 && treeCali->noCutNumCluster[9]>0 && treeCali->noCutNumCluster[11]>0) treeCali->noCutRightFlagFirst=1;

  for(int i=0;i<6;i++){
    for(UInt_t m=0;m<treeCali->noCutNumCluster[i/2*4+i%2];m++){
      for(UInt_t n=0;n<treeCali->noCutNumCluster[i/2*4+i%2+2];n++){
        treeCali->noCutPairClusterNum[i/2*4+i%2].push_back(m);
        treeCali->noCutPairClusterNum[i/2*4+i%2+2].push_back(n);
        treeCali->noCutPairClusterTC[i/2*4+i%2].push_back(treeCali->noCutClusterTC[i/2*4+i%2][m]);
        treeCali->noCutPairClusterTC[i/2*4+i%2+2].push_back(treeCali->noCutClusterTC[i/2*4+i%2+2][n]);
        treeCali->noCutPairClusterNumSt[i/2*4+i%2].push_back(treeCali->noCutClusterNumSt[i/2*4+i%2][m]);
        treeCali->noCutPairClusterNumSt[i/2*4+i%2+2].push_back(treeCali->noCutClusterNumSt[i/2*4+i%2+2][n]);
        treeCali->noCutPairClusterStNum[i/2*4+i%2].push_back(treeCali->noCutClusterStNum[i/2*4+i%2][m]);
        treeCali->noCutPairClusterStNum[i/2*4+i%2+2].push_back(treeCali->noCutClusterStNum[i/2*4+i%2+2][n]);
        treeCali->noCutPairClusterStCharge[i/2*4+i%2].push_back(treeCali->noCutClusterStCharge[i/2*4+i%2][m]);
        treeCali->noCutPairClusterStCharge[i/2*4+i%2+2].push_back(treeCali->noCutClusterStCharge[i/2*4+i%2+2][n]);
        treeCali->noCutPairPosCentroid[i/2*4+i%2].push_back(event[i/2*4+i%2]->getEachClusterPositionCentroid(treeCali->noCutClusterStNum[i/2*4+i%2][m],treeCali->noCutClusterStCharge[i/2*4+i%2][m]));
        treeCali->noCutPairPosCentroid[i/2*4+i%2+2].push_back(event[i/2*4+i%2+2]->getEachClusterPositionCentroid(treeCali->noCutClusterStNum[i/2*4+i%2+2][n],treeCali->noCutClusterStCharge[i/2*4+i%2+2][n]));
        treeCali->noCutPairPosErrCentroid[i/2*4+i%2].push_back(event[i/2*4+i%2]->getEachClusterPositionErrorCentroid(treeCali->noCutClusterStCharge[i/2*4+i%2][m]));
        treeCali->noCutPairPosErrCentroid[i/2*4+i%2+2].push_back(event[i/2*4+i%2+2]->getEachClusterPositionErrorCentroid(treeCali->noCutClusterStCharge[i/2*4+i%2+2][n]));
      }
    }
  treeCali->noCutNumPair[i/2*4+i%2]=treeCali->noCutPairClusterNum[i/2*4+i%2].size();
  treeCali->noCutNumPair[i/2*4+i%2+2]=treeCali->noCutPairClusterNum[i/2*4+i%2+2].size();
  }
  if(treeCali->noCutNumPair[0]>0 && treeCali->noCutNumPair[4]>0 && treeCali->noCutNumPair[8]>0) treeCali->noCutLeftFlagSecond=1;
  if(treeCali->noCutNumPair[1]>0 && treeCali->noCutNumPair[5]>0 && treeCali->noCutNumPair[9]>0) treeCali->noCutRightFlagSecond=1;

  //// Specified variables for threshold and total-charge cuts
  for(int i=0;i<12;i++){
    treeCali->thTCCutNumCluster[i]=thTCCutClusterStNum[i].size();
    treeCali->thTCCutClusterStNum[i]=thTCCutClusterStNum[i];
    treeCali->thTCCutClusterStCharge[i]=thTCCutClusterStCharge[i];
    for(UInt_t j=0;j<thTCCutClusterStNum[i].size();j++){
      treeCali->thTCCutClusterTC[i].push_back(thTCCutClusterTC[i][j]);
      treeCali->thTCCutClusterNumSt[i].push_back(thTCCutClusterStNum[i].size());
      treeCali->thTCCutClusterPosCentroid[i].push_back(event[i]->getEachClusterPositionCentroid(treeCali->thTCCutClusterStNum[i][j],treeCali->thTCCutClusterStCharge[i][j]));
      treeCali->thTCCutClusterPosErrCentroid[i].push_back(event[i]->getEachClusterPositionErrorCentroid(treeCali->thTCCutClusterStCharge[i][j]));
    }
  }

  if(treeCali->thTCCutNumCluster[0]>0 && treeCali->thTCCutNumCluster[2]>0 && treeCali->thTCCutNumCluster[4]>0 && treeCali->thTCCutNumCluster[6]>0 && treeCali->thTCCutNumCluster[8]>0 && treeCali->thTCCutNumCluster[10]>0) treeCali->thTCCutLeftFlagFirst=1;
  if(treeCali->thTCCutNumCluster[1]>0 && treeCali->thTCCutNumCluster[3]>0 && treeCali->thTCCutNumCluster[5]>0 && treeCali->thTCCutNumCluster[7]>0 && treeCali->thTCCutNumCluster[9]>0 && treeCali->thTCCutNumCluster[11]>0) treeCali->thTCCutRightFlagFirst=1;

  for(int i=0;i<6;i++){
    for(UInt_t m=0;m<treeCali->thTCCutNumCluster[i/2*4+i%2];m++){
      for(UInt_t n=0;n<treeCali->thTCCutNumCluster[i/2*4+i%2+2];n++){
        treeCali->thTCCutPairClusterNum[i/2*4+i%2].push_back(m);
        treeCali->thTCCutPairClusterNum[i/2*4+i%2+2].push_back(n);
        treeCali->thTCCutPairClusterTC[i/2*4+i%2].push_back(treeCali->thTCCutClusterTC[i/2*4+i%2][m]);
        treeCali->thTCCutPairClusterTC[i/2*4+i%2+2].push_back(treeCali->thTCCutClusterTC[i/2*4+i%2+2][n]);
        treeCali->thTCCutPairClusterNumSt[i/2*4+i%2].push_back(treeCali->thTCCutClusterNumSt[i/2*4+i%2][m]);
        treeCali->thTCCutPairClusterNumSt[i/2*4+i%2+2].push_back(treeCali->thTCCutClusterNumSt[i/2*4+i%2+2][n]);
        treeCali->thTCCutPairClusterStNum[i/2*4+i%2].push_back(treeCali->thTCCutClusterStNum[i/2*4+i%2][m]);
        treeCali->thTCCutPairClusterStNum[i/2*4+i%2+2].push_back(treeCali->thTCCutClusterStNum[i/2*4+i%2+2][n]);
        treeCali->thTCCutPairClusterStCharge[i/2*4+i%2].push_back(treeCali->thTCCutClusterStCharge[i/2*4+i%2][m]);
        treeCali->thTCCutPairClusterStCharge[i/2*4+i%2+2].push_back(treeCali->thTCCutClusterStCharge[i/2*4+i%2+2][n]);
        treeCali->thTCCutPairPosCentroid[i/2*4+i%2].push_back(event[i/2*4+i%2]->getEachClusterPositionCentroid(treeCali->thTCCutClusterStNum[i/2*4+i%2][m],treeCali->thTCCutClusterStCharge[i/2*4+i%2][m]));
        treeCali->thTCCutPairPosCentroid[i/2*4+i%2+2].push_back(event[i/2*4+i%2+2]->getEachClusterPositionCentroid(treeCali->thTCCutClusterStNum[i/2*4+i%2+2][n],treeCali->thTCCutClusterStCharge[i/2*4+i%2+2][n]));
        treeCali->thTCCutPairPosErrCentroid[i/2*4+i%2].push_back(event[i/2*4+i%2]->getEachClusterPositionErrorCentroid(treeCali->thTCCutClusterStCharge[i/2*4+i%2][m]));
        treeCali->thTCCutPairPosErrCentroid[i/2*4+i%2+2].push_back(event[i/2*4+i%2+2]->getEachClusterPositionErrorCentroid(treeCali->thTCCutClusterStCharge[i/2*4+i%2+2][n]));
      }
    }
  treeCali->thTCCutNumPair[i/2*4+i%2]=treeCali->thTCCutPairClusterNum[i/2*4+i%2].size();
  treeCali->thTCCutNumPair[i/2*4+i%2+2]=treeCali->thTCCutPairClusterNum[i/2*4+i%2+2].size();
  }
  if(treeCali->thTCCutNumPair[0]>0 && treeCali->thTCCutNumPair[4]>0 && treeCali->thTCCutNumPair[8]>0) treeCali->thTCCutLeftFlagSecond=1;
  if(treeCali->thTCCutNumPair[1]>0 && treeCali->thTCCutNumPair[5]>0 && treeCali->thTCCutNumPair[9]>0) treeCali->thTCCutRightFlagSecond=1;

  ////Specified variables for pair cut
  for(int i=0;i<12;i++){
    treeCali->pairCutNumCluster[i]=noCutClusterStNum[i].size();
    treeCali->pairCutClusterStNum[i]=noCutClusterStNum[i];
    treeCali->pairCutClusterStCharge[i]=noCutClusterStCharge[i];
    for(UInt_t j=0;j<noCutClusterStNum[i].size();j++){
      treeCali->pairCutClusterTC[i].push_back(noCutClusterTC[i][j]);
      treeCali->pairCutClusterNumSt[i].push_back(noCutClusterStNum[i].size());
      treeCali->pairCutClusterPosCentroid[i].push_back(event[i]->getEachClusterPositionCentroid(treeCali->pairCutClusterStNum[i][j],treeCali->pairCutClusterStCharge[i][j]));
      treeCali->pairCutClusterPosErrCentroid[i].push_back(event[i]->getEachClusterPositionErrorCentroid(treeCali->pairCutClusterStCharge[i][j]));
    }
  }

  if(treeCali->pairCutNumCluster[0]>0 && treeCali->pairCutNumCluster[2]>0 && treeCali->pairCutNumCluster[4]>0 && treeCali->pairCutNumCluster[6]>0 && treeCali->pairCutNumCluster[8]>0 && treeCali->pairCutNumCluster[10]>0) treeCali->pairCutLeftFlagFirst=1;
  if(treeCali->pairCutNumCluster[1]>0 && treeCali->pairCutNumCluster[3]>0 && treeCali->pairCutNumCluster[5]>0 && treeCali->pairCutNumCluster[7]>0 && treeCali->pairCutNumCluster[9]>0 && treeCali->pairCutNumCluster[11]>0) treeCali->pairCutRightFlagFirst=1;

  for(int i=0;i<6;i++){
    for(UInt_t m=0;m<treeCali->pairCutNumCluster[i/2*4+i%2];m++){
      for(UInt_t n=0;n<treeCali->pairCutNumCluster[i/2*4+i%2+2];n++){
        if(treeCali->pairCutClusterTC[i/2*4+i%2][m]-treeCali->pairCutClusterTC[i/2*4+i%2+2][n]>cutFuncLow[i]->Eval(treeCali->pairCutClusterTC[i/2*4+i%2][m]+treeCali->pairCutClusterTC[i/2*4+i%2+2][n]) && treeCali->pairCutClusterTC[i/2*4+i%2][m]-treeCali->pairCutClusterTC[i/2*4+i%2+2][n]<cutFuncUp[i]->Eval(treeCali->pairCutClusterTC[i/2*4+i%2][m]+treeCali->pairCutClusterTC[i/2*4+i%2+2][n])){
          treeCali->pairCutPairClusterNum[i/2*4+i%2].push_back(m);
          treeCali->pairCutPairClusterNum[i/2*4+i%2+2].push_back(n);
          treeCali->pairCutPairClusterTC[i/2*4+i%2].push_back(treeCali->pairCutClusterTC[i/2*4+i%2][m]);
          treeCali->pairCutPairClusterTC[i/2*4+i%2+2].push_back(treeCali->pairCutClusterTC[i/2*4+i%2+2][n]);
          treeCali->pairCutPairClusterNumSt[i/2*4+i%2].push_back(treeCali->pairCutClusterNumSt[i/2*4+i%2][m]);
          treeCali->pairCutPairClusterNumSt[i/2*4+i%2+2].push_back(treeCali->pairCutClusterNumSt[i/2*4+i%2+2][n]);
          treeCali->pairCutPairClusterStNum[i/2*4+i%2].push_back(treeCali->pairCutClusterStNum[i/2*4+i%2][m]);
          treeCali->pairCutPairClusterStNum[i/2*4+i%2+2].push_back(treeCali->pairCutClusterStNum[i/2*4+i%2+2][n]);
          treeCali->pairCutPairClusterStCharge[i/2*4+i%2].push_back(treeCali->pairCutClusterStCharge[i/2*4+i%2][m]);
          treeCali->pairCutPairClusterStCharge[i/2*4+i%2+2].push_back(treeCali->pairCutClusterStCharge[i/2*4+i%2+2][n]);
          treeCali->pairCutPairPosCentroid[i/2*4+i%2].push_back(event[i/2*4+i%2]->getEachClusterPositionCentroid(treeCali->pairCutClusterStNum[i/2*4+i%2][m],treeCali->pairCutClusterStCharge[i/2*4+i%2][m]));
          treeCali->pairCutPairPosCentroid[i/2*4+i%2+2].push_back(event[i/2*4+i%2+2]->getEachClusterPositionCentroid(treeCali->pairCutClusterStNum[i/2*4+i%2+2][n],treeCali->pairCutClusterStCharge[i/2*4+i%2+2][n]));
          treeCali->pairCutPairPosErrCentroid[i/2*4+i%2].push_back(event[i/2*4+i%2]->getEachClusterPositionErrorCentroid(treeCali->pairCutClusterStCharge[i/2*4+i%2][m]));
          treeCali->pairCutPairPosErrCentroid[i/2*4+i%2+2].push_back(event[i/2*4+i%2+2]->getEachClusterPositionErrorCentroid(treeCali->pairCutClusterStCharge[i/2*4+i%2+2][n]));
        }
      }
    }
  treeCali->pairCutNumPair[i/2*4+i%2]=treeCali->pairCutPairClusterNum[i/2*4+i%2].size();
  treeCali->pairCutNumPair[i/2*4+i%2+2]=treeCali->pairCutPairClusterNum[i/2*4+i%2+2].size();
  }
  if(treeCali->pairCutNumPair[0]>0 && treeCali->pairCutNumPair[4]>0 && treeCali->pairCutNumPair[8]>0) treeCali->pairCutLeftFlagSecond=1;
  if(treeCali->pairCutNumPair[1]>0 && treeCali->pairCutNumPair[5]>0 && treeCali->pairCutNumPair[9]>0) treeCali->pairCutRightFlagSecond=1;

  ////Specified variables for threshold, total-charge and pair cuts
  for(int i=0;i<12;i++){
    treeCali->thTCPairCutNumCluster[i]=thTCCutClusterStNum[i].size();
    treeCali->thTCPairCutClusterStNum[i]=thTCCutClusterStNum[i];
    treeCali->thTCPairCutClusterStCharge[i]=thTCCutClusterStCharge[i];
    for(UInt_t j=0;j<thTCCutClusterStNum[i].size();j++){
      treeCali->thTCPairCutClusterTC[i].push_back(thTCCutClusterTC[i][j]);
      treeCali->thTCPairCutClusterNumSt[i].push_back(thTCCutClusterStNum[i].size());
      treeCali->thTCPairCutClusterPosCentroid[i].push_back(event[i]->getEachClusterPositionCentroid(treeCali->thTCPairCutClusterStNum[i][j],treeCali->thTCPairCutClusterStCharge[i][j]));
      treeCali->thTCPairCutClusterPosErrCentroid[i].push_back(event[i]->getEachClusterPositionErrorCentroid(treeCali->thTCPairCutClusterStCharge[i][j]));
    }
  }
  if(treeCali->thTCPairCutNumCluster[0]>0 && treeCali->thTCPairCutNumCluster[2]>0 && treeCali->thTCPairCutNumCluster[4]>0 && treeCali->thTCPairCutNumCluster[6]>0 && treeCali->thTCPairCutNumCluster[8]>0 && treeCali->thTCPairCutNumCluster[10]>0) treeCali->thTCPairCutLeftFlagFirst=1;
  if(treeCali->thTCPairCutNumCluster[1]>0 && treeCali->thTCPairCutNumCluster[3]>0 && treeCali->thTCPairCutNumCluster[5]>0 && treeCali->thTCPairCutNumCluster[7]>0 && treeCali->thTCPairCutNumCluster[9]>0 && treeCali->thTCPairCutNumCluster[11]>0) treeCali->thTCPairCutRightFlagFirst=1;

  for(int i=0;i<6;i++){
    for(UInt_t m=0;m<treeCali->thTCPairCutNumCluster[i/2*4+i%2];m++){
      for(UInt_t n=0;n<treeCali->thTCPairCutNumCluster[i/2*4+i%2+2];n++){
        if(treeCali->thTCPairCutClusterTC[i/2*4+i%2][m]-treeCali->thTCPairCutClusterTC[i/2*4+i%2+2][n]>cutFuncLow[i]->Eval(treeCali->thTCPairCutClusterTC[i/2*4+i%2][m]+treeCali->thTCPairCutClusterTC[i/2*4+i%2+2][n]) && treeCali->thTCPairCutClusterTC[i/2*4+i%2][m]-treeCali->thTCPairCutClusterTC[i/2*4+i%2+2][n]<cutFuncUp[i]->Eval(treeCali->thTCPairCutClusterTC[i/2*4+i%2][m]+treeCali->thTCPairCutClusterTC[i/2*4+i%2+2][n])){
          treeCali->thTCPairCutPairClusterNum[i/2*4+i%2].push_back(m);
          treeCali->thTCPairCutPairClusterNum[i/2*4+i%2+2].push_back(n);
          treeCali->thTCPairCutPairClusterTC[i/2*4+i%2].push_back(treeCali->thTCPairCutClusterTC[i/2*4+i%2][m]);
          treeCali->thTCPairCutPairClusterTC[i/2*4+i%2+2].push_back(treeCali->thTCPairCutClusterTC[i/2*4+i%2+2][n]);
          treeCali->thTCPairCutPairClusterNumSt[i/2*4+i%2].push_back(treeCali->thTCPairCutClusterNumSt[i/2*4+i%2][m]);
          treeCali->thTCPairCutPairClusterNumSt[i/2*4+i%2+2].push_back(treeCali->thTCPairCutClusterNumSt[i/2*4+i%2+2][n]);
          treeCali->thTCPairCutPairClusterStNum[i/2*4+i%2].push_back(treeCali->thTCPairCutClusterStNum[i/2*4+i%2][m]);
          treeCali->thTCPairCutPairClusterStNum[i/2*4+i%2+2].push_back(treeCali->thTCPairCutClusterStNum[i/2*4+i%2+2][n]);
          treeCali->thTCPairCutPairClusterStCharge[i/2*4+i%2].push_back(treeCali->thTCPairCutClusterStCharge[i/2*4+i%2][m]);
          treeCali->thTCPairCutPairClusterStCharge[i/2*4+i%2+2].push_back(treeCali->thTCPairCutClusterStCharge[i/2*4+i%2+2][n]);
          treeCali->thTCPairCutPairPosCentroid[i/2*4+i%2].push_back(event[i/2*4+i%2]->getEachClusterPositionCentroid(treeCali->thTCPairCutClusterStNum[i/2*4+i%2][m],treeCali->thTCPairCutClusterStCharge[i/2*4+i%2][m]));
          treeCali->thTCPairCutPairPosCentroid[i/2*4+i%2+2].push_back(event[i/2*4+i%2+2]->getEachClusterPositionCentroid(treeCali->thTCPairCutClusterStNum[i/2*4+i%2+2][n],treeCali->thTCPairCutClusterStCharge[i/2*4+i%2+2][n]));
          treeCali->thTCPairCutPairPosErrCentroid[i/2*4+i%2].push_back(event[i/2*4+i%2]->getEachClusterPositionErrorCentroid(treeCali->thTCPairCutClusterStCharge[i/2*4+i%2][m]));
          treeCali->thTCPairCutPairPosErrCentroid[i/2*4+i%2+2].push_back(event[i/2*4+i%2+2]->getEachClusterPositionErrorCentroid(treeCali->thTCPairCutClusterStCharge[i/2*4+i%2+2][n]));
        }
      }
    }
  treeCali->thTCPairCutNumPair[i/2*4+i%2]=treeCali->thTCPairCutPairClusterNum[i/2*4+i%2].size();
  treeCali->thTCPairCutNumPair[i/2*4+i%2+2]=treeCali->thTCPairCutPairClusterNum[i/2*4+i%2+2].size();
  }
  if(treeCali->thTCPairCutNumPair[0]>0 && treeCali->thTCPairCutNumPair[4]>0 && treeCali->thTCPairCutNumPair[8]>0) treeCali->thTCPairCutLeftFlagSecond=1;
  if(treeCali->thTCPairCutNumPair[1]>0 && treeCali->thTCPairCutNumPair[5]>0 && treeCali->thTCPairCutNumPair[9]>0) treeCali->thTCPairCutRightFlagSecond=1;

  return 0;
}

Long_t Det_MWPC::done_extractCluster(){
  TCanvas *canv=new TCanvas("test","test",1000,600);
  gPad->SetLogz();
  test->Draw("colz");  
  canv->SaveAs("test.pdf");
  return 0;
};

