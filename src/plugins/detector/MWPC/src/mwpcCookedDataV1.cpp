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
#include "MwpcGainCali.h"
#include "mwpcLocalGlobal.h"
using namespace std;

static string inputFilePath(getenv("mwpc_inputFiles"));
//static string envRunNum(getenv("envRunNum"));
//static string inputFile=inputFilePath+"/gapNumTable/run"+envRunNum.c_str();
static string inputFile=inputFilePath+"/gapNumTable/run3994";
static ifstream input(inputFile.c_str());
static int runNum, eventNum, gapNumTOF1, gapNumTOF2;
static string rootFileDeadStrip=inputFilePath+"/MwpcAdjacentDeadStripCutFunction.root";
static TFile *rootDeadStrip=new TFile(rootFileDeadStrip.c_str());
static TF1 *cutFunc0313[6],*cutFunc2429[6];

Long_t Det_MWPC::histos_mwpcCookedDataV1(){
  ostringstream sstr_name,sstr_title;
  return 0;
}

static vector <double> gain[12][6];
Long_t Det_MWPC::startup_mwpcCookedDataV1(){
  getBranchObject("RawMwpcInfo",(TObject **) &treeRaw);
  treeCookV1=0;
  out->Branch("CookMwpcInfoV1",&treeCookV1,160000,0);
  gStyle->SetOptStat(0000);
  gStyle->SetOptFit(1111);

  for(int i=0;i<6;i++){
    for(int j=0;j<16;j++){
      gain[2][i].push_back(mwpcGainWeightC2YL[i][j]);
      gain[3][i].push_back(mwpcGainWeightC2YR[i][j]);
      gain[6][i].push_back(mwpcGainWeightC3YL[i][j]);
      gain[7][i].push_back(mwpcGainWeightC3YR[i][j]);
      gain[10][i].push_back(mwpcGainWeightC4YL[i][j]);
      gain[11][i].push_back(mwpcGainWeightC4YR[i][j]);
    }
  }

  for(int i=0;i<6;i++){
    for(int j=0;j<56;j++){
      gain[0][i].push_back(mwpcGainWeightC2XL[i][j]);
      gain[1][i].push_back(mwpcGainWeightC2XR[i][j]);
    }
  }

  for(int i=0;i<6;i++){
    for(int j=0;j<64;j++){
      gain[4][i].push_back(mwpcGainWeightC3XL[i][j]);
      gain[5][i].push_back(mwpcGainWeightC3XR[i][j]);
    }
  }

  for(int i=0;i<6;i++){
    for(int j=0;j<72;j++){
      gain[8][i].push_back(mwpcGainWeightC4XL[i][j]);
      gain[9][i].push_back(mwpcGainWeightC4XR[i][j]);
    }
  }

  ostringstream sstr;
  for(int i=0;i<6;i++){
    sstr<<"St0313_Gap"<<i+1<<"_low";
    cutFunc0313[i]=(TF1*)rootDeadStrip->Get(sstr.str().c_str());
    sstr.str("");
    sstr<<"St2429_Gap"<<i+1<<"_low";
    cutFunc2429[i]=(TF1*)rootDeadStrip->Get(sstr.str().c_str());
    sstr.str("");
  }


  return 0;
};

static const double* pedValue[12]={pedValueC2XL,pedValueC2XR,pedValueC2YL,pedValueC2YR,pedValueC3XL,pedValueC3XR,pedValueC3YL,pedValueC3YR,pedValueC4XL,pedValueC4XR,pedValueC4YL,pedValueC4YR};
static const double* pedSigma[12]={pedSigmaC2XL,pedSigmaC2XR,pedSigmaC2YL,pedSigmaC2YR,pedSigmaC3XL,pedSigmaC3XR,pedSigmaC3YL,pedSigmaC3YR,pedSigmaC4XL,pedSigmaC4XR,pedSigmaC4YL,pedSigmaC4YR};
static const double* totalChargeCut[12]={totChargeCutC2XL,totChargeCutC2XR,totChargeCutC2YL,totChargeCutC2YR,totChargeCutC3XL,totChargeCutC3XR,totChargeCutC3YL,totChargeCutC3YR,totChargeCutC4XL,totChargeCutC4XR,totChargeCutC4YL,totChargeCutC4YR};

static vector<double> stCharge;
static double totalCharge;
Long_t Det_MWPC::process_mwpcCookedDataV1(){
  input>>runNum>>eventNum>>gapNumTOF1>>gapNumTOF2;
  if(runNum!=treeRaw->run || eventNum!=treeRaw->event){
    cout<<"Error: Run number or event number are not consistent between two files!"<<endl;
    cout<<"For MWPC Raw file: "<<"Run number - "<<treeCookV1->run<<"and  Event number - "<<treeCookV1->event<<endl;
    cout<<"For TRIUMF file: "<<"Run number - "<<runNum<<"and  Event number - "<<eventNum<<endl;
  }

  treeCookV1->init(); // Init data memebers
  treeCookV1->run=treeRaw->run;
  treeCookV1->event=treeRaw->event;

  int* data[12]={treeRaw->ADC_C2X_L,treeRaw->ADC_C2X_R,treeRaw->ADC_C2Y_L,treeRaw->ADC_C2Y_R,treeRaw->ADC_C3X_L,treeRaw->ADC_C3X_R,treeRaw->ADC_C3Y_L,treeRaw->ADC_C3Y_R,treeRaw->ADC_C4X_L,treeRaw->ADC_C4X_R,treeRaw->ADC_C4Y_L,treeRaw->ADC_C4Y_R};
  MwpcEvent *event[12];
  vector<vector<int> > clusterStNum[12];
  vector<vector<double> > clusterStCharge[12];
  double *clusterTC[12];
  for(int i=0;i<12;i++){
    event[i]=new MwpcEvent(((i%4)/2*(-5)+(i/4)*(!((i%4)/2))+7)*8,!((i/4==0)&&((i%4/2))==0));
    clusterStNum[i]=event[i]->clusterStripNumber(event[i],data[i],pedValue[i],pedSigma[i],1,totalChargeCut[i],1,singMiniSigmaClusterSecondApproach[i],1,doubMiniSigmaXClusterSecondApproach[i],doubMiniTotalChargeXCluster[i]);
    clusterStCharge[i]=event[i]->clusterStripCharge(event[i],data[i],pedValue[i],pedSigma[i],1,totalChargeCut[i],1,singMiniSigmaClusterSecondApproach[i],1,doubMiniSigmaXClusterSecondApproach[i],doubMiniTotalChargeXCluster[i]);
    clusterTC[i]=event[i]->getClusterTotCharge(event[i],data[i],pedValue[i],pedSigma[i],1,totalChargeCut[i],1,singMiniSigmaClusterSecondApproach[i],1,doubMiniSigmaXClusterSecondApproach[i],doubMiniTotalChargeXCluster[i]);
  }

  //Cable Swapping: For gap 7, C3X3<->C4X3 & C3X4<->C4X4; For gap 11, C4X3<->C4X4
  int dataC3XLGap7Swap[64],dataC4XLGap7Swap[72],dataC4XLGap11Swap[72];
  double pedValueC3XLGap7Swap[64], pedValueC4XLGap7Swap[72], pedValueC4XLGap11Swap[72];
  double pedSigmaC3XLGap7Swap[64], pedSigmaC4XLGap7Swap[72], pedSigmaC4XLGap11Swap[72];
  for(int i=0;i<32;i++){
    dataC3XLGap7Swap[i]=data[4][i];
    dataC4XLGap7Swap[i]=data[8][i];
    dataC4XLGap11Swap[i]=data[8][i];
    pedValueC3XLGap7Swap[i]=pedValue[4][i];
    pedValueC4XLGap7Swap[i]=pedValue[8][i];
    pedValueC4XLGap11Swap[i]=pedValue[8][i];
    pedSigmaC3XLGap7Swap[i]=pedSigma[4][i];
    pedSigmaC4XLGap7Swap[i]=pedSigma[8][i];
    pedSigmaC4XLGap11Swap[i]=pedSigma[8][i];
  }
  for(int i=32;i<64;i++){
    dataC3XLGap7Swap[i]=data[8][i];
    dataC4XLGap7Swap[i]=data[4][i];
    pedValueC3XLGap7Swap[i]=pedValue[8][i];
    pedValueC4XLGap7Swap[i]=pedValue[4][i];
    pedSigmaC3XLGap7Swap[i]=pedSigma[8][i];
    pedSigmaC4XLGap7Swap[i]=pedSigma[4][i];
  }
  for(int i=32;i<48;i++){
    dataC4XLGap11Swap[i]=data[8][i+16];
    pedValueC4XLGap11Swap[i]=pedValue[8][i+16];
    pedSigmaC4XLGap11Swap[i]=pedSigma[8][i+16];
  }
  for(int i=48;i<64;i++){
    dataC4XLGap11Swap[i]=data[8][i-16];
    pedValueC4XLGap11Swap[i]=pedValue[8][i-16];
    pedSigmaC4XLGap11Swap[i]=pedSigma[8][i-16];
  }
  for(int i=64;i<72;i++){
    dataC4XLGap7Swap[i]=data[8][i];
    dataC4XLGap11Swap[i]=data[8][i];
    pedValueC4XLGap7Swap[i]=pedValue[8][i];
    pedValueC4XLGap11Swap[i]=pedValue[8][i];
    pedSigmaC4XLGap7Swap[i]=pedSigma[8][i];
    pedSigmaC4XLGap11Swap[i]=pedSigma[8][i];
  }

  vector<vector<int> > clusterStNumC3XLGap7Swap, clusterStNumC4XLGap7Swap, clusterStNumC4XLGap11Swap;
  vector<vector<double> > clusterStChargeC3XLGap7Swap, clusterStChargeC4XLGap7Swap, clusterStChargeC4XLGap11Swap;
  double *clusterTCC3XLGap7Swap, *clusterTCC4XLGap7Swap, *clusterTCC4XLGap11Swap;
  clusterStNumC3XLGap7Swap=event[4]->clusterStripNumber(event[4],dataC3XLGap7Swap,pedValueC3XLGap7Swap,pedSigmaC3XLGap7Swap,1,totalChargeCut[4],1,3.3,1,2.65,180);
  clusterStChargeC3XLGap7Swap=event[4]->clusterStripCharge(event[4],dataC3XLGap7Swap,pedValueC3XLGap7Swap,pedSigmaC3XLGap7Swap,1,totalChargeCut[4],1,3.3,1,2.65,180);
  clusterTCC3XLGap7Swap=event[4]->getClusterTotCharge(event[4],dataC3XLGap7Swap,pedValueC3XLGap7Swap,pedSigmaC3XLGap7Swap,1,totalChargeCut[4],1,3.3,1,2.65,180);
  clusterStNumC4XLGap7Swap=event[8]->clusterStripNumber(event[8],dataC4XLGap7Swap,pedValueC4XLGap7Swap,pedSigmaC4XLGap7Swap,1,totalChargeCut[8],1,3.3,1,2.65,180);
  clusterStChargeC4XLGap7Swap=event[8]->clusterStripCharge(event[8],dataC4XLGap7Swap,pedValueC4XLGap7Swap,pedSigmaC4XLGap7Swap,1,totalChargeCut[8],1,3.3,1,2.65,180);
  clusterTCC4XLGap7Swap=event[8]->getClusterTotCharge(event[8],dataC4XLGap7Swap,pedValueC4XLGap7Swap,pedSigmaC4XLGap7Swap,1,totalChargeCut[8],1,3.3,1,2.65,180);
  clusterStNumC4XLGap11Swap=event[8]->clusterStripNumber(event[8],dataC4XLGap11Swap,pedValueC4XLGap11Swap,pedSigmaC4XLGap11Swap,1,totalChargeCut[8],1,3.4,1,2.7,180);
  clusterStChargeC4XLGap11Swap=event[8]->clusterStripCharge(event[8],dataC4XLGap11Swap,pedValueC4XLGap11Swap,pedSigmaC4XLGap11Swap,1,totalChargeCut[8],1,3.4,1,2.7,180);
  clusterTCC4XLGap11Swap=event[8]->getClusterTotCharge(event[8],dataC4XLGap11Swap,pedValueC4XLGap11Swap,pedSigmaC4XLGap11Swap,1,totalChargeCut[8],1,3.4,1,2.7,180);
  //Cable Swapping

  double* adcValuePedSub;
  for(int i=0;i<12;i++){
    adcValuePedSub=event[i]->getAdcValue(data[i],pedValue[i]);
    for(int j=0;j<((i%4)/2*(-5)+(i/4)*(!((i%4)/2))+7)*8;j++){
      treeCookV1->adcValueOriginal[i].push_back(data[i][j]);
      treeCookV1->adcValuePedSub[i].push_back(adcValuePedSub[j]);
    }
  }

 //Cable Swapping
  treeCookV1->numClusterC3XLGap7Swap=clusterStNumC3XLGap7Swap.size();
  treeCookV1->clusterStNumC3XLGap7Swap=clusterStNumC3XLGap7Swap;
  treeCookV1->clusterStChargeC3XLGap7Swap=clusterStChargeC3XLGap7Swap;
  for(UInt_t j=0;j<clusterStNumC3XLGap7Swap.size();j++){
    treeCookV1->clusterTCC3XLGap7Swap.push_back(clusterTCC3XLGap7Swap[j]);
    treeCookV1->clusterNumStC3XLGap7Swap.push_back(clusterStNumC3XLGap7Swap.size());
  }
  treeCookV1->numClusterC4XLGap7Swap=clusterStNumC4XLGap7Swap.size();
  treeCookV1->clusterStNumC4XLGap7Swap=clusterStNumC4XLGap7Swap;
  treeCookV1->clusterStChargeC4XLGap7Swap=clusterStChargeC4XLGap7Swap;
  for(UInt_t j=0;j<clusterStNumC4XLGap7Swap.size();j++){
    treeCookV1->clusterTCC4XLGap7Swap.push_back(clusterTCC4XLGap7Swap[j]);
    treeCookV1->clusterNumStC4XLGap7Swap.push_back(clusterStNumC4XLGap7Swap.size());
  }
  treeCookV1->numClusterC4XLGap11Swap=clusterStNumC4XLGap11Swap.size();
  treeCookV1->clusterStNumC4XLGap11Swap=clusterStNumC4XLGap11Swap;
  treeCookV1->clusterStChargeC4XLGap11Swap=clusterStChargeC4XLGap11Swap;
  for(UInt_t j=0;j<clusterStNumC4XLGap11Swap.size();j++){
    treeCookV1->clusterTCC4XLGap11Swap.push_back(clusterTCC4XLGap11Swap[j]);
    treeCookV1->clusterNumStC4XLGap11Swap.push_back(clusterStNumC4XLGap11Swap.size());
  }
  //Cabel Swapping

  for(int i=0;i<12;i++){
    treeCookV1->numCluster[i]=clusterStNum[i].size();
    treeCookV1->clusterStNum[i]=clusterStNum[i];
    treeCookV1->clusterStCharge[i]=clusterStCharge[i];
    for(UInt_t j=0;j<clusterStNum[i].size();j++){
      treeCookV1->clusterTC[i].push_back(clusterTC[i][j]);
      treeCookV1->clusterNumSt[i].push_back(clusterStNum[i].size());
    }
  }

  //After gap-dependent cable swapping and gain calibration, positions are determined
  treeCookV1->gapNum=gapNumTOF2;
  for(int i=0;i<12;i++){
    if(gapNumTOF2!=-1){
      if((i==4 || i==6) && gapNumTOF2==7){
        if(treeCookV1->numClusterC3XLGap7Swap==1 && treeCookV1->numCluster[6]==1 && i==4){
          stCharge.clear();
          totalCharge=0;
          for(UInt_t j=0;j<treeCookV1->clusterStNumC3XLGap7Swap[0].size();j++){
            stCharge.push_back(treeCookV1->clusterStChargeC3XLGap7Swap[0][j]*gain[4][0][treeCookV1->clusterStNumC3XLGap7Swap[0][j]]);
            totalCharge+=treeCookV1->clusterStChargeC3XLGap7Swap[0][j]*gain[4][0][treeCookV1->clusterStNumC3XLGap7Swap[0][j]];
          }
          treeCookV1->totalCharge[4]=totalCharge;
          treeCookV1->numSt[4]=treeCookV1->clusterStNumC3XLGap7Swap[0].size();
          treeCookV1->stNum[4]=treeCookV1->clusterStNumC3XLGap7Swap[0];
          treeCookV1->stCharge[4]=stCharge;
          treeCookV1->posCentroid[4]=event[4]->getEachClusterPositionCentroidUpto3St(treeCookV1->clusterStNumC3XLGap7Swap[0],stCharge);
          treeCookV1->posErrCentroid[4]=event[4]->getEachClusterPositionErrorCentroidUpto3St(stCharge);
          treeCookV1->posChargeRatio[4]=event[4]->getEachClusterPositionChargeRatio(treeCookV1->clusterStNumC3XLGap7Swap[0],stCharge);
          treeCookV1->posChargeRatioLocal[4]=treeCookV1->posChargeRatio[4]*10-320;
          treeCookV1->flagLeftRight=-1;
          stCharge.clear();
          totalCharge=0;
          for(UInt_t j=0;j<treeCookV1->clusterStNum[6][0].size();j++){
            stCharge.push_back(treeCookV1->clusterStCharge[6][0][j]*gain[6][0][treeCookV1->clusterStNum[6][0][j]]);
            totalCharge+=treeCookV1->clusterStCharge[6][0][j]*gain[6][0][treeCookV1->clusterStNum[6][0][j]];
          }
          treeCookV1->totalCharge[6]=totalCharge;
          treeCookV1->numSt[6]=treeCookV1->clusterStNum[6][0].size();
          treeCookV1->stNum[6]=treeCookV1->clusterStNum[6][0];
          treeCookV1->stCharge[6]=stCharge;
          treeCookV1->posCentroid[6]=event[6]->getEachClusterPositionCentroidUpto3St(treeCookV1->clusterStNum[6][0],stCharge);
          treeCookV1->posErrCentroid[6]=event[6]->getEachClusterPositionErrorCentroidUpto3St(stCharge);
          treeCookV1->posChargeRatio[6]=event[6]->getEachClusterPositionChargeRatio(treeCookV1->clusterStNum[6][0],stCharge);
          treeCookV1->posChargeRatioLocal[6]=treeCookV1->posChargeRatio[6]*12.5-100;

          treeCookV1->posC3X=treeCookV1->posChargeRatioLocal[4]+c3XCenter[gapNumTOF2-1];
          treeCookV1->posC3Y=treeCookV1->posChargeRatioLocal[6]+c3YCenter[gapNumTOF2-1];
          treeCookV1->posC3Z=c3ZCenter[gapNumTOF2-1];

          treeCookV1->flagLeftRight=-1;
        }
      }
      else if((i==8 || i==10) && gapNumTOF2==7){
        if(treeCookV1->numClusterC4XLGap7Swap==1 && treeCookV1->numCluster[10]==1 && i==8){
          stCharge.clear();
          totalCharge=0;
          for(UInt_t j=0;j<treeCookV1->clusterStNumC4XLGap7Swap[0].size();j++){
            stCharge.push_back(treeCookV1->clusterStChargeC4XLGap7Swap[0][j]*gain[8][0][treeCookV1->clusterStNumC4XLGap7Swap[0][j]]);
            totalCharge+=treeCookV1->clusterStChargeC4XLGap7Swap[0][j]*gain[8][0][treeCookV1->clusterStNumC4XLGap7Swap[0][j]];
          }
          treeCookV1->totalCharge[8]=totalCharge;
          treeCookV1->numSt[8]=treeCookV1->clusterStNumC4XLGap7Swap[0].size();
          treeCookV1->stNum[8]=treeCookV1->clusterStNumC4XLGap7Swap[0];
          treeCookV1->stCharge[8]=stCharge;
          treeCookV1->posCentroid[8]=event[8]->getEachClusterPositionCentroidUpto3St(treeCookV1->clusterStNumC4XLGap7Swap[0],stCharge);
          treeCookV1->posErrCentroid[8]=event[8]->getEachClusterPositionErrorCentroidUpto3St(stCharge);
          treeCookV1->posChargeRatio[8]=event[8]->getEachClusterPositionChargeRatio(treeCookV1->clusterStNumC4XLGap7Swap[0],stCharge);
          treeCookV1->posChargeRatioLocal[8]=treeCookV1->posChargeRatio[8]*10-360;
          treeCookV1->flagLeftRight=-1;
          stCharge.clear();
          totalCharge=0;
          for(UInt_t j=0;j<treeCookV1->clusterStNum[10][0].size();j++){
            stCharge.push_back(treeCookV1->clusterStCharge[10][0][j]*gain[10][0][treeCookV1->clusterStNum[10][0][j]]);
            totalCharge+=treeCookV1->clusterStCharge[10][0][j]*gain[10][0][treeCookV1->clusterStNum[10][0][j]];
          }
          treeCookV1->totalCharge[10]=totalCharge;
          treeCookV1->numSt[10]=treeCookV1->clusterStNum[10][0].size();
          treeCookV1->stNum[10]=treeCookV1->clusterStNum[10][0];
          treeCookV1->stCharge[10]=stCharge;
          treeCookV1->posCentroid[10]=16-event[10]->getEachClusterPositionCentroidUpto3St(treeCookV1->clusterStNum[10][0],stCharge);
          treeCookV1->posErrCentroid[10]=event[10]->getEachClusterPositionErrorCentroidUpto3St(stCharge);
          treeCookV1->posChargeRatio[10]=16-event[10]->getEachClusterPositionChargeRatio(treeCookV1->clusterStNum[10][0],stCharge);
          treeCookV1->posChargeRatioLocal[10]=treeCookV1->posChargeRatio[10]*12.5-100;

          treeCookV1->posC4X=treeCookV1->posChargeRatioLocal[8]+c4XCenter[gapNumTOF2-1];
          treeCookV1->posC4Y=treeCookV1->posChargeRatioLocal[10]+c4YCenter[gapNumTOF2-1];
          treeCookV1->posC4Z=c4ZCenter[gapNumTOF2-1];

          treeCookV1->flagLeftRight=-1;
        }
      }
      else if((i==8 || i==10) && gapNumTOF2==11){
        if(treeCookV1->numClusterC4XLGap11Swap==1 && treeCookV1->numCluster[10]==1 && i==8){
          stCharge.clear();
          totalCharge=0;
          for(UInt_t j=0;j<treeCookV1->clusterStNumC4XLGap11Swap[0].size();j++){
            stCharge.push_back(treeCookV1->clusterStChargeC4XLGap11Swap[0][j]*gain[8][4][treeCookV1->clusterStNumC4XLGap11Swap[0][j]]);
            totalCharge+=treeCookV1->clusterStChargeC4XLGap11Swap[0][j]*gain[8][4][treeCookV1->clusterStNumC4XLGap11Swap[0][j]];
          }
          treeCookV1->totalCharge[8]=totalCharge;
          treeCookV1->numSt[8]=treeCookV1->clusterStNumC4XLGap11Swap[0].size();
          treeCookV1->stNum[8]=treeCookV1->clusterStNumC4XLGap11Swap[0];
          treeCookV1->stCharge[8]=stCharge;
          treeCookV1->posCentroid[8]=event[8]->getEachClusterPositionCentroidUpto3St(treeCookV1->clusterStNumC4XLGap11Swap[0],stCharge);
          treeCookV1->posErrCentroid[8]=event[8]->getEachClusterPositionErrorCentroidUpto3St(stCharge);
          treeCookV1->posChargeRatio[8]=event[8]->getEachClusterPositionChargeRatio(treeCookV1->clusterStNumC4XLGap11Swap[0],stCharge);
          treeCookV1->posChargeRatioLocal[8]=treeCookV1->posChargeRatio[8]*10-360;
          treeCookV1->flagLeftRight=-1;
          stCharge.clear();
          totalCharge=0;
          for(UInt_t j=0;j<treeCookV1->clusterStNum[10][0].size();j++){
            stCharge.push_back(treeCookV1->clusterStCharge[10][0][j]*gain[10][4][treeCookV1->clusterStNum[10][0][j]]);
            totalCharge+=treeCookV1->clusterStCharge[10][0][j]*gain[10][4][treeCookV1->clusterStNum[10][0][j]];
          }
          treeCookV1->totalCharge[10]=totalCharge;
          treeCookV1->numSt[10]=treeCookV1->clusterStNum[10][0].size();
          treeCookV1->stNum[10]=treeCookV1->clusterStNum[10][0];
          treeCookV1->stCharge[10]=stCharge;
          treeCookV1->posCentroid[10]=16-event[10]->getEachClusterPositionCentroidUpto3St(treeCookV1->clusterStNum[10][0],stCharge);
          treeCookV1->posErrCentroid[10]=event[10]->getEachClusterPositionErrorCentroidUpto3St(stCharge);
          treeCookV1->posChargeRatio[10]=16-event[10]->getEachClusterPositionChargeRatio(treeCookV1->clusterStNum[10][0],stCharge);
          treeCookV1->posChargeRatioLocal[10]=treeCookV1->posChargeRatio[10]*12.5-100;

          treeCookV1->posC4X=treeCookV1->posChargeRatioLocal[8]+c4XCenter[gapNumTOF2-1];
          treeCookV1->posC4Y=treeCookV1->posChargeRatioLocal[10]+c4YCenter[gapNumTOF2-1];
          treeCookV1->posC4Z=c4ZCenter[gapNumTOF2-1];

          treeCookV1->flagLeftRight=-1;
        }
      }
      else{
        if(treeCookV1->numCluster[i]==1 && treeCookV1->numCluster[i+2-i%4/2*4]==1){
          if(gapNumTOF2>=(!(i%2))*6+1 && gapNumTOF2<=(!(i%2))*6+6){
            stCharge.clear();
            totalCharge=0;
            for(UInt_t j=0;j<treeCookV1->clusterStNum[i][0].size();j++){
              stCharge.push_back(treeCookV1->clusterStCharge[i][0][j]*gain[i][gapNumTOF2-(!(i%2))*6-1][treeCookV1->clusterStNum[i][0][j]]);
              totalCharge+=treeCookV1->clusterStCharge[i][0][j]*gain[i][gapNumTOF2-(!(i%2))*6-1][treeCookV1->clusterStNum[i][0][j]];
            }
            treeCookV1->totalCharge[i]=totalCharge;
            treeCookV1->numSt[i]=treeCookV1->clusterStNum[i][0].size();
            treeCookV1->stNum[i]=treeCookV1->clusterStNum[i][0];
            treeCookV1->stCharge[i]=stCharge;
            treeCookV1->posErrCentroid[i]=event[i]->getEachClusterPositionErrorCentroidUpto3St(stCharge);
            if(i==6 || i==7){
              treeCookV1->posCentroid[i]=event[i]->getEachClusterPositionCentroidUpto3St(treeCookV1->clusterStNum[i][0],stCharge);
              treeCookV1->posChargeRatio[i]=event[i]->getEachClusterPositionChargeRatio(treeCookV1->clusterStNum[i][0],stCharge);
              treeCookV1->posChargeRatioLocal[i]=treeCookV1->posChargeRatio[i]*12.5-100;

              treeCookV1->posC3Y=treeCookV1->posChargeRatioLocal[i]+c3YCenter[gapNumTOF2-1];
              treeCookV1->posC3Z=c3ZCenter[gapNumTOF2-1];
            } 
            else if(i==10 || i==11){
              treeCookV1->posCentroid[i]=16-event[i]->getEachClusterPositionCentroidUpto3St(treeCookV1->clusterStNum[i][0],stCharge);
              treeCookV1->posChargeRatio[i]=16-event[i]->getEachClusterPositionChargeRatio(treeCookV1->clusterStNum[i][0],stCharge);
              treeCookV1->posChargeRatioLocal[i]=treeCookV1->posChargeRatio[i]*12.5-100;

              treeCookV1->posC4Y=treeCookV1->posChargeRatioLocal[i]+c4YCenter[gapNumTOF2-1];
              treeCookV1->posC4Z=c4ZCenter[gapNumTOF2-1];

            }
            else if(i==4 || i==5){
              treeCookV1->posCentroid[i]=event[i]->getEachClusterPositionCentroidUpto3St(treeCookV1->clusterStNum[i][0],stCharge);
              treeCookV1->posChargeRatio[i]=event[i]->getEachClusterPositionChargeRatio(treeCookV1->clusterStNum[i][0],stCharge);
              treeCookV1->posChargeRatioLocal[i]=treeCookV1->posChargeRatio[i]*10-((i%4)/2*(-5)+(i/4)*(!((i%4)/2))+7)*40;

              treeCookV1->posC3X=treeCookV1->posChargeRatioLocal[i]+c3XCenter[gapNumTOF2-1];
            }
            else if(i==8 || i==9){
              treeCookV1->posCentroid[i]=event[i]->getEachClusterPositionCentroidUpto3St(treeCookV1->clusterStNum[i][0],stCharge);
              treeCookV1->posChargeRatio[i]=event[i]->getEachClusterPositionChargeRatio(treeCookV1->clusterStNum[i][0],stCharge);
              treeCookV1->posChargeRatioLocal[i]=treeCookV1->posChargeRatio[i]*10-((i%4)/2*(-5)+(i/4)*(!((i%4)/2))+7)*40;

              treeCookV1->posC4X=treeCookV1->posChargeRatioLocal[i]+c4XCenter[gapNumTOF2-1];
            }
            else if(i==2 || i==3){
              treeCookV1->posCentroid[i]=16-event[i]->getEachClusterPositionCentroidUpto3St(treeCookV1->clusterStNum[i][0],stCharge);
              treeCookV1->posChargeRatio[i]=16-event[i]->getEachClusterPositionChargeRatio(treeCookV1->clusterStNum[i][0],stCharge);
              treeCookV1->posChargeRatioLocal[i]=treeCookV1->posChargeRatio[i]*10-((i%4)/2*(-5)+(i/4)*(!((i%4)/2))+7)*40;

              treeCookV1->posC2Y=treeCookV1->posChargeRatioLocal[i]+c2YCenter[gapNumTOF2-1];
            }
            else if(i==0 || i==1){
              treeCookV1->posCentroid[i]=event[i]->getEachClusterPositionCentroidUpto3St(treeCookV1->clusterStNum[i][0],stCharge);
              treeCookV1->posChargeRatio[i]=event[i]->getEachClusterPositionChargeRatio(treeCookV1->clusterStNum[i][0],stCharge);
              treeCookV1->posChargeRatioLocal[i]=treeCookV1->posChargeRatio[i]*10-((i%4)/2*(-5)+(i/4)*(!((i%4)/2))+7)*40;
              
              treeCookV1->posC2Z=treeCookV1->posChargeRatioLocal[i]*sqrt(1-c2XCorrPar[gapNumTOF2-1]*c2XCorrPar[gapNumTOF2-1])+c2ZCenter[gapNumTOF2-1];
              treeCookV1->posC2X=treeCookV1->posChargeRatioLocal[i]*c2XCorrPar[gapNumTOF2-1]+c2XCenter[gapNumTOF2-1];
            }
            if(gapNumTOF2>=7) treeCookV1->flagLeftRight=-1;
            else treeCookV1->flagLeftRight=1;
            if(i==9){
              if(event[9]->getEachClusterStNumWithHighestADC(treeCookV1->clusterStNum[9][0],treeCookV1->clusterStCharge[9][0])==15 && cutFunc0313[gapNumTOF2-1]->Eval(treeCookV1->clusterTC[9][0]+treeCookV1->clusterTC[11][0])>treeCookV1->clusterTC[9][0]-treeCookV1->clusterTC[11][0]) treeCookV1->flagDeadStripCluster=1;
              else if((event[9]->getEachClusterStNumWithHighestADC(treeCookV1->clusterStNum[9][0],treeCookV1->clusterStCharge[9][0])==17 || event[9]->getEachClusterStNumWithHighestADC(treeCookV1->clusterStNum[9][0],treeCookV1->clusterStCharge[9][0])==20 || event[9]->getEachClusterStNumWithHighestADC(treeCookV1->clusterStNum[9][0],treeCookV1->clusterStCharge[9][0])==22) && cutFunc2429[gapNumTOF2-1]->Eval(treeCookV1->clusterTC[9][0]+treeCookV1->clusterTC[11][0])>treeCookV1->clusterTC[9][0]-treeCookV1->clusterTC[11][0]) treeCookV1->flagDeadStripCluster=1;
            }
          }
        }
      }
    }
  }

  return 0;
}

Long_t Det_MWPC::done_mwpcCookedDataV1(){

  return 0;
};

