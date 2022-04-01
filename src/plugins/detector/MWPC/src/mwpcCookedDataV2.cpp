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

//static string inputFilePath(getenv("mwpc_inputFiles"));
//static string inputFilePath="/data/trek/E36/gautam/crnRoot/cooker/inputFiles/MWPC";


//static string envRunNum(getenv("envRunNum"));
//static string inputFile=inputFilePath+"/gapNumTable/run"+envRunNum.c_str();
//static string rootFileDeadStrip=inputFilePath+"/MwpcAdjacentDeadStripCutFunction.root";
//static string rootFileDeadStrip="/data/trek/E36/gautam/crnRoot/cooker/inputFiles/MWPC/MwpcAdjacentDeadStripCutFunction.root";
//static TFile *rootDeadStrip=new TFile(rootFileDeadStrip.c_str());
//static TFile *rootDeadStrip=new TFile("/data/trek/E36/gautam/crnRoot/cooker/inputFiles/MWPC/MwpcAdjacentDeadStripCutFunction.root");
static TF1 *cutFunc0313[6],*cutFunc2429[6];

Long_t Det_MWPC::histos_mwpcCookedDataV2(){
  ostringstream sstr_name,sstr_title;
//cout<<"  check:  "<<endl;
//cout<<" INPUTU:  "<<inputFilePath<<endl;

  return 0;
}

static vector <double> gain[12][6];
Long_t Det_MWPC::startup_mwpcCookedDataV2(){
  getBranchObject("RawMwpcInfo",(TObject **) &treeRaw);
  treeCookV2=0;
  out->Branch("CookMwpcInfoV2",&treeCookV2,160000,0);
  gStyle->SetOptStat(0000);
  gStyle->SetOptFit(1111);
TFile *rootDeadStrip=new TFile("/data/trek/E36/gautam/crnRoot/cooker/inputFiles/MWPC/MwpcAdjacentDeadStripCutFunction.root");

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
Long_t Det_MWPC::process_mwpcCookedDataV2(){

  treeCookV2->init(); // Init data memebers
  treeCookV2->run=treeRaw->run;
  treeCookV2->event=treeRaw->event;

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



  //After gap-dependent cable swapping and gain calibration, positions are determined
  for(int j=0;j<6;j++){

    //C2XL
    for(UInt_t k=0;k<clusterStNum[0].size();k++){
      stCharge.clear();
      for(UInt_t t=0;t<clusterStNum[0][k].size();t++){
        stCharge.push_back(clusterStCharge[0][k][t]*gain[0][j][clusterStNum[0][k][t]]);
      }
      double posChargeRatio=event[0]->getEachClusterPositionChargeRatio(clusterStNum[0][k],stCharge);
      double posChargeRatioLocal=posChargeRatio*10-280;
      double c2Z=posChargeRatioLocal*sqrt(1-c2XCorrPar[j+6]*c2XCorrPar[j+6])+c2ZCenter[j+6];
      double c2X=posChargeRatioLocal*c2XCorrPar[j+6]+c2XCenter[j+6];
      treeCookV2->posC2Z[j+6].push_back(c2Z);
      treeCookV2->posC2X[j+6].push_back(c2X);
   }

    //C2YL
    for(UInt_t k=0;k<clusterStNum[2].size();k++){
      stCharge.clear();
      for(UInt_t t=0;t<clusterStNum[2][k].size();t++){
        stCharge.push_back(clusterStCharge[2][k][t]*gain[2][j][clusterStNum[2][k][t]]);
      }
      double posChargeRatio=16-event[2]->getEachClusterPositionChargeRatio(clusterStNum[2][k],stCharge);
      double posChargeRatioLocal=posChargeRatio*10-80;
      double c2Y=posChargeRatioLocal+c2YCenter[j+6];
      treeCookV2->posC2Y[j+6].push_back(c2Y);
   }  

    //C2XR
    for(UInt_t k=0;k<clusterStNum[1].size();k++){
      stCharge.clear();
      for(UInt_t t=0;t<clusterStNum[1][k].size();t++){
        stCharge.push_back(clusterStCharge[1][k][t]*gain[1][j][clusterStNum[1][k][t]]);
      }
      double posChargeRatio=event[1]->getEachClusterPositionChargeRatio(clusterStNum[1][k],stCharge);
      double posChargeRatioLocal=posChargeRatio*10-280;
      double c2Z=posChargeRatioLocal*sqrt(1-c2XCorrPar[j]*c2XCorrPar[j])+c2ZCenter[j];
      double c2X=posChargeRatioLocal*c2XCorrPar[j]+c2XCenter[j];
      treeCookV2->posC2Z[j].push_back(c2Z);
      treeCookV2->posC2X[j].push_back(c2X);
   }

    //C2YR
    for(UInt_t k=0;k<clusterStNum[3].size();k++){
      stCharge.clear();
      for(UInt_t t=0;t<clusterStNum[3][k].size();t++){
        stCharge.push_back(clusterStCharge[3][k][t]*gain[3][j][clusterStNum[3][k][t]]);
      }
      double posChargeRatio=16-event[3]->getEachClusterPositionChargeRatio(clusterStNum[3][k],stCharge);
      double posChargeRatioLocal=posChargeRatio*10-80;
      double c2Y=posChargeRatioLocal+c2YCenter[j];
      treeCookV2->posC2Y[j].push_back(c2Y);
   }


    //C3XL
    if(j==0){
      for(UInt_t k=0;k<clusterStNumC3XLGap7Swap.size();k++){
        stCharge.clear();
        for(UInt_t t=0;t<clusterStNumC3XLGap7Swap[k].size();t++){
          stCharge.push_back(clusterStChargeC3XLGap7Swap[k][t]*gain[4][0][clusterStNumC3XLGap7Swap[k][t]]);
        }
        double posChargeRatio=event[4]->getEachClusterPositionChargeRatio(clusterStNumC3XLGap7Swap[k],stCharge);
        double posChargeRatioLocal=posChargeRatio*10-320;
        double c3X=posChargeRatioLocal+c3XCenter[6];
        treeCookV2->posC3X[6].push_back(c3X);
       }
    }
    else{
      for(UInt_t k=0;k<clusterStNum[4].size();k++){
        stCharge.clear();
        for(UInt_t t=0;t<clusterStNum[4][k].size();t++){
          stCharge.push_back(clusterStCharge[4][k][t]*gain[4][j][clusterStNum[4][k][t]]);
        }
        double posChargeRatio=event[4]->getEachClusterPositionChargeRatio(clusterStNum[4][k],stCharge);
        double posChargeRatioLocal=posChargeRatio*10-320;
        double c3X=posChargeRatioLocal+c3XCenter[j+6];
        treeCookV2->posC3X[j+6].push_back(c3X);
      }
    }

    //C3YL
    for(UInt_t k=0;k<clusterStNum[6].size();k++){
      stCharge.clear();
      for(UInt_t t=0;t<clusterStNum[6][k].size();t++){
        stCharge.push_back(clusterStCharge[6][k][t]*gain[6][j][clusterStNum[6][k][t]]);
      }
      double posChargeRatio=event[6]->getEachClusterPositionChargeRatio(clusterStNum[6][k],stCharge);
      double posChargeRatioLocal=posChargeRatio*12.5-100;
      double c3Y=posChargeRatioLocal+c3YCenter[j+6];
      treeCookV2->posC3Y[j+6].push_back(c3Y);
    }

    treeCookV2->posC3Z[j+6].push_back(c3ZCenter[j+6]);

    //C3XR
    for(UInt_t k=0;k<clusterStNum[5].size();k++){
      stCharge.clear();
      for(UInt_t t=0;t<clusterStNum[5][k].size();t++){
        stCharge.push_back(clusterStCharge[5][k][t]*gain[5][j][clusterStNum[5][k][t]]);
      }
      double posChargeRatio=event[5]->getEachClusterPositionChargeRatio(clusterStNum[5][k],stCharge);
      double posChargeRatioLocal=posChargeRatio*10-320;
      double c3X=posChargeRatioLocal+c3XCenter[j];
      treeCookV2->posC3X[j].push_back(c3X);
    }
  
    //C3YR
    for(UInt_t k=0;k<clusterStNum[7].size();k++){
      stCharge.clear();
      for(UInt_t t=0;t<clusterStNum[7][k].size();t++){
        stCharge.push_back(clusterStCharge[7][k][t]*gain[7][j][clusterStNum[7][k][t]]);
      }
      double posChargeRatio=event[7]->getEachClusterPositionChargeRatio(clusterStNum[7][k],stCharge);
      double posChargeRatioLocal=posChargeRatio*12.5-100;
      double c3Y=posChargeRatioLocal+c3YCenter[j];
      treeCookV2->posC3Y[j].push_back(c3Y);
    }

    treeCookV2->posC3Z[j].push_back(c3ZCenter[j]);

    //C4XL
    if(j==0){
      for(UInt_t k=0;k<clusterStNumC4XLGap7Swap.size();k++){
        stCharge.clear();
        for(UInt_t t=0;t<clusterStNumC4XLGap7Swap[k].size();t++){
          stCharge.push_back(clusterStChargeC4XLGap7Swap[k][t]*gain[8][0][clusterStNumC4XLGap7Swap[k][t]]);
        }
        double posChargeRatio=event[8]->getEachClusterPositionChargeRatio(clusterStNumC4XLGap7Swap[k],stCharge);
        double posChargeRatioLocal=posChargeRatio*10-360;
        double c4X=posChargeRatioLocal+c4XCenter[6];
        treeCookV2->posC4X[6].push_back(c4X);
       }
    }
    else if(j==4){
      for(UInt_t k=0;k<clusterStNumC4XLGap11Swap.size();k++){
        stCharge.clear();
        for(UInt_t t=0;t<clusterStNumC4XLGap11Swap[k].size();t++){
          stCharge.push_back(clusterStChargeC4XLGap11Swap[k][t]*gain[8][4][clusterStNumC4XLGap11Swap[k][t]]);
        }
        double posChargeRatio=event[8]->getEachClusterPositionChargeRatio(clusterStNumC4XLGap11Swap[k],stCharge);
        double posChargeRatioLocal=posChargeRatio*10-360;
        double c4X=posChargeRatioLocal+c4XCenter[10];
        treeCookV2->posC4X[10].push_back(c4X);
       }
    }
    else{
      for(UInt_t k=0;k<clusterStNum[8].size();k++){
        stCharge.clear();
        for(UInt_t t=0;t<clusterStNum[8][k].size();t++){
          stCharge.push_back(clusterStCharge[8][k][t]*gain[8][j][clusterStNum[8][k][t]]);
        }
        double posChargeRatio=event[8]->getEachClusterPositionChargeRatio(clusterStNum[8][k],stCharge);
        double posChargeRatioLocal=posChargeRatio*10-360;
        double c4X=posChargeRatioLocal+c4XCenter[j+6];
        treeCookV2->posC4X[j+6].push_back(c4X);
      }
    }

    //C4YL
    for(UInt_t k=0;k<clusterStNum[10].size();k++){
      stCharge.clear();
      for(UInt_t t=0;t<clusterStNum[10][k].size();t++){
        stCharge.push_back(clusterStCharge[10][k][t]*gain[10][j][clusterStNum[10][k][t]]);
      }
      double posChargeRatio=16-event[10]->getEachClusterPositionChargeRatio(clusterStNum[10][k],stCharge);
      double posChargeRatioLocal=posChargeRatio*12.5-100;
      double c4Y=posChargeRatioLocal+c4YCenter[j+6];
      treeCookV2->posC4Y[j+6].push_back(c4Y);
    }

    treeCookV2->posC4Z[j+6].push_back(c4ZCenter[j+6]);

    //C4XR
    for(UInt_t k=0;k<clusterStNum[9].size();k++){
      stCharge.clear();
      for(UInt_t t=0;t<clusterStNum[9][k].size();t++){
        stCharge.push_back(clusterStCharge[9][k][t]*gain[9][j][clusterStNum[9][k][t]]);
      }
      double posChargeRatio=event[9]->getEachClusterPositionChargeRatio(clusterStNum[9][k],stCharge);
      double posChargeRatioLocal=posChargeRatio*10-360;
      double c4X=posChargeRatioLocal+c4XCenter[j];
      treeCookV2->posC4X[j].push_back(c4X);
    }

    //C4YR
    for(UInt_t k=0;k<clusterStNum[11].size();k++){
      stCharge.clear();
      for(UInt_t t=0;t<clusterStNum[11][k].size();t++){
        stCharge.push_back(clusterStCharge[11][k][t]*gain[11][j][clusterStNum[11][k][t]]);
      }
      double posChargeRatio=16-event[11]->getEachClusterPositionChargeRatio(clusterStNum[11][k],stCharge);
      double posChargeRatioLocal=posChargeRatio*12.5-100;
      double c4Y=posChargeRatioLocal+c4YCenter[j];
      treeCookV2->posC4Y[j].push_back(c4Y);
    }

    treeCookV2->posC4Z[j].push_back(c4ZCenter[j]);



    //Dead strip
    if(clusterStNum[9].size()==1 && clusterStNum[11].size()==1){
      if(event[9]->getEachClusterStNumWithHighestADC(clusterStNum[9][0],clusterStCharge[9][0])==15 && cutFunc0313[j]->Eval(clusterTC[9][0]+clusterTC[11][0])>clusterTC[9][0]-clusterTC[11][0]) treeCookV2->flagDeadStripClusterC4XR[j].push_back(1);
      else if((event[9]->getEachClusterStNumWithHighestADC(clusterStNum[9][0],clusterStCharge[9][0])==17 || event[9]->getEachClusterStNumWithHighestADC(clusterStNum[9][0],clusterStCharge[9][0])==20 || event[9]->getEachClusterStNumWithHighestADC(clusterStNum[9][0],clusterStCharge[9][0])==22) && cutFunc2429[j]->Eval(clusterTC[9][0]+clusterTC[11][0])>clusterTC[9][0]-clusterTC[11][0]) treeCookV2->flagDeadStripClusterC4XR[j].push_back(1);
      else treeCookV2->flagDeadStripClusterC4XR[j].push_back(0);
    }
    else{
      for(UInt_t k=0;k<clusterStNum[9].size();k++){
        if(event[9]->getEachClusterStNumWithHighestADC(clusterStNum[9][k],clusterStCharge[9][k])==15 || event[9]->getEachClusterStNumWithHighestADC(clusterStNum[9][k],clusterStCharge[9][k])==17 || event[9]->getEachClusterStNumWithHighestADC(clusterStNum[9][k],clusterStCharge[9][k])==20 || event[9]->getEachClusterStNumWithHighestADC(clusterStNum[9][k],clusterStCharge[9][k])==22) treeCookV2->flagDeadStripClusterC4XR[j].push_back(1);
        else treeCookV2->flagDeadStripClusterC4XR[j].push_back(0);
      }
    }

  }
  return 0;
}

Long_t Det_MWPC::done_mwpcCookedDataV2(){

  return 0;
};

