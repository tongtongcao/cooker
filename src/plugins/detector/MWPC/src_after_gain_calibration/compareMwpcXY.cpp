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
#include "MwpcPedestalValue.h"
#include "MwpcClusterThresholdSecondApproach.h"
#include "MwpcClusterThresholdFromTotalCharge.h"
using namespace std;

static  TFile *fileCutFunc=new TFile("/home/caot/trek-e36/cooker/src/plugins/detector/MWPC/src/MwpcDiffVsSumTotalChargeCutFunction.root");
static  TF1 *cutFuncUp[6],*cutFuncLow[6];

static TH1D *histTotalCharge_allcut[12];
static TH1D *histTotalCharge_nosingletcut[12];
static TH2D *histXYNClusNoCut[6];
static TH2D *histXYNClusBeFlagCut[6];
static TH2D *histXYNClusAfFlagCut[6];
static TH1D *histDiffXYTotalChargeX1Y1[6];
static TH1D *histSumXYTotalChargeX1Y1[6];
static TH2D *histDiffVsSumXYTotalChargeX1Y1[6];
static TH2D *histDiffVsSumXYTotalChargeX1Y1OnlyFlagCut[6];
static TH2D *histDiffVsSumXYTotalChargeX1Y1RNoAAll[6];
static TH2D *histDiffVsSumXYTotalChargeX1Y1ANoRAll[6];
static TH2D *histDiffVsSumXYTotalChargeX1Y1ANoRAllNoFlagCut[6];
static TH1D *histX1Y2DiffY2TotalCharge[6]; // X1Y2 case, difference of total charge of 2 clusters in Y
static TH2D *histDiffX1Y2TotalChargeFurtherVsCloser[6];
static TH2D *histDiffVsSumXYTotalChargeX1Y2Closer[6];
static TH2D *histDiffVsSumXYTotalChargeX1Y2Further[6];
static TH2D *histDiffVsSumXYTotalChargeX1Y2[6];
static TH1D *histX2Y1DiffX2TotalCharge[6]; // X2Y1 case, difference of total charge of 2 clusters in X
static TH2D *histDiffX2Y1TotalChargeFurtherVsCloser[6];
static TH2D *histDiffVsSumXYTotalChargeX2Y1Closer[6];
static TH2D *histDiffVsSumXYTotalChargeX2Y1Further[6];
static TH2D *histDiffVsSumXYTotalChargeX2Y1[6];

Long_t Det_MWPC::histos_comparexy(){
  ostringstream sstr;
  for(int i=0;i<6;i++){
    sstr<<"C"<<i/2+2<<(char)((i%2)*6+76)<<"_up";
    cutFuncUp[i]=(TF1*)fileCutFunc->Get(sstr.str().c_str());
    sstr.str("");
    sstr<<"C"<<i/2+2<<(char)((i%2)*6+76)<<"_low";
    cutFuncLow[i]=(TF1*)fileCutFunc->Get(sstr.str().c_str());
    sstr.str("");
  }

  ostringstream sstr_name, sstr_title;
  for(int i=0;i<12;i++){
    sstr_name<<"Distribution of Total charge for all cut for "<<"MwpcC"<<i/4+2<<(char)((i%4)/2+88)<<(char)((i%4%2)*6+76);
    sstr_title<<sstr_name.str()<<";Total Charge;Counts";
    histTotalCharge_allcut[i]=new TH1D(sstr_name.str().c_str(),sstr_title.str().c_str(),100,0,3000);
    sstr_name.str("");
    sstr_title.str("");

    sstr_name<<"Distribution of Total charge for no singlet cut for "<<"MwpcC"<<i/4+2<<(char)((i%4)/2+88)<<(char)((i%4%2)*6+76);
    sstr_title<<sstr_name.str()<<";Total Charge;Counts";
    histTotalCharge_nosingletcut[i]=new TH1D(sstr_name.str().c_str(),sstr_title.str().c_str(),100,0,3000);
    sstr_name.str("");
    sstr_title.str("");
  }

  for(int i=0;i<6;i++){
    sstr_name<<"Comparison of number of clusters between X and Y without any cuts for  "<<"MwpcC"<<i/2+2<<(char)((i%2)*6+76);
    sstr_title<<sstr_name.str()<<";Number of clusters in X;Number of clusters in Y";
    histXYNClusNoCut[i]=new TH2D(sstr_name.str().c_str(),sstr_title.str().c_str(),15,-0.5,14.5,6,-0.5,5.5);
    sstr_name.str("");
    sstr_title.str("");

    sstr_name<<"Comparison of number of clusters between X and Y before flag cut for  "<<"MwpcC"<<i/2+2<<(char)((i%2)*6+76);
    sstr_title<<sstr_name.str()<<";Number of clusters in X;Number of clusters in Y";
    histXYNClusBeFlagCut[i]=new TH2D(sstr_name.str().c_str(),sstr_title.str().c_str(),15,-0.5,14.5,6,-0.5,5.5);
    sstr_name.str("");
    sstr_title.str("");

    sstr_name<<"Comparison of number of clusters between X and Y after flag cut for  "<<"MwpcC"<<i/2+2<<(char)((i%2)*6+76);
    sstr_title<<sstr_name.str()<<";Number of clusters in X;Number of clusters in Y";
    histXYNClusAfFlagCut[i]=new TH2D(sstr_name.str().c_str(),sstr_title.str().c_str(),15,-0.5,14.5,6,-0.5,5.5);
    sstr_name.str("");
    sstr_title.str("");

    sstr_name<<"Difference of total charge between X and Y for X1Y1 for  "<<"MwpcC"<<i/2+2<<(char)((i%2)*6+76);
    sstr_title<<sstr_name.str()<<";TC_{X} - TC_{Y};Counts";
    histDiffXYTotalChargeX1Y1[i]=new TH1D(sstr_name.str().c_str(),sstr_title.str().c_str(),200,-1000,1500);
    sstr_name.str("");
    sstr_title.str("");

    sstr_name<<"Sum of total charge of X and Y for X1Y1 for  "<<"MwpcC"<<i/2+2<<(char)((i%2)*6+76);
    sstr_title<<sstr_name.str()<<";TC_{X} + TC_{Y};Counts";
    histSumXYTotalChargeX1Y1[i]=new TH1D(sstr_name.str().c_str(),sstr_title.str().c_str(),200,0,5000);
    sstr_name.str("");
    sstr_title.str("");

    sstr_name<<"Difference vs Sum for total charge of X and Y for X1Y1 for  "<<"MwpcC"<<i/2+2<<(char)((i%2)*6+76);
    sstr_title<<sstr_name.str()<<";TC_{X} + TC_{Y};TC_{X} - TC_{Y}";
    histDiffVsSumXYTotalChargeX1Y1[i]=new TH2D(sstr_name.str().c_str(),sstr_title.str().c_str(),200,0,5000,200,-1000,1500);
    sstr_name.str("");
    sstr_title.str("");

    sstr_name<<"Difference vs Sum for total charge of X and Y for X1Y1 for  "<<"MwpcC"<<i/2+2<<(char)((i%2)*6+76)<<" for only flag cut";
    sstr_title<<sstr_name.str()<<";TC_{X} + TC_{Y};TC_{X} - TC_{Y}";
    histDiffVsSumXYTotalChargeX1Y1OnlyFlagCut[i]=new TH2D(sstr_name.str().c_str(),sstr_title.str().c_str(),200,0,5000,200,-1000,1500);
    sstr_name.str("");
    sstr_title.str("");

    sstr_name<<"Comparison without and with all cut for Difference vs Sum for rejection of no cut and accptance of all cut for  "<<"MwpcC"<<i/2+2<<(char)((i%2)*6+76)<<" after only flag cut";
    sstr_title<<sstr_name.str()<<";TC_{X} + TC_{Y};TC_{X} - TC_{Y}";
    histDiffVsSumXYTotalChargeX1Y1RNoAAll[i]=new TH2D(sstr_name.str().c_str(),sstr_title.str().c_str(),200,0,5000,200,-1000,1500);
    sstr_name.str("");
    sstr_title.str("");

    sstr_name<<"Comparison without and with all cut for Difference vs Sum for acceptance of no cut and rejection of all cut for  "<<"MwpcC"<<i/2+2<<(char)((i%2)*6+76)<<" after only flag cut";
    sstr_title<<sstr_name.str()<<";TC_{X} + TC_{Y};TC_{X} - TC_{Y}";
    histDiffVsSumXYTotalChargeX1Y1ANoRAll[i]=new TH2D(sstr_name.str().c_str(),sstr_title.str().c_str(),200,0,5000,200,-1000,1500);
    sstr_name.str("");
    sstr_title.str("");

    sstr_name<<"Comparison without and with all cut for Difference vs Sum for acceptance of no cut and rejection of all cut without flag cut for  "<<"MwpcC"<<i/2+2<<(char)((i%2)*6+76)<<" after only flag cut";
    sstr_title<<sstr_name.str()<<";TC_{X} + TC_{Y};TC_{X} - TC_{Y}";
    histDiffVsSumXYTotalChargeX1Y1ANoRAllNoFlagCut[i]=new TH2D(sstr_name.str().c_str(),sstr_title.str().c_str(),200,0,5000,200,-1000,1500);
    sstr_name.str("");
    sstr_title.str("");

    sstr_name<<"Difference of total charge of 2 clusters in Y for X1Y2 for  "<<"MwpcC"<<i/2+2<<(char)((i%2)*6+76);
    sstr_title<<sstr_name.str()<<";TC_{Closer} - TC_{Further};Counts";
    histX1Y2DiffY2TotalCharge[i]=new TH1D(sstr_name.str().c_str(),sstr_title.str().c_str(),50,-1000,2000);
    sstr_name.str("");
    sstr_title.str("");


    sstr_name<<"Difference of total charge for further vs closer of X1Y2 for  "<<"MwpcC"<<i/2+2<<(char)((i%2)*6+76);
    sstr_title<<sstr_name.str()<<";TC_{X} - TC_{Y} for closer;TC_{X} - TC_{Y} for further";
    histDiffX1Y2TotalChargeFurtherVsCloser[i]=new TH2D(sstr_name.str().c_str(),sstr_title.str().c_str(),50,-1000,1000,50,-1000,1000);
    sstr_name.str("");
    sstr_title.str("");

    sstr_name<<"Difference vs Sum for total charge of X and Y for closer X1Y2 for  "<<"MwpcC"<<i/2+2<<(char)((i%2)*6+76);
    sstr_title<<sstr_name.str()<<";TC_{X} + TC_{Y};TC_{X} - TC_{Y}";
    histDiffVsSumXYTotalChargeX1Y2Closer[i]=new TH2D(sstr_name.str().c_str(),sstr_title.str().c_str(),200,0,5000,200,-1000,1500);
    sstr_name.str("");
    sstr_title.str("");

    sstr_name<<"Difference vs Sum for total charge of X and Y for further X1Y2 for  "<<"MwpcC"<<i/2+2<<(char)((i%2)*6+76);
    sstr_title<<sstr_name.str()<<";TC_{X} + TC_{Y};TC_{X} - TC_{Y}";
    histDiffVsSumXYTotalChargeX1Y2Further[i]=new TH2D(sstr_name.str().c_str(),sstr_title.str().c_str(),200,0,5000,200,-1000,1500);
    sstr_name.str("");
    sstr_title.str("");

    sstr_name<<"Difference vs Sum for total charge of X and Y for X1Y2 for  "<<"MwpcC"<<i/2+2<<(char)((i%2)*6+76);
    sstr_title<<sstr_name.str()<<";TC_{X} + TC_{Y};TC_{X} - TC_{Y}";
    histDiffVsSumXYTotalChargeX1Y2[i]=new TH2D(sstr_name.str().c_str(),sstr_title.str().c_str(),200,0,5000,200,-1000,1500);
    sstr_name.str("");
    sstr_title.str("");

    sstr_name<<"Difference of total charge of 2 clusters in X for X2Y1 for  "<<"MwpcC"<<i/2+2<<(char)((i%2)*6+76);
    sstr_title<<sstr_name.str()<<";TC_{Closer} - TC_{Further};Counts";
    histX2Y1DiffX2TotalCharge[i]=new TH1D(sstr_name.str().c_str(),sstr_title.str().c_str(),50,-1000,2000);
    sstr_name.str("");
    sstr_title.str("");


    sstr_name<<"Difference of total charge for further vs closer of X2Y1 for  "<<"MwpcC"<<i/2+2<<(char)((i%2)*6+76);
    sstr_title<<sstr_name.str()<<";TC_{X} - TC_{Y} for closer;TC_{X} - TC_{Y} for further";
    histDiffX2Y1TotalChargeFurtherVsCloser[i]=new TH2D(sstr_name.str().c_str(),sstr_title.str().c_str(),50,-1000,1000,50,-1000,1000);
    sstr_name.str("");
    sstr_title.str("");

    sstr_name<<"Difference vs Sum for total charge of X and Y for closer X2Y1 for  "<<"MwpcC"<<i/2+2<<(char)((i%2)*6+76);
    sstr_title<<sstr_name.str()<<";TC_{X} + TC_{Y};TC_{X} - TC_{Y}";
    histDiffVsSumXYTotalChargeX2Y1Closer[i]=new TH2D(sstr_name.str().c_str(),sstr_title.str().c_str(),200,0,5000,200,-1000,1500);
    sstr_name.str("");
    sstr_title.str("");

    sstr_name<<"Difference vs Sum for total charge of X and Y for further X2Y1 for  "<<"MwpcC"<<i/2+2<<(char)((i%2)*6+76);
    sstr_title<<sstr_name.str()<<";TC_{X} + TC_{Y};TC_{X} - TC_{Y}";
    histDiffVsSumXYTotalChargeX2Y1Further[i]=new TH2D(sstr_name.str().c_str(),sstr_title.str().c_str(),200,0,5000,200,-1000,1500);
    sstr_name.str("");
    sstr_title.str("");

    sstr_name<<"Difference vs Sum for total charge of X and Y for X2Y1 for  "<<"MwpcC"<<i/2+2<<(char)((i%2)*6+76);
    sstr_title<<sstr_name.str()<<";TC_{X} + TC_{Y};TC_{X} - TC_{Y}";
    histDiffVsSumXYTotalChargeX2Y1[i]=new TH2D(sstr_name.str().c_str(),sstr_title.str().c_str(),200,0,5000,200,-1000,1500);
    sstr_name.str("");
    sstr_title.str("");
  }
     
  return 0;
}

Long_t Det_MWPC::startup_comparexy(){
  getBranchObject("RawMwpcInfo",(TObject **) &treeRaw);
  gStyle->SetOptStat(0000);
  gStyle->SetOptFit(1111);
  return 0;
};


Long_t Det_MWPC::process_comparexy(){
  int* data[12]={treeRaw->ADC_C2X_L,treeRaw->ADC_C2X_R,treeRaw->ADC_C2Y_L,treeRaw->ADC_C2Y_R,treeRaw->ADC_C3X_L,treeRaw->ADC_C3X_R,treeRaw->ADC_C3Y_L,treeRaw->ADC_C3Y_R,treeRaw->ADC_C4X_L,treeRaw->ADC_C4X_R,treeRaw->ADC_C4Y_L,treeRaw->ADC_C4Y_R};
  const double* pedValue[12]={pedValueC2XL,pedValueC2XR,pedValueC2YL,pedValueC2YR,pedValueC3XL,pedValueC3XR,pedValueC3YL,pedValueC3YR,pedValueC4XL,pedValueC4XR,pedValueC4YL,pedValueC4YR};
  const double* pedSigma[12]={pedSigmaC2XL,pedSigmaC2XR,pedSigmaC2YL,pedSigmaC2YR,pedSigmaC3XL,pedSigmaC3XR,pedSigmaC3YL,pedSigmaC3YR,pedSigmaC4XL,pedSigmaC4XR,pedSigmaC4YL,pedSigmaC4YR}; 
  const double* totalChargeCut[12]={totChargeCutC2XL,totChargeCutC2XR,totChargeCutC2YL,totChargeCutC2YR,totChargeCutC3XL,totChargeCutC3XR,totChargeCutC3YL,totChargeCutC3YR,totChargeCutC4XL,totChargeCutC4XR,totChargeCutC4YL,totChargeCutC4YR};
  MwpcEvent *event[12];
  vector<vector<int> > clusStNum_allcut[12],clusStNum_nosingletcut[12],clusStNum_nocut[12];
  
  double *clusTotalCharge_allcut[12],*clusTotalCharge_nosingletcut[12],*clusTotalCharge_nocut[12];
  for(int i=0;i<12;i++){
    event[i]=new MwpcEvent(((i%4)/2*(-5)+(i/4)*(!((i%4)/2))+7)*8,!((i/4==0)&&((i%4/2))==0));
    clusStNum_allcut[i]=event[i]->clusterStripNumber(event[i],data[i],pedValue[i],pedSigma[i],1,totalChargeCut[i],1,singMiniSigmaClusterSecondApproach[i],1,doubMiniSigmaXClusterSecondApproach[i],doubMiniTotalChargeXCluster[i]);
    clusTotalCharge_allcut[i]=event[i]->getClusterTotCharge(event[i],data[i],pedValue[i],pedSigma[i],1,totalChargeCut[i],1,singMiniSigmaClusterSecondApproach[i],1,doubMiniSigmaXClusterSecondApproach[i],doubMiniTotalChargeXCluster[i]);
    clusStNum_nosingletcut[i]=event[i]->clusterStripNumber(event[i],data[i],pedValue[i],pedSigma[i],1,totalChargeCut[i],0,singMiniSigmaClusterSecondApproach[i],1,doubMiniSigmaXClusterSecondApproach[i],doubMiniTotalChargeXCluster[i]);
    clusTotalCharge_nosingletcut[i]=event[i]->getClusterTotCharge(event[i],data[i],pedValue[i],pedSigma[i],1,totalChargeCut[i],0,singMiniSigmaClusterSecondApproach[i],1,doubMiniSigmaXClusterSecondApproach[i],doubMiniTotalChargeXCluster[i]); 
    clusStNum_nocut[i]=event[i]->clusterStripNumber(event[i],data[i],pedValue[i],pedSigma[i],0,totalChargeCut[i],0,singMiniSigmaClusterSecondApproach[i],0,doubMiniSigmaXClusterSecondApproach[i],doubMiniTotalChargeXCluster[i]);
    clusTotalCharge_nocut[i]=event[i]->getClusterTotCharge(event[i],data[i],pedValue[i],pedSigma[i],0,totalChargeCut[i],0,singMiniSigmaClusterSecondApproach[i],0,doubMiniSigmaXClusterSecondApproach[i],doubMiniTotalChargeXCluster[i]);
    for(int j=0;j<(int)clusStNum_allcut[i].size();j++){
      histTotalCharge_allcut[i]->Fill(clusTotalCharge_allcut[i][j]);
    }
    for(int j=0;j<(int)clusStNum_nosingletcut[i].size();j++){
      histTotalCharge_nosingletcut[i]->Fill(clusTotalCharge_nosingletcut[i][j]);
    }
  }
  
  int flag;
  for(int i=0;i<2;i++){
    for(int j=0;j<3;j++){
       histXYNClusNoCut[2*j+i]->Fill(clusStNum_nocut[j*4+i].size(),clusStNum_nocut[j*4+i+2].size());
       histXYNClusBeFlagCut[2*j+i]->Fill(clusStNum_allcut[j*4+i].size(),clusStNum_allcut[j*4+i+2].size());
       histDiffVsSumXYTotalChargeX1Y1ANoRAllNoFlagCut[2*j+i]->Fill(clusTotalCharge_nocut[j*4+i][0]+clusTotalCharge_nocut[j*4+i+2][0],clusTotalCharge_nocut[j*4+i][0]-clusTotalCharge_nocut[j*4+i+2][0]);
    }

    flag=0;
    for(int j=0;j<6;j++){
      if(clusStNum_allcut[i+2*j].size()>0) flag++;
    }
    if(flag==6){
      for(int j=0;j<3;j++){
        histXYNClusAfFlagCut[2*j+i]->Fill(clusStNum_allcut[j*4+i].size(),clusStNum_allcut[j*4+i+2].size()); 
        if(clusStNum_allcut[j*4+i].size()==1 && clusStNum_allcut[j*4+i+2].size()==1){
          histDiffXYTotalChargeX1Y1[2*j+i]->Fill(clusTotalCharge_allcut[j*4+i][0]-clusTotalCharge_allcut[j*4+i+2][0]);
          histSumXYTotalChargeX1Y1[2*j+i]->Fill(clusTotalCharge_allcut[j*4+i][0]+clusTotalCharge_allcut[j*4+i+2][0]);
          histDiffVsSumXYTotalChargeX1Y1[2*j+i]->Fill(clusTotalCharge_allcut[j*4+i][0]+clusTotalCharge_allcut[j*4+i+2][0],clusTotalCharge_allcut[j*4+i][0]-clusTotalCharge_allcut[j*4+i+2][0]);
        }
        if(clusStNum_nocut[j*4+i].size()==1 && clusStNum_nocut[j*4+i+2].size()==1){
          histDiffVsSumXYTotalChargeX1Y1OnlyFlagCut[2*j+i]->Fill(clusTotalCharge_nocut[j*4+i][0]+clusTotalCharge_nocut[j*4+i+2][0],clusTotalCharge_nocut[j*4+i][0]-clusTotalCharge_nocut[j*4+i+2][0]);
        }
        if((clusStNum_allcut[j*4+i].size()==1 && clusStNum_allcut[j*4+i+2].size()==1) && (clusStNum_nocut[j*4+i].size()!=1 || clusStNum_nocut[j*4+i+2].size()!=1)){
          histDiffVsSumXYTotalChargeX1Y1RNoAAll[2*j+i]->Fill(clusTotalCharge_allcut[j*4+i][0]+clusTotalCharge_allcut[j*4+i+2][0],clusTotalCharge_allcut[j*4+i][0]-clusTotalCharge_allcut[j*4+i+2][0]);
        }
        if((clusStNum_nocut[j*4+i].size()==1 && clusStNum_nocut[j*4+i+2].size()==1) && ((clusStNum_allcut[j*4+i].size()==1 && clusStNum_allcut[j*4+i+2].size()==0) || (clusStNum_allcut[j*4+i].size()==0 && clusStNum_allcut[j*4+i+2].size()==1) || (clusStNum_allcut[j*4+i].size()==0 && clusStNum_allcut[j*4+i+2].size()==0))){
          histDiffVsSumXYTotalChargeX1Y1ANoRAll[2*j+i]->Fill(clusTotalCharge_nocut[j*4+i][0]+clusTotalCharge_nocut[j*4+i+2][0],clusTotalCharge_nocut[j*4+i][0]-clusTotalCharge_nocut[j*4+i+2][0]);
        }




        if(clusStNum_allcut[j*4+i].size()==1 && clusStNum_allcut[j*4+i+2].size()==2){
          histDiffVsSumXYTotalChargeX1Y2[2*j+i]->Fill(clusTotalCharge_allcut[j*4+i][0]+clusTotalCharge_allcut[j*4+i+2][0],clusTotalCharge_allcut[j*4+i][0]-clusTotalCharge_allcut[j*4+i+2][0]);
          histDiffVsSumXYTotalChargeX1Y2[2*j+i]->Fill(clusTotalCharge_allcut[j*4+i][0]+clusTotalCharge_allcut[j*4+i+2][1],clusTotalCharge_allcut[j*4+i][0]-clusTotalCharge_allcut[j*4+i+2][1]);

          if(fabs(clusTotalCharge_allcut[j*4+i][0]-clusTotalCharge_allcut[j*4+i+2][0]) < fabs(clusTotalCharge_allcut[j*4+i][0]-clusTotalCharge_allcut[j*4+i+2][1])){
            histDiffX1Y2TotalChargeFurtherVsCloser[2*j+i]->Fill(clusTotalCharge_allcut[j*4+i][0]-clusTotalCharge_allcut[j*4+i+2][0],clusTotalCharge_allcut[j*4+i][0]-clusTotalCharge_allcut[j*4+i+2][1]);
            histX1Y2DiffY2TotalCharge[2*j+i]->Fill(clusTotalCharge_allcut[j*4+i+2][0]-clusTotalCharge_allcut[j*4+i+2][1]);
            histDiffVsSumXYTotalChargeX1Y2Closer[2*j+i]->Fill(clusTotalCharge_allcut[j*4+i][0]+clusTotalCharge_allcut[j*4+i+2][0],clusTotalCharge_allcut[j*4+i][0]-clusTotalCharge_allcut[j*4+i+2][0]);
            histDiffVsSumXYTotalChargeX1Y2Further[2*j+i]->Fill(clusTotalCharge_allcut[j*4+i][0]+clusTotalCharge_allcut[j*4+i+2][1],clusTotalCharge_allcut[j*4+i][0]-clusTotalCharge_allcut[j*4+i+2][1]);
          }
          else{
            histDiffX1Y2TotalChargeFurtherVsCloser[2*j+i]->Fill(clusTotalCharge_allcut[j*4+i][0]-clusTotalCharge_allcut[j*4+i+2][1],clusTotalCharge_allcut[j*4+i][0]-clusTotalCharge_allcut[j*4+i+2][0]);
            histX1Y2DiffY2TotalCharge[2*j+i]->Fill(clusTotalCharge_allcut[j*4+i+2][1]-clusTotalCharge_allcut[j*4+i+2][0]);
            histDiffVsSumXYTotalChargeX1Y2Closer[2*j+i]->Fill(clusTotalCharge_allcut[j*4+i][0]+clusTotalCharge_allcut[j*4+i+2][1],clusTotalCharge_allcut[j*4+i][0]-clusTotalCharge_allcut[j*4+i+2][1]);
            histDiffVsSumXYTotalChargeX1Y2Further[2*j+i]->Fill(clusTotalCharge_allcut[j*4+i][0]+clusTotalCharge_allcut[j*4+i+2][0],clusTotalCharge_allcut[j*4+i][0]-clusTotalCharge_allcut[j*4+i+2][0]);
          }
        }




        if(clusStNum_allcut[j*4+i].size()==2 && clusStNum_allcut[j*4+i+2].size()==1){
          histDiffVsSumXYTotalChargeX2Y1[2*j+i]->Fill(clusTotalCharge_allcut[j*4+i][0]+clusTotalCharge_allcut[j*4+i+2][0],clusTotalCharge_allcut[j*4+i][0]-clusTotalCharge_allcut[j*4+i+2][0]);
          histDiffVsSumXYTotalChargeX2Y1[2*j+i]->Fill(clusTotalCharge_allcut[j*4+i][1]+clusTotalCharge_allcut[j*4+i+2][0],clusTotalCharge_allcut[j*4+i][1]-clusTotalCharge_allcut[j*4+i+2][0]);
          if(fabs(clusTotalCharge_allcut[j*4+i][0]-clusTotalCharge_allcut[j*4+i+2][0]) < fabs(clusTotalCharge_allcut[j*4+i][1]-clusTotalCharge_allcut[j*4+i+2][0])){
            histDiffX2Y1TotalChargeFurtherVsCloser[2*j+i]->Fill(clusTotalCharge_allcut[j*4+i][0]-clusTotalCharge_allcut[j*4+i+2][0],clusTotalCharge_allcut[j*4+i][1]-clusTotalCharge_allcut[j*4+i+2][0]);
            histX2Y1DiffX2TotalCharge[2*j+i]->Fill(clusTotalCharge_allcut[j*4+i][0]-clusTotalCharge_allcut[j*4+i][1]);
            histDiffVsSumXYTotalChargeX2Y1Closer[2*j+i]->Fill(clusTotalCharge_allcut[j*4+i][0]+clusTotalCharge_allcut[j*4+i+2][0],clusTotalCharge_allcut[j*4+i][0]-clusTotalCharge_allcut[j*4+i+2][0]);
            histDiffVsSumXYTotalChargeX2Y1Further[2*j+i]->Fill(clusTotalCharge_allcut[j*4+i][1]+clusTotalCharge_allcut[j*4+i+2][0],clusTotalCharge_allcut[j*4+i][1]-clusTotalCharge_allcut[j*4+i+2][0]);
          }
          else{
            histDiffX2Y1TotalChargeFurtherVsCloser[2*j+i]->Fill(clusTotalCharge_allcut[j*4+i][1]-clusTotalCharge_allcut[j*4+i+2][0],clusTotalCharge_allcut[j*4+i][0]-clusTotalCharge_allcut[j*4+i+2][0]);
            histX2Y1DiffX2TotalCharge[2*j+i]->Fill(clusTotalCharge_allcut[j*4+i][1]-clusTotalCharge_allcut[j*4+i][0]);
            histDiffVsSumXYTotalChargeX2Y1Closer[2*j+i]->Fill(clusTotalCharge_allcut[j*4+i][1]+clusTotalCharge_allcut[j*4+i+2][0],clusTotalCharge_allcut[j*4+i][1]-clusTotalCharge_allcut[j*4+i+2][0]);
            histDiffVsSumXYTotalChargeX2Y1Further[2*j+i]->Fill(clusTotalCharge_allcut[j*4+i][0]+clusTotalCharge_allcut[j*4+i+2][0],clusTotalCharge_allcut[j*4+i][0]-clusTotalCharge_allcut[j*4+i+2][0]);

          }
        }
      }
    }
  }
  return 0;
}

Long_t Det_MWPC::done_comparexy(){
  ostringstream sstr;

  TCanvas *canvTotalCharge=new TCanvas("Distribution of total charge","Distribution of total charge",2000,1200);
  canvTotalCharge->Divide(4,3);
  for(int i=0;i<12;i++){ 
    canvTotalCharge->cd(i+1);
    histTotalCharge_nosingletcut[i]->SetTitle("Distribution of total charge");
    histTotalCharge_nosingletcut[i]->SetMinimum(0);
    histTotalCharge_nosingletcut[i]->SetMaximum(6000);
    histTotalCharge_nosingletcut[i]->Draw();
    histTotalCharge_nosingletcut[i]->SetLineWidth(2);
    histTotalCharge_allcut[i]->SetLineColor(kRed);
    histTotalCharge_allcut[i]->Draw("same");
  }
  canvTotalCharge->cd(1);
  TLegend *legend=new TLegend(0.5,0.7,0.9,0.9);
  legend->AddEntry(histTotalCharge_allcut[0],"All cut applied","l");
  legend->AddEntry(histTotalCharge_nosingletcut[0],"No singlet threshold cut","l");
  legend->Draw();
  canvTotalCharge->SaveAs("Distribution of total charge.pdf");

  TCanvas *canvXYNClusNoCut=new TCanvas("Comparison of number of clusters between X and Y without any cuts","Comparison of number of clusters between X and Y without any cuts",600,360);
  canvXYNClusNoCut->Divide(2,3);
  for(int i=0;i<6;i++){
    canvXYNClusNoCut->cd(i+1);
    histXYNClusNoCut[i]->Draw("colz");
  }
  canvXYNClusNoCut->SaveAs("Comparison of number of clusters between X and Y without any cuts.pdf");

  TCanvas *canvXYNClusBeFlagCut=new TCanvas("Comparison of number of clusters between X and Y before flag cut","Comparison of number of clusters between X and Y before flag cut",600,360);
  canvXYNClusBeFlagCut->Divide(2,3);
  for(int i=0;i<6;i++){
    canvXYNClusBeFlagCut->cd(i+1);
    histXYNClusBeFlagCut[i]->Draw("colz");
  }
  canvXYNClusBeFlagCut->SaveAs("Comparison of number of clusters between X and Y before flag cut.pdf");

  TCanvas *canvXYNClusAfFlagCut=new TCanvas("Comparison of number of clusters between X and Y after flag cut","Comparison of number of clusters between X and Y after flag cut",600,360);
  canvXYNClusAfFlagCut->Divide(2,3);
  for(int i=0;i<6;i++){
    canvXYNClusAfFlagCut->cd(i+1);
    histXYNClusAfFlagCut[i]->Draw("colz");
  }
  canvXYNClusAfFlagCut->SaveAs("Comparison of number of clusters between X and Y after flag cut.pdf");

  double fitLow[6],fitUp[6];
  TF1 *fitGaus[6],*fitLandau[6];
  for(int i=0;i<6;i++){
    fitLow[i]=histDiffXYTotalChargeX1Y1[i]->GetBinLowEdge(histDiffXYTotalChargeX1Y1[i]->GetMaximumBin()-10);
    fitUp[i]=histDiffXYTotalChargeX1Y1[i]->GetBinLowEdge(histDiffXYTotalChargeX1Y1[i]->GetMaximumBin()+16);
    sstr<<"fitGaus"<<i;
    fitGaus[i]=new TF1(sstr.str().c_str(),"gaus",fitLow[i],fitUp[i]);
    sstr.str("");
    sstr<<"fitLandau"<<i;
    fitLandau[i]=new TF1(sstr.str().c_str(),"[0]*TMath::Landau(x,[1],[2],0)",fitLow[i],fitUp[i]);
    sstr.str("");
  }
  
  TCanvas *canvDiffXYTotalChargeX1Y1=new TCanvas("Difference of total charge between X and Y for X1Y1","Difference of total charge between X and Y for X1Y1",1000,900);
  canvDiffXYTotalChargeX1Y1->Divide(2,3);
  for(int i=0;i<6;i++){
    canvDiffXYTotalChargeX1Y1->cd(i+1);
    sstr<<"fitGaus"<<i;
    histDiffXYTotalChargeX1Y1[i]->Fit(sstr.str().c_str(),"R");
    sstr.str("");
    sstr<<"fitLandau"<<i;
    fitLandau[i]->SetParameter(1,fitGaus[i]->GetParameter(1));
    fitLandau[i]->SetParameter(2,fitGaus[i]->GetParameter(2));
    histDiffXYTotalChargeX1Y1[i]->Fit(sstr.str().c_str(),"R");
    sstr.str("");
    histDiffXYTotalChargeX1Y1[i]->Draw();
    fitGaus[i]->SetLineColor(kBlack);
    fitGaus[i]->Draw("same");
  }
  TLegend *diffLegend=new TLegend(0.1,0.7,0.3,0.9);
  diffLegend->AddEntry(fitGaus[0],"Gaussian","l");
  diffLegend->AddEntry(fitLandau[0],"Landau","l");
  canvDiffXYTotalChargeX1Y1->cd(1);
  diffLegend->Draw();
  canvDiffXYTotalChargeX1Y1->SaveAs("Difference of total charge between X and Y for X1Y1.pdf");

  double fitLowSum[6],fitUpSum[6];
  TF1 *fitGausSum[6],*fitLandauSum[6];
  for(int i=0;i<6;i++){
    fitLowSum[i]=histSumXYTotalChargeX1Y1[i]->GetBinLowEdge(histSumXYTotalChargeX1Y1[i]->GetMaximumBin()-30);
    fitUpSum[i]=histSumXYTotalChargeX1Y1[i]->GetBinLowEdge(histSumXYTotalChargeX1Y1[i]->GetMaximumBin()+51);
    sstr<<"fitGausSum"<<i;
    fitGausSum[i]=new TF1(sstr.str().c_str(),"gaus",fitLowSum[i],fitUpSum[i]);
    sstr.str("");
    sstr<<"fitLandauSum"<<i;
    fitLandauSum[i]=new TF1(sstr.str().c_str(),"[0]*TMath::Landau(x,[1],[2],0)",fitLowSum[i],fitUpSum[i]);
    sstr.str("");
  }

  TCanvas *canvSumXYTotalChargeX1Y1=new TCanvas("Sum of total charge of X and Y for X1Y1","Sum of total charge of X and Y for X1Y1",1000,900);
  canvSumXYTotalChargeX1Y1->Divide(2,3);
  for(int i=0;i<6;i++){
    canvSumXYTotalChargeX1Y1->cd(i+1);
    sstr<<"fitGausSum"<<i;
    histSumXYTotalChargeX1Y1[i]->Fit(sstr.str().c_str(),"R");
    sstr.str("");
    sstr<<"fitLandauSum"<<i;
    fitLandauSum[i]->SetParameter(1,fitGausSum[i]->GetParameter(1));
    fitLandauSum[i]->SetParameter(2,fitGausSum[i]->GetParameter(2));
    histSumXYTotalChargeX1Y1[i]->Fit(sstr.str().c_str(),"R");
    sstr.str("");
    histSumXYTotalChargeX1Y1[i]->Draw();
  }
  canvSumXYTotalChargeX1Y1->SaveAs("Sum of total charge of X and Y for X1Y1.pdf");

  TCanvas *canvDiffVsSumXYTotalChargeX1Y1=new TCanvas("Difference vs Sum for total charge of X and Y for X1Y1","Difference vs Sum for total charge of X and Y for X1Y1",1000,900);
  canvDiffVsSumXYTotalChargeX1Y1->Divide(2,3);
  for(int i=0;i<6;i++){
    canvDiffVsSumXYTotalChargeX1Y1->cd(i+1);
    gPad->SetLogz();
    histDiffVsSumXYTotalChargeX1Y1[i]->Draw("colz");
  }
  canvDiffVsSumXYTotalChargeX1Y1->SaveAs("Difference vs Sum for total charge of X and Y for X1Y1.pdf");

  TCanvas *canvDiffVsSumXYTotalChargeX1Y1OnlyFlagCut=new TCanvas("Difference vs Sum for total charge of X and Y for X1Y1 after only flag cut","Difference vs Sum for total charge of X and Y for X1Y1 after only flag cut",1000,900);
  canvDiffVsSumXYTotalChargeX1Y1OnlyFlagCut->Divide(2,3);
  for(int i=0;i<6;i++){
    canvDiffVsSumXYTotalChargeX1Y1OnlyFlagCut->cd(i+1);
    gPad->SetLogz();
    histDiffVsSumXYTotalChargeX1Y1OnlyFlagCut[i]->Draw("colz");
  }
  canvDiffVsSumXYTotalChargeX1Y1OnlyFlagCut->SaveAs("Difference vs Sum for total charge of X and Y for X1Y1 after only flag cut.pdf");

  TCanvas *canvDiffVsSumXYTotalChargeX1Y1RNoAAll=new TCanvas("Difference vs Sum for total charge of X and Y for rejection of no cut and acceptance of all cut for X1Y1 after only flag cut","Difference vs Sum for total charge of X and Y for rejection of no cut and acceptance of all cut for X1Y1 after only flag cut",1000,900);
  canvDiffVsSumXYTotalChargeX1Y1RNoAAll->Divide(2,3);
  for(int i=0;i<6;i++){
    canvDiffVsSumXYTotalChargeX1Y1RNoAAll->cd(i+1);
    gPad->SetLogz();
    histDiffVsSumXYTotalChargeX1Y1RNoAAll[i]->Draw("colz");
  }
  canvDiffVsSumXYTotalChargeX1Y1RNoAAll->SaveAs("Difference vs Sum for total charge of X and Y for rejection of no cut and acceptance of all cut for X1Y1 after only flag cut.pdf");

  TCanvas *canvDiffVsSumXYTotalChargeX1Y1ANoRAll=new TCanvas("Difference vs Sum for total charge of X and Y for acceptance of no cut and acceptance of all cut for X1Y1 after only flag cut","Difference vs Sum for total charge of X and Y for acceptance of no cut and acceptance of all cut for X1Y1 after only flag cut",1000,900);
  canvDiffVsSumXYTotalChargeX1Y1ANoRAll->Divide(2,3);
  for(int i=0;i<6;i++){
    canvDiffVsSumXYTotalChargeX1Y1ANoRAll->cd(i+1);
    gPad->SetLogz();
    histDiffVsSumXYTotalChargeX1Y1ANoRAll[i]->Draw("colz");
  }
  canvDiffVsSumXYTotalChargeX1Y1ANoRAll->SaveAs("Difference vs Sum for total charge of X and Y for acceptance of no cut and acceptance of all cut for X1Y1 after only flag cut.pdf");

  TCanvas *canvDiffVsSumXYTotalChargeX1Y1ANoRAllNoFlagCut=new TCanvas("Difference vs Sum for total charge of X and Y for acceptance of no cut and acceptance of all cut without flag cut for X1Y1 after only flag cut","Difference vs Sum for total charge of X and Y for acceptance of no cut and acceptance of all cut without flag cut for X1Y1 after only flag cut",1000,900);
  canvDiffVsSumXYTotalChargeX1Y1ANoRAllNoFlagCut->Divide(2,3);
  for(int i=0;i<6;i++){
    canvDiffVsSumXYTotalChargeX1Y1ANoRAllNoFlagCut->cd(i+1);
    gPad->SetLogz();
    histDiffVsSumXYTotalChargeX1Y1ANoRAllNoFlagCut[i]->Draw("colz");
  }
  canvDiffVsSumXYTotalChargeX1Y1ANoRAllNoFlagCut->SaveAs("Difference vs Sum for total charge of X and Y for acceptance of no cut and acceptance of all cut without flag cut for X1Y1 after only flag cut.pdf");

  TCanvas *canvX1Y2DiffY2TotalCharge=new TCanvas("Difference of total charge of 2 clusters in Y for X1Y2","Difference of total charge of 2 clusters in Y for X1Y2",1000,1500);
  canvX1Y2DiffY2TotalCharge->Divide(2,3);
  for(int i=0;i<6;i++){
    canvX1Y2DiffY2TotalCharge->cd(i+1);
    gPad->SetLogz();
    histX1Y2DiffY2TotalCharge[i]->Draw("colz");
  }
  canvX1Y2DiffY2TotalCharge->SaveAs("Difference of total charge of 2 clusters in Y for X1Y2.pdf");

  TCanvas *canvDiffX1Y2TotalChargeFurtherVsCloser=new TCanvas("Difference of total charge for further vs closer of X1Y2","Difference of total charge for further vs closer of X1Y2",1000,1500);
  canvDiffX1Y2TotalChargeFurtherVsCloser->Divide(2,3);
  for(int i=0;i<6;i++){
    canvDiffX1Y2TotalChargeFurtherVsCloser->cd(i+1);
    gPad->SetLogz();
    histDiffX1Y2TotalChargeFurtherVsCloser[i]->Draw("colz");
  }
  canvDiffX1Y2TotalChargeFurtherVsCloser->SaveAs("Difference of total charge for further vs closer of X1Y2.pdf");

  TCanvas *canvDiffVsSumXYTotalChargeX1Y2Closer=new TCanvas("Difference vs Sum for total charge of X and Y for closer X1Y2","Difference vs Sum for total charge of X and Y for closer X1Y2",1000,900);
  canvDiffVsSumXYTotalChargeX1Y2Closer->Divide(2,3);
  for(int i=0;i<6;i++){
    canvDiffVsSumXYTotalChargeX1Y2Closer->cd(i+1);
    gPad->SetLogz();
    histDiffVsSumXYTotalChargeX1Y2Closer[i]->Draw("colz");
    cutFuncUp[i]->Draw("same");
    cutFuncLow[i]->Draw("same");
  }
  canvDiffVsSumXYTotalChargeX1Y2Closer->SaveAs("Difference vs Sum for total charge of X and Y for closer of X1Y2.pdf");

  TCanvas *canvDiffVsSumXYTotalChargeX1Y2Further=new TCanvas("Difference vs Sum for total charge of X and Y for further X1Y2","Difference vs Sum for total charge of X and Y for further X1Y2",1000,900);
  canvDiffVsSumXYTotalChargeX1Y2Further->Divide(2,3);
  for(int i=0;i<6;i++){
    canvDiffVsSumXYTotalChargeX1Y2Further->cd(i+1);
    gPad->SetLogz();
    histDiffVsSumXYTotalChargeX1Y2Further[i]->Draw("colz");
    cutFuncUp[i]->Draw("same");
    cutFuncLow[i]->Draw("same");
  }
  canvDiffVsSumXYTotalChargeX1Y2Further->SaveAs("Difference vs Sum for total charge of X and Y for further of X1Y2.pdf");

  TCanvas *canvDiffVsSumXYTotalChargeX1Y2=new TCanvas("Difference vs Sum for total charge of X and Y for X1Y2","Difference vs Sum for total charge of X and Y for X1Y2",1000,900);
  canvDiffVsSumXYTotalChargeX1Y2->Divide(2,3);
  for(int i=0;i<6;i++){
    canvDiffVsSumXYTotalChargeX1Y2->cd(i+1);
    gPad->SetLogz();
    histDiffVsSumXYTotalChargeX1Y2[i]->Draw("colz");
    cutFuncUp[i]->Draw("same");
    cutFuncLow[i]->Draw("same");
  }
  canvDiffVsSumXYTotalChargeX1Y2->SaveAs("Difference vs Sum for total charge of X and Y for X1Y2.pdf");

  TCanvas *canvX2Y1DiffX2TotalCharge=new TCanvas("Difference of total charge of 2 clusters in X for X2Y1","Difference of total charge of 2 clusters in X for X2Y1",1000,1500);
  canvX2Y1DiffX2TotalCharge->Divide(2,3);
  for(int i=0;i<6;i++){
    canvX2Y1DiffX2TotalCharge->cd(i+1);
    gPad->SetLogz();
    histX2Y1DiffX2TotalCharge[i]->Draw("colz");
  }
  canvX2Y1DiffX2TotalCharge->SaveAs("Difference of total charge of 2 clusters in X for X2Y1.pdf");

  TCanvas *canvDiffX2Y1TotalChargeFurtherVsCloser=new TCanvas("Difference of total charge for further vs closer of X2Y1","Difference of total charge for further vs closer of X2Y1",1000,1500);
  canvDiffX2Y1TotalChargeFurtherVsCloser->Divide(2,3);
  for(int i=0;i<6;i++){
    canvDiffX2Y1TotalChargeFurtherVsCloser->cd(i+1);
    gPad->SetLogz();
    histDiffX2Y1TotalChargeFurtherVsCloser[i]->Draw("colz");
  }
  canvDiffX2Y1TotalChargeFurtherVsCloser->SaveAs("Difference of total charge for further vs closer of X2Y1.pdf");

  TCanvas *canvDiffVsSumXYTotalChargeX2Y1Closer=new TCanvas("Difference vs Sum for total charge of X and Y for closer X2Y1","Difference vs Sum for total charge of X and Y for closer X2Y1",1000,900);
  canvDiffVsSumXYTotalChargeX2Y1Closer->Divide(2,3);
  for(int i=0;i<6;i++){
    canvDiffVsSumXYTotalChargeX2Y1Closer->cd(i+1);
    gPad->SetLogz();
    histDiffVsSumXYTotalChargeX2Y1Closer[i]->Draw("colz");
    cutFuncUp[i]->Draw("same");
    cutFuncLow[i]->Draw("same");
  }
  canvDiffVsSumXYTotalChargeX2Y1Closer->SaveAs("Difference vs Sum for total charge of X and Y for closer of X2Y1.pdf");

  TCanvas *canvDiffVsSumXYTotalChargeX2Y1Further=new TCanvas("Difference vs Sum for total charge of X and Y for further X2Y1","Difference vs Sum for total charge of X and Y for further X2Y1",1000,900);
  canvDiffVsSumXYTotalChargeX2Y1Further->Divide(2,3);
  for(int i=0;i<6;i++){
    canvDiffVsSumXYTotalChargeX2Y1Further->cd(i+1);
    gPad->SetLogz();
    histDiffVsSumXYTotalChargeX2Y1Further[i]->Draw("colz");
    cutFuncUp[i]->Draw("same");
    cutFuncLow[i]->Draw("same");
  }
  canvDiffVsSumXYTotalChargeX2Y1Further->SaveAs("Difference vs Sum for total charge of X and Y for further of X2Y1.pdf");

  TCanvas *canvDiffVsSumXYTotalChargeX2Y1=new TCanvas("Difference vs Sum for total charge of X and Y for X2Y1","Difference vs Sum for total charge of X and Y for X2Y1",1000,900);
  canvDiffVsSumXYTotalChargeX2Y1->Divide(2,3);
  for(int i=0;i<6;i++){
    canvDiffVsSumXYTotalChargeX2Y1->cd(i+1);
    gPad->SetLogz();
    histDiffVsSumXYTotalChargeX2Y1[i]->Draw("colz");
    cutFuncUp[i]->Draw("same");
    cutFuncLow[i]->Draw("same");
  }
  canvDiffVsSumXYTotalChargeX2Y1->SaveAs("Difference vs Sum for total charge of X and Y for X2Y1.pdf");

  return 0;
};

