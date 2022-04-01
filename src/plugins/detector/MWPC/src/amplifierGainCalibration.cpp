#include <Det_MWPC.h>
#include <cmath>
#include <vector>
#include <fstream>
#include <iostream>
#include <TStyle.h>
#include <stdlib.h>
#include <sstream>
#include "TLine.h"
#include <string>
#include "TH1D.h"
#include "TCanvas.h"
#include "TText.h"
#include "TLegend.h"
#include "TMatrix.h"
#include "TLine.h"
#include "TLatex.h"
#include "TTree.h"
#include "TF1.h"
#include "C2X_Strip_Trans.h"
using namespace std;

static string inputFilePath(getenv("mwpc_inputFiles"));
static string inputFile=inputFilePath+"/gapNumTable/run3990-4059";
static ifstream input(inputFile.c_str());
static int runNum, eventNum, gapNumTOF1, gapNumTOF2;

static string rootFileDeadStrip=inputFilePath+"/MwpcAdjacentDeadStripCutFunction.root";
static TFile *rootDeadStrip=new TFile(rootFileDeadStrip.c_str());
static TF1 *cutFunc0313[6],*cutFunc2429[6];
//static  string start(getenv("start"));

static const int c3YRNumSt=16; // Strip number starts from 0
static int c3YRNumClusSing[6][c3YRNumSt]={{0}}, c3YRNumClusDoub[6][c3YRNumSt]={{0}},c3YRNumClusTrip[6][c3YRNumSt]={{0}}; 
static double c3YRSumChargeSing[6][c3YRNumSt]={{0}}, c3YRSumChargeDoubM1[6][c3YRNumSt]={{0}}, c3YRSumChargeDoub0[6][c3YRNumSt]={{0}}, c3YRSumChargeDoubP1[6][c3YRNumSt]={{0}}, c3YRSumChargeTripM2[6][c3YRNumSt]={{0}}, c3YRSumChargeTripM1[6][c3YRNumSt]={{0}}, c3YRSumChargeTrip0[6][c3YRNumSt]={{0}}, c3YRSumChargeTripP1[6][c3YRNumSt]={{0}},c3YRSumChargeTripP2[6][c3YRNumSt]={{0}};
 
static const int c3XRNumSt=64; // Strip number starts from 0
static int c3XRNumClusSing[6][c3XRNumSt]={{0}}, c3XRNumClusDoub[6][c3XRNumSt]={{0}},c3XRNumClusTrip[6][c3XRNumSt]={{0}};
static double c3XRSumChargeSing[6][c3XRNumSt]={{0}}, c3XRSumChargeDoubM1[6][c3XRNumSt]={{0}}, c3XRSumChargeDoub0[6][c3XRNumSt]={{0}}, c3XRSumChargeDoubP1[6][c3XRNumSt]={{0}}, c3XRSumChargeTripM2[6][c3XRNumSt]={{0}}, c3XRSumChargeTripM1[6][c3XRNumSt]={{0}}, c3XRSumChargeTrip0[6][c3XRNumSt]={{0}}, c3XRSumChargeTripP1[6][c3XRNumSt]={{0}},c3XRSumChargeTripP2[6][c3XRNumSt]={{0}};

static const int c3YLNumSt=16; // Strip number starts from 0
static int c3YLNumClusSing[6][c3YLNumSt]={{0}}, c3YLNumClusDoub[6][c3YLNumSt]={{0}},c3YLNumClusTrip[6][c3YLNumSt]={{0}};
static double c3YLSumChargeSing[6][c3YLNumSt]={{0}}, c3YLSumChargeDoubM1[6][c3YLNumSt]={{0}}, c3YLSumChargeDoub0[6][c3YLNumSt]={{0}}, c3YLSumChargeDoubP1[6][c3YLNumSt]={{0}}, c3YLSumChargeTripM2[6][c3YLNumSt]={{0}}, c3YLSumChargeTripM1[6][c3YLNumSt]={{0}}, c3YLSumChargeTrip0[6][c3YLNumSt]={{0}}, c3YLSumChargeTripP1[6][c3YLNumSt]={{0}},c3YLSumChargeTripP2[6][c3YLNumSt]={{0}};

static const int c3XLNumSt=64; // Strip number starts from 0
static int c3XLNumClusSing[6][c3XLNumSt]={{0}}, c3XLNumClusDoub[6][c3XLNumSt]={{0}},c3XLNumClusTrip[6][c3XLNumSt]={{0}};
static double c3XLSumChargeSing[6][c3XLNumSt]={{0}}, c3XLSumChargeDoubM1[6][c3XLNumSt]={{0}}, c3XLSumChargeDoub0[6][c3XLNumSt]={{0}}, c3XLSumChargeDoubP1[6][c3XLNumSt]={{0}}, c3XLSumChargeTripM2[6][c3XLNumSt]={{0}}, c3XLSumChargeTripM1[6][c3XLNumSt]={{0}}, c3XLSumChargeTrip0[6][c3XLNumSt]={{0}}, c3XLSumChargeTripP1[6][c3XLNumSt]={{0}},c3XLSumChargeTripP2[6][c3XLNumSt]={{0}};

static const int c4YRNumSt=16; // Strip number starts from 0
static int c4YRNumClusSing[6][c4YRNumSt]={{0}}, c4YRNumClusDoub[6][c4YRNumSt]={{0}},c4YRNumClusTrip[6][c4YRNumSt]={{0}};
static double c4YRSumChargeSing[6][c4YRNumSt]={{0}}, c4YRSumChargeDoubM1[6][c4YRNumSt]={{0}}, c4YRSumChargeDoub0[6][c4YRNumSt]={{0}}, c4YRSumChargeDoubP1[6][c4YRNumSt]={{0}}, c4YRSumChargeTripM2[6][c4YRNumSt]={{0}}, c4YRSumChargeTripM1[6][c4YRNumSt]={{0}}, c4YRSumChargeTrip0[6][c4YRNumSt]={{0}}, c4YRSumChargeTripP1[6][c4YRNumSt]={{0}},c4YRSumChargeTripP2[6][c4YRNumSt]={{0}};

static const int c4XRNumSt=72; // Strip number starts from 0
static int c4XRNumClusSing[6][c4XRNumSt]={{0}}, c4XRNumClusDoub[6][c4XRNumSt]={{0}},c4XRNumClusTrip[6][c4XRNumSt]={{0}};
static double c4XRSumChargeSing[6][c4XRNumSt]={{0}}, c4XRSumChargeDoubM1[6][c4XRNumSt]={{0}}, c4XRSumChargeDoub0[6][c4XRNumSt]={{0}}, c4XRSumChargeDoubP1[6][c4XRNumSt]={{0}}, c4XRSumChargeTripM2[6][c4XRNumSt]={{0}}, c4XRSumChargeTripM1[6][c4XRNumSt]={{0}}, c4XRSumChargeTrip0[6][c4XRNumSt]={{0}}, c4XRSumChargeTripP1[6][c4XRNumSt]={{0}},c4XRSumChargeTripP2[6][c4XRNumSt]={{0}};

static const int c4YLNumSt=16; // Strip number starts from 0
static int c4YLNumClusSing[6][c4YLNumSt]={{0}}, c4YLNumClusDoub[6][c4YLNumSt]={{0}},c4YLNumClusTrip[6][c4YLNumSt]={{0}};
static double c4YLSumChargeSing[6][c4YLNumSt]={{0}}, c4YLSumChargeDoubM1[6][c4YLNumSt]={{0}}, c4YLSumChargeDoub0[6][c4YLNumSt]={{0}}, c4YLSumChargeDoubP1[6][c4YLNumSt]={{0}}, c4YLSumChargeTripM2[6][c4YLNumSt]={{0}}, c4YLSumChargeTripM1[6][c4YLNumSt]={{0}}, c4YLSumChargeTrip0[6][c4YLNumSt]={{0}}, c4YLSumChargeTripP1[6][c4YLNumSt]={{0}},c4YLSumChargeTripP2[6][c4YLNumSt]={{0}};

static const int c4XLNumSt=72; // Strip number starts from 0
static int c4XLNumClusSing[6][c4XLNumSt]={{0}}, c4XLNumClusDoub[6][c4XLNumSt]={{0}},c4XLNumClusTrip[6][c4XLNumSt]={{0}};
static double c4XLSumChargeSing[6][c4XLNumSt]={{0}}, c4XLSumChargeDoubM1[6][c4XLNumSt]={{0}}, c4XLSumChargeDoub0[6][c4XLNumSt]={{0}}, c4XLSumChargeDoubP1[6][c4XLNumSt]={{0}}, c4XLSumChargeTripM2[6][c4XLNumSt]={{0}}, c4XLSumChargeTripM1[6][c4XLNumSt]={{0}}, c4XLSumChargeTrip0[6][c4XLNumSt]={{0}}, c4XLSumChargeTripP1[6][c4XLNumSt]={{0}},c4XLSumChargeTripP2[6][c4XLNumSt]={{0}};

static const int c2YRNumSt=16; // Strip number starts from 0
static int c2YRNumClusSing[6][c2YRNumSt]={{0}}, c2YRNumClusDoub[6][c2YRNumSt]={{0}},c2YRNumClusTrip[6][c2YRNumSt]={{0}};
static double c2YRSumChargeSing[6][c2YRNumSt]={{0}}, c2YRSumChargeDoubM1[6][c2YRNumSt]={{0}}, c2YRSumChargeDoub0[6][c2YRNumSt]={{0}}, c2YRSumChargeDoubP1[6][c2YRNumSt]={{0}}, c2YRSumChargeTripM2[6][c2YRNumSt]={{0}}, c2YRSumChargeTripM1[6][c2YRNumSt]={{0}}, c2YRSumChargeTrip0[6][c2YRNumSt]={{0}}, c2YRSumChargeTripP1[6][c2YRNumSt]={{0}},c2YRSumChargeTripP2[6][c2YRNumSt]={{0}};

static const int c2XRNumSt=56; // Strip number starts from 0
static int c2XRNumClusSing[6][c2XRNumSt]={{0}}, c2XRNumClusDoub[6][c2XRNumSt]={{0}},c2XRNumClusTrip[6][c2XRNumSt]={{0}};
static double c2XRSumChargeSing[6][c2XRNumSt]={{0}}, c2XRSumChargeDoubM1[6][c2XRNumSt]={{0}}, c2XRSumChargeDoub0[6][c2XRNumSt]={{0}}, c2XRSumChargeDoubP1[6][c2XRNumSt]={{0}}, c2XRSumChargeTripM2[6][c2XRNumSt]={{0}}, c2XRSumChargeTripM1[6][c2XRNumSt]={{0}}, c2XRSumChargeTrip0[6][c2XRNumSt]={{0}}, c2XRSumChargeTripP1[6][c2XRNumSt]={{0}},c2XRSumChargeTripP2[6][c2XRNumSt]={{0}};

static const int c2YLNumSt=16; // Strip number starts from 0
static int c2YLNumClusSing[6][c2YLNumSt]={{0}}, c2YLNumClusDoub[6][c2YLNumSt]={{0}},c2YLNumClusTrip[6][c2YLNumSt]={{0}};
static double c2YLSumChargeSing[6][c2YLNumSt]={{0}}, c2YLSumChargeDoubM1[6][c2YLNumSt]={{0}}, c2YLSumChargeDoub0[6][c2YLNumSt]={{0}}, c2YLSumChargeDoubP1[6][c2YLNumSt]={{0}}, c2YLSumChargeTripM2[6][c2YLNumSt]={{0}}, c2YLSumChargeTripM1[6][c2YLNumSt]={{0}}, c2YLSumChargeTrip0[6][c2YLNumSt]={{0}}, c2YLSumChargeTripP1[6][c2YLNumSt]={{0}},c2YLSumChargeTripP2[6][c2YLNumSt]={{0}};

static const int c2XLNumSt=56; // Strip number starts from 0
static int c2XLNumClusSing[6][c2XLNumSt]={{0}}, c2XLNumClusDoub[6][c2XLNumSt]={{0}},c2XLNumClusTrip[6][c2XLNumSt]={{0}};
static double c2XLSumChargeSing[6][c2XLNumSt]={{0}}, c2XLSumChargeDoubM1[6][c2XLNumSt]={{0}}, c2XLSumChargeDoub0[6][c2XLNumSt]={{0}}, c2XLSumChargeDoubP1[6][c2XLNumSt]={{0}}, c2XLSumChargeTripM2[6][c2XLNumSt]={{0}}, c2XLSumChargeTripM1[6][c2XLNumSt]={{0}}, c2XLSumChargeTrip0[6][c2XLNumSt]={{0}}, c2XLSumChargeTripP1[6][c2XLNumSt]={{0}},c2XLSumChargeTripP2[6][c2XLNumSt]={{0}};

Long_t Det_MWPC::histos_amplifierGainCalibration(){
  ostringstream sstr_name,sstr_title;

  return 0;
}

Long_t Det_MWPC::startup_amplifierGainCalibration(){
  getBranchObject("CaliMwpcInfo",(TObject **) &treeCali);
  gStyle->SetOptStat(0000);
  gStyle->SetOptFit(1111);

  ostringstream sstr;
  for(int i=0;i<6;i++){
    sstr<<"St0313_Gap"<<i+1<<"_low";
    cutFunc0313[i]=(TF1*)rootDeadStrip->Get(sstr.str().c_str());
    sstr.str("");
    sstr<<"St2429_Gap"<<i+1<<"_low";
    cutFunc2429[i]=(TF1*)rootDeadStrip->Get(sstr.str().c_str());
    sstr.str("");
  }

/*
  for(int i=0;i<atoi(start.c_str());i++){
    input>>runNum>>eventNum>>gapNumTOF1>>gapNumTOF2;
  }
*/
  return 0;
};


Long_t Det_MWPC::process_amplifierGainCalibration(){
  input>>runNum>>eventNum>>gapNumTOF1>>gapNumTOF2;
  if((UInt_t)runNum!=treeCali->run || (UInt_t)eventNum!=treeCali->event){
    cout<<"Error: Run number or event number are not consistent between two files!"<<endl;
    cout<<"For MWPC file: "<<"Run number - "<<treeCali->run<<"and  Event number - "<<treeCali->event<<endl;
    cout<<"For TRIUMF file: "<<"Run number - "<<runNum<<"and  Event number - "<<eventNum<<endl;
  }

  if(treeCali->thTCCutNumCluster[5]==1 && treeCali->thTCCutNumCluster[7]==1 && gapNumTOF2!=-1){
    for(int iGap=0;iGap<6;iGap++){
      if(gapNumTOF2==iGap+1){
        if(treeCali->thTCCutClusterStNum[7][0].size()==1){
          c3YRNumClusSing[iGap][treeCali->thTCCutClusterStNum[7][0][0]]++;
          c3YRSumChargeSing[iGap][treeCali->thTCCutClusterStNum[7][0][0]]+=treeCali->thTCCutClusterStCharge[7][0][0];
        }
        if(treeCali->thTCCutClusterStNum[7][0].size()==2){
          if(treeCali->thTCCutClusterStCharge[7][0][0]>treeCali->thTCCutClusterStCharge[7][0][1]){
            c3YRNumClusDoub[iGap][treeCali->thTCCutClusterStNum[7][0][0]]++;
            c3YRSumChargeDoub0[iGap][treeCali->thTCCutClusterStNum[7][0][0]]+=treeCali->thTCCutClusterStCharge[7][0][0];
            c3YRSumChargeDoubP1[iGap][treeCali->thTCCutClusterStNum[7][0][0]]+=treeCali->thTCCutClusterStCharge[7][0][1];
          }
          else{
            c3YRNumClusDoub[iGap][treeCali->thTCCutClusterStNum[7][0][1]]++;
            c3YRSumChargeDoub0[iGap][treeCali->thTCCutClusterStNum[7][0][1]]+=treeCali->thTCCutClusterStCharge[7][0][1];
            c3YRSumChargeDoubM1[iGap][treeCali->thTCCutClusterStNum[7][0][1]]+=treeCali->thTCCutClusterStCharge[7][0][0];
          }
        } 
        if(treeCali->thTCCutClusterStNum[7][0].size()==3){
          if(treeCali->thTCCutClusterStCharge[7][0][0]>treeCali->thTCCutClusterStCharge[7][0][1] && treeCali->thTCCutClusterStCharge[7][0][0]>treeCali->thTCCutClusterStCharge[7][0][2]){
            c3YRNumClusTrip[iGap][treeCali->thTCCutClusterStNum[7][0][0]]++;
            c3YRSumChargeTrip0[iGap][treeCali->thTCCutClusterStNum[7][0][0]]+=treeCali->thTCCutClusterStCharge[7][0][0];
            c3YRSumChargeTripP1[iGap][treeCali->thTCCutClusterStNum[7][0][0]]+=treeCali->thTCCutClusterStCharge[7][0][1];
            c3YRSumChargeTripP2[iGap][treeCali->thTCCutClusterStNum[7][0][0]]+=treeCali->thTCCutClusterStCharge[7][0][2];
          }
          if(treeCali->thTCCutClusterStCharge[7][0][1]>treeCali->thTCCutClusterStCharge[7][0][0] && treeCali->thTCCutClusterStCharge[7][0][1]>treeCali->thTCCutClusterStCharge[7][0][2]){
            c3YRNumClusTrip[iGap][treeCali->thTCCutClusterStNum[7][0][1]]++;
            c3YRSumChargeTrip0[iGap][treeCali->thTCCutClusterStNum[7][0][1]]+=treeCali->thTCCutClusterStCharge[7][0][1];
            c3YRSumChargeTripM1[iGap][treeCali->thTCCutClusterStNum[7][0][1]]+=treeCali->thTCCutClusterStCharge[7][0][0];
            c3YRSumChargeTripP1[iGap][treeCali->thTCCutClusterStNum[7][0][1]]+=treeCali->thTCCutClusterStCharge[7][0][2];
          }
          if(treeCali->thTCCutClusterStCharge[7][0][2]>treeCali->thTCCutClusterStCharge[7][0][0] && treeCali->thTCCutClusterStCharge[7][0][2]>treeCali->thTCCutClusterStCharge[7][0][1]){
            c3YRNumClusTrip[iGap][treeCali->thTCCutClusterStNum[7][0][2]]++;
            c3YRSumChargeTrip0[iGap][treeCali->thTCCutClusterStNum[7][0][2]]+=treeCali->thTCCutClusterStCharge[7][0][2];
            c3YRSumChargeTripM2[iGap][treeCali->thTCCutClusterStNum[7][0][2]]+=treeCali->thTCCutClusterStCharge[7][0][0];
            c3YRSumChargeTripM1[iGap][treeCali->thTCCutClusterStNum[7][0][2]]+=treeCali->thTCCutClusterStCharge[7][0][1];
          }
        }

        if(treeCali->thTCCutClusterStNum[5][0].size()==1){
          c3XRNumClusSing[iGap][treeCali->thTCCutClusterStNum[5][0][0]]++;
          c3XRSumChargeSing[iGap][treeCali->thTCCutClusterStNum[5][0][0]]+=treeCali->thTCCutClusterStCharge[5][0][0];
        }
        if(treeCali->thTCCutClusterStNum[5][0].size()==2){
          if(treeCali->thTCCutClusterStCharge[5][0][0]>treeCali->thTCCutClusterStCharge[5][0][1]){
            c3XRNumClusDoub[iGap][treeCali->thTCCutClusterStNum[5][0][0]]++;
            c3XRSumChargeDoub0[iGap][treeCali->thTCCutClusterStNum[5][0][0]]+=treeCali->thTCCutClusterStCharge[5][0][0];
            c3XRSumChargeDoubP1[iGap][treeCali->thTCCutClusterStNum[5][0][0]]+=treeCali->thTCCutClusterStCharge[5][0][1];
          }
          else{
            c3XRNumClusDoub[iGap][treeCali->thTCCutClusterStNum[5][0][1]]++;
            c3XRSumChargeDoub0[iGap][treeCali->thTCCutClusterStNum[5][0][1]]+=treeCali->thTCCutClusterStCharge[5][0][1];
            c3XRSumChargeDoubM1[iGap][treeCali->thTCCutClusterStNum[5][0][1]]+=treeCali->thTCCutClusterStCharge[5][0][0];
          }
        }
        if(treeCali->thTCCutClusterStNum[5][0].size()==3){
          if(treeCali->thTCCutClusterStCharge[5][0][0]>treeCali->thTCCutClusterStCharge[5][0][1] && treeCali->thTCCutClusterStCharge[5][0][0]>treeCali->thTCCutClusterStCharge[5][0][2]){
            c3XRNumClusTrip[iGap][treeCali->thTCCutClusterStNum[5][0][0]]++;
            c3XRSumChargeTrip0[iGap][treeCali->thTCCutClusterStNum[5][0][0]]+=treeCali->thTCCutClusterStCharge[5][0][0];
            c3XRSumChargeTripP1[iGap][treeCali->thTCCutClusterStNum[5][0][0]]+=treeCali->thTCCutClusterStCharge[5][0][1];
            c3XRSumChargeTripP2[iGap][treeCali->thTCCutClusterStNum[5][0][0]]+=treeCali->thTCCutClusterStCharge[5][0][2];
          }
          if(treeCali->thTCCutClusterStCharge[5][0][1]>treeCali->thTCCutClusterStCharge[5][0][0] && treeCali->thTCCutClusterStCharge[5][0][1]>treeCali->thTCCutClusterStCharge[5][0][2]){
            c3XRNumClusTrip[iGap][treeCali->thTCCutClusterStNum[5][0][1]]++;
            c3XRSumChargeTrip0[iGap][treeCali->thTCCutClusterStNum[5][0][1]]+=treeCali->thTCCutClusterStCharge[5][0][1];
            c3XRSumChargeTripM1[iGap][treeCali->thTCCutClusterStNum[5][0][1]]+=treeCali->thTCCutClusterStCharge[5][0][0];
            c3XRSumChargeTripP1[iGap][treeCali->thTCCutClusterStNum[5][0][1]]+=treeCali->thTCCutClusterStCharge[5][0][2];
          }
          if(treeCali->thTCCutClusterStCharge[5][0][2]>treeCali->thTCCutClusterStCharge[5][0][0] && treeCali->thTCCutClusterStCharge[5][0][2]>treeCali->thTCCutClusterStCharge[5][0][1]){
            c3XRNumClusTrip[iGap][treeCali->thTCCutClusterStNum[5][0][2]]++;
            c3XRSumChargeTrip0[iGap][treeCali->thTCCutClusterStNum[5][0][2]]+=treeCali->thTCCutClusterStCharge[5][0][2];
            c3XRSumChargeTripM2[iGap][treeCali->thTCCutClusterStNum[5][0][2]]+=treeCali->thTCCutClusterStCharge[5][0][0];
            c3XRSumChargeTripM1[iGap][treeCali->thTCCutClusterStNum[5][0][2]]+=treeCali->thTCCutClusterStCharge[5][0][1];
          }
        }
      }
    }
  }

  if(gapNumTOF2!=7 && treeCali->thTCCutNumCluster[4]==1 && treeCali->thTCCutNumCluster[6]==1 && gapNumTOF2!=-1){
    for(int iGap=0;iGap<6;iGap++){
      if(gapNumTOF2==iGap+7){
        if(treeCali->thTCCutClusterStNum[6][0].size()==1){
          c3YLNumClusSing[iGap][treeCali->thTCCutClusterStNum[6][0][0]]++;
          c3YLSumChargeSing[iGap][treeCali->thTCCutClusterStNum[6][0][0]]+=treeCali->thTCCutClusterStCharge[6][0][0];
       }
        if(treeCali->thTCCutClusterStNum[6][0].size()==2){
          if(treeCali->thTCCutClusterStCharge[6][0][0]>treeCali->thTCCutClusterStCharge[6][0][1]){
            c3YLNumClusDoub[iGap][treeCali->thTCCutClusterStNum[6][0][0]]++;
            c3YLSumChargeDoub0[iGap][treeCali->thTCCutClusterStNum[6][0][0]]+=treeCali->thTCCutClusterStCharge[6][0][0];
            c3YLSumChargeDoubP1[iGap][treeCali->thTCCutClusterStNum[6][0][0]]+=treeCali->thTCCutClusterStCharge[6][0][1];
          }
          else{
            c3YLNumClusDoub[iGap][treeCali->thTCCutClusterStNum[6][0][1]]++;
            c3YLSumChargeDoub0[iGap][treeCali->thTCCutClusterStNum[6][0][1]]+=treeCali->thTCCutClusterStCharge[6][0][1];
            c3YLSumChargeDoubM1[iGap][treeCali->thTCCutClusterStNum[6][0][1]]+=treeCali->thTCCutClusterStCharge[6][0][0];
          }
        }
        if(treeCali->thTCCutClusterStNum[6][0].size()==3){
          if(treeCali->thTCCutClusterStCharge[6][0][0]>treeCali->thTCCutClusterStCharge[6][0][1] && treeCali->thTCCutClusterStCharge[6][0][0]>treeCali->thTCCutClusterStCharge[6][0][2]){
            c3YLNumClusTrip[iGap][treeCali->thTCCutClusterStNum[6][0][0]]++;
            c3YLSumChargeTrip0[iGap][treeCali->thTCCutClusterStNum[6][0][0]]+=treeCali->thTCCutClusterStCharge[6][0][0];
            c3YLSumChargeTripP1[iGap][treeCali->thTCCutClusterStNum[6][0][0]]+=treeCali->thTCCutClusterStCharge[6][0][1];
            c3YLSumChargeTripP2[iGap][treeCali->thTCCutClusterStNum[6][0][0]]+=treeCali->thTCCutClusterStCharge[6][0][2];
          }
          if(treeCali->thTCCutClusterStCharge[6][0][1]>treeCali->thTCCutClusterStCharge[6][0][0] && treeCali->thTCCutClusterStCharge[6][0][1]>treeCali->thTCCutClusterStCharge[6][0][2]){
            c3YLNumClusTrip[iGap][treeCali->thTCCutClusterStNum[6][0][1]]++;
            c3YLSumChargeTrip0[iGap][treeCali->thTCCutClusterStNum[6][0][1]]+=treeCali->thTCCutClusterStCharge[6][0][1];
            c3YLSumChargeTripM1[iGap][treeCali->thTCCutClusterStNum[6][0][1]]+=treeCali->thTCCutClusterStCharge[6][0][0];
            c3YLSumChargeTripP1[iGap][treeCali->thTCCutClusterStNum[6][0][1]]+=treeCali->thTCCutClusterStCharge[6][0][2];
          }
          if(treeCali->thTCCutClusterStCharge[6][0][2]>treeCali->thTCCutClusterStCharge[6][0][0] && treeCali->thTCCutClusterStCharge[6][0][2]>treeCali->thTCCutClusterStCharge[6][0][1]){
            c3YLNumClusTrip[iGap][treeCali->thTCCutClusterStNum[6][0][2]]++;
            c3YLSumChargeTrip0[iGap][treeCali->thTCCutClusterStNum[6][0][2]]+=treeCali->thTCCutClusterStCharge[6][0][2];
            c3YLSumChargeTripM2[iGap][treeCali->thTCCutClusterStNum[6][0][2]]+=treeCali->thTCCutClusterStCharge[6][0][0];
            c3YLSumChargeTripM1[iGap][treeCali->thTCCutClusterStNum[6][0][2]]+=treeCali->thTCCutClusterStCharge[6][0][1];
          }
        }
        if(treeCali->thTCCutClusterStNum[4][0].size()==1){
          c3XLNumClusSing[iGap][treeCali->thTCCutClusterStNum[4][0][0]]++;
          c3XLSumChargeSing[iGap][treeCali->thTCCutClusterStNum[4][0][0]]+=treeCali->thTCCutClusterStCharge[4][0][0];
        }
        if(treeCali->thTCCutClusterStNum[4][0].size()==2){
          if(treeCali->thTCCutClusterStCharge[4][0][0]>treeCali->thTCCutClusterStCharge[4][0][1]){
            c3XLNumClusDoub[iGap][treeCali->thTCCutClusterStNum[4][0][0]]++;
            c3XLSumChargeDoub0[iGap][treeCali->thTCCutClusterStNum[4][0][0]]+=treeCali->thTCCutClusterStCharge[4][0][0];
            c3XLSumChargeDoubP1[iGap][treeCali->thTCCutClusterStNum[4][0][0]]+=treeCali->thTCCutClusterStCharge[4][0][1];
          }
          else{
            c3XLNumClusDoub[iGap][treeCali->thTCCutClusterStNum[4][0][1]]++;
            c3XLSumChargeDoub0[iGap][treeCali->thTCCutClusterStNum[4][0][1]]+=treeCali->thTCCutClusterStCharge[4][0][1];
            c3XLSumChargeDoubM1[iGap][treeCali->thTCCutClusterStNum[4][0][1]]+=treeCali->thTCCutClusterStCharge[4][0][0];
          }
        }
        if(treeCali->thTCCutClusterStNum[4][0].size()==3){
          if(treeCali->thTCCutClusterStCharge[4][0][0]>treeCali->thTCCutClusterStCharge[4][0][1] && treeCali->thTCCutClusterStCharge[4][0][0]>treeCali->thTCCutClusterStCharge[4][0][2]){
            c3XLNumClusTrip[iGap][treeCali->thTCCutClusterStNum[4][0][0]]++;
            c3XLSumChargeTrip0[iGap][treeCali->thTCCutClusterStNum[4][0][0]]+=treeCali->thTCCutClusterStCharge[4][0][0];
            c3XLSumChargeTripP1[iGap][treeCali->thTCCutClusterStNum[4][0][0]]+=treeCali->thTCCutClusterStCharge[4][0][1];
            c3XLSumChargeTripP2[iGap][treeCali->thTCCutClusterStNum[4][0][0]]+=treeCali->thTCCutClusterStCharge[4][0][2];
          }
          if(treeCali->thTCCutClusterStCharge[4][0][1]>treeCali->thTCCutClusterStCharge[4][0][0] && treeCali->thTCCutClusterStCharge[4][0][1]>treeCali->thTCCutClusterStCharge[4][0][2]){
            c3XLNumClusTrip[iGap][treeCali->thTCCutClusterStNum[4][0][1]]++;
            c3XLSumChargeTrip0[iGap][treeCali->thTCCutClusterStNum[4][0][1]]+=treeCali->thTCCutClusterStCharge[4][0][1];
            c3XLSumChargeTripM1[iGap][treeCali->thTCCutClusterStNum[4][0][1]]+=treeCali->thTCCutClusterStCharge[4][0][0];
            c3XLSumChargeTripP1[iGap][treeCali->thTCCutClusterStNum[4][0][1]]+=treeCali->thTCCutClusterStCharge[4][0][2];
          }
          if(treeCali->thTCCutClusterStCharge[4][0][2]>treeCali->thTCCutClusterStCharge[4][0][0] && treeCali->thTCCutClusterStCharge[4][0][2]>treeCali->thTCCutClusterStCharge[4][0][1]){
            c3XLNumClusTrip[iGap][treeCali->thTCCutClusterStNum[4][0][2]]++;
            c3XLSumChargeTrip0[iGap][treeCali->thTCCutClusterStNum[4][0][2]]+=treeCali->thTCCutClusterStCharge[4][0][2];
            c3XLSumChargeTripM2[iGap][treeCali->thTCCutClusterStNum[4][0][2]]+=treeCali->thTCCutClusterStCharge[4][0][0];
            c3XLSumChargeTripM1[iGap][treeCali->thTCCutClusterStNum[4][0][2]]+=treeCali->thTCCutClusterStCharge[4][0][1];
          }
        }
      }
    }
  }
  else if(gapNumTOF2==7 && treeCali->thTCCutNumClusterC3XLGap7Swap==1 && treeCali->thTCCutNumCluster[6]==1){
    if(treeCali->thTCCutClusterStNum[6][0].size()==1){
      c3YLNumClusSing[0][treeCali->thTCCutClusterStNum[6][0][0]]++;
      c3YLSumChargeSing[0][treeCali->thTCCutClusterStNum[6][0][0]]+=treeCali->thTCCutClusterStCharge[6][0][0];
    }
    if(treeCali->thTCCutClusterStNum[6][0].size()==2){
      if(treeCali->thTCCutClusterStCharge[6][0][0]>treeCali->thTCCutClusterStCharge[6][0][1]){
        c3YLNumClusDoub[0][treeCali->thTCCutClusterStNum[6][0][0]]++;
        c3YLSumChargeDoub0[0][treeCali->thTCCutClusterStNum[6][0][0]]+=treeCali->thTCCutClusterStCharge[6][0][0];
        c3YLSumChargeDoubP1[0][treeCali->thTCCutClusterStNum[6][0][0]]+=treeCali->thTCCutClusterStCharge[6][0][1];
      }
      else{
        c3YLNumClusDoub[0][treeCali->thTCCutClusterStNum[6][0][1]]++;
        c3YLSumChargeDoub0[0][treeCali->thTCCutClusterStNum[6][0][1]]+=treeCali->thTCCutClusterStCharge[6][0][1];
        c3YLSumChargeDoubM1[0][treeCali->thTCCutClusterStNum[6][0][1]]+=treeCali->thTCCutClusterStCharge[6][0][0];
      }
    }
    if(treeCali->thTCCutClusterStNum[6][0].size()==3){
      if(treeCali->thTCCutClusterStCharge[6][0][0]>treeCali->thTCCutClusterStCharge[6][0][1] && treeCali->thTCCutClusterStCharge[6][0][0]>treeCali->thTCCutClusterStCharge[6][0][2]){
        c3YLNumClusTrip[0][treeCali->thTCCutClusterStNum[6][0][0]]++;
        c3YLSumChargeTrip0[0][treeCali->thTCCutClusterStNum[6][0][0]]+=treeCali->thTCCutClusterStCharge[6][0][0];
        c3YLSumChargeTripP1[0][treeCali->thTCCutClusterStNum[6][0][0]]+=treeCali->thTCCutClusterStCharge[6][0][1];
        c3YLSumChargeTripP2[0][treeCali->thTCCutClusterStNum[6][0][0]]+=treeCali->thTCCutClusterStCharge[6][0][2];
      }
      if(treeCali->thTCCutClusterStCharge[6][0][1]>treeCali->thTCCutClusterStCharge[6][0][0] && treeCali->thTCCutClusterStCharge[6][0][1]>treeCali->thTCCutClusterStCharge[6][0][2]){
        c3YLNumClusTrip[0][treeCali->thTCCutClusterStNum[6][0][1]]++;
        c3YLSumChargeTrip0[0][treeCali->thTCCutClusterStNum[6][0][1]]+=treeCali->thTCCutClusterStCharge[6][0][1];
        c3YLSumChargeTripM1[0][treeCali->thTCCutClusterStNum[6][0][1]]+=treeCali->thTCCutClusterStCharge[6][0][0];
        c3YLSumChargeTripP1[0][treeCali->thTCCutClusterStNum[6][0][1]]+=treeCali->thTCCutClusterStCharge[6][0][2];
      }
      if(treeCali->thTCCutClusterStCharge[6][0][2]>treeCali->thTCCutClusterStCharge[6][0][0] && treeCali->thTCCutClusterStCharge[6][0][2]>treeCali->thTCCutClusterStCharge[6][0][1]){
        c3YLNumClusTrip[0][treeCali->thTCCutClusterStNum[6][0][2]]++;
        c3YLSumChargeTrip0[0][treeCali->thTCCutClusterStNum[6][0][2]]+=treeCali->thTCCutClusterStCharge[6][0][2];
        c3YLSumChargeTripM2[0][treeCali->thTCCutClusterStNum[6][0][2]]+=treeCali->thTCCutClusterStCharge[6][0][0];
        c3YLSumChargeTripM1[0][treeCali->thTCCutClusterStNum[6][0][2]]+=treeCali->thTCCutClusterStCharge[6][0][1];
      }
    }
    if(treeCali->thTCCutClusterStNumC3XLGap7Swap[0].size()==1){
      c3XLNumClusSing[0][treeCali->thTCCutClusterStNumC3XLGap7Swap[0][0]]++;
      c3XLSumChargeSing[0][treeCali->thTCCutClusterStNumC3XLGap7Swap[0][0]]+=treeCali->thTCCutClusterStChargeC3XLGap7Swap[0][0];
    }
    if(treeCali->thTCCutClusterStNumC3XLGap7Swap[0].size()==2){
      if(treeCali->thTCCutClusterStChargeC3XLGap7Swap[0][0]>treeCali->thTCCutClusterStChargeC3XLGap7Swap[0][1]){
      c3XLNumClusDoub[0][treeCali->thTCCutClusterStNumC3XLGap7Swap[0][0]]++;
      c3XLSumChargeDoub0[0][treeCali->thTCCutClusterStNumC3XLGap7Swap[0][0]]+=treeCali->thTCCutClusterStChargeC3XLGap7Swap[0][0];
      c3XLSumChargeDoubP1[0][treeCali->thTCCutClusterStNumC3XLGap7Swap[0][0]]+=treeCali->thTCCutClusterStChargeC3XLGap7Swap[0][1];
      }
      else{
        c3XLNumClusDoub[0][treeCali->thTCCutClusterStNumC3XLGap7Swap[0][1]]++;
        c3XLSumChargeDoub0[0][treeCali->thTCCutClusterStNumC3XLGap7Swap[0][1]]+=treeCali->thTCCutClusterStChargeC3XLGap7Swap[0][1];
        c3XLSumChargeDoubM1[0][treeCali->thTCCutClusterStNumC3XLGap7Swap[0][1]]+=treeCali->thTCCutClusterStChargeC3XLGap7Swap[0][0];
      }
    }
    if(treeCali->thTCCutClusterStNumC3XLGap7Swap[0].size()==3){
      if(treeCali->thTCCutClusterStChargeC3XLGap7Swap[0][0]>treeCali->thTCCutClusterStChargeC3XLGap7Swap[0][1] && treeCali->thTCCutClusterStChargeC3XLGap7Swap[0][0]>treeCali->thTCCutClusterStChargeC3XLGap7Swap[0][2]){
        c3XLNumClusTrip[0][treeCali->thTCCutClusterStNumC3XLGap7Swap[0][0]]++;
        c3XLSumChargeTrip0[0][treeCali->thTCCutClusterStNumC3XLGap7Swap[0][0]]+=treeCali->thTCCutClusterStChargeC3XLGap7Swap[0][0];
        c3XLSumChargeTripP1[0][treeCali->thTCCutClusterStNumC3XLGap7Swap[0][0]]+=treeCali->thTCCutClusterStChargeC3XLGap7Swap[0][1];
         c3XLSumChargeTripP2[0][treeCali->thTCCutClusterStNumC3XLGap7Swap[0][0]]+=treeCali->thTCCutClusterStChargeC3XLGap7Swap[0][2];
      }
      if(treeCali->thTCCutClusterStChargeC3XLGap7Swap[0][1]>treeCali->thTCCutClusterStChargeC3XLGap7Swap[0][0] && treeCali->thTCCutClusterStChargeC3XLGap7Swap[0][1]>treeCali->thTCCutClusterStChargeC3XLGap7Swap[0][2]){
        c3XLNumClusTrip[0][treeCali->thTCCutClusterStNumC3XLGap7Swap[0][1]]++;
        c3XLSumChargeTrip0[0][treeCali->thTCCutClusterStNumC3XLGap7Swap[0][1]]+=treeCali->thTCCutClusterStChargeC3XLGap7Swap[0][1];
        c3XLSumChargeTripM1[0][treeCali->thTCCutClusterStNumC3XLGap7Swap[0][1]]+=treeCali->thTCCutClusterStChargeC3XLGap7Swap[0][0];
        c3XLSumChargeTripP1[0][treeCali->thTCCutClusterStNumC3XLGap7Swap[0][1]]+=treeCali->thTCCutClusterStChargeC3XLGap7Swap[0][2];
      }
      if(treeCali->thTCCutClusterStChargeC3XLGap7Swap[0][2]>treeCali->thTCCutClusterStChargeC3XLGap7Swap[0][0] && treeCali->thTCCutClusterStChargeC3XLGap7Swap[0][2]>treeCali->thTCCutClusterStChargeC3XLGap7Swap[0][1]){
        c3XLNumClusTrip[0][treeCali->thTCCutClusterStNumC3XLGap7Swap[0][2]]++;
        c3XLSumChargeTrip0[0][treeCali->thTCCutClusterStNumC3XLGap7Swap[0][2]]+=treeCali->thTCCutClusterStChargeC3XLGap7Swap[0][2];
        c3XLSumChargeTripM2[0][treeCali->thTCCutClusterStNumC3XLGap7Swap[0][2]]+=treeCali->thTCCutClusterStChargeC3XLGap7Swap[0][0];
        c3XLSumChargeTripM1[0][treeCali->thTCCutClusterStNumC3XLGap7Swap[0][2]]+=treeCali->thTCCutClusterStChargeC3XLGap7Swap[0][1];
      }
    }
  }

  if(treeCali->thTCCutNumCluster[9]==1 && treeCali->thTCCutNumCluster[11]==1 && gapNumTOF2!=-1){
    for(int iGap=0;iGap<6;iGap++){
      if(gapNumTOF2==iGap+1){
        if(treeCali->thTCCutClusterStNum[11][0].size()==1){
          c4YRNumClusSing[iGap][treeCali->thTCCutClusterStNum[11][0][0]]++;
          c4YRSumChargeSing[iGap][treeCali->thTCCutClusterStNum[11][0][0]]+=treeCali->thTCCutClusterStCharge[11][0][0];
        }
        if(treeCali->thTCCutClusterStNum[11][0].size()==2){
          if(treeCali->thTCCutClusterStCharge[11][0][0]>treeCali->thTCCutClusterStCharge[11][0][1]){
            c4YRNumClusDoub[iGap][treeCali->thTCCutClusterStNum[11][0][0]]++;
            c4YRSumChargeDoub0[iGap][treeCali->thTCCutClusterStNum[11][0][0]]+=treeCali->thTCCutClusterStCharge[11][0][0];
            c4YRSumChargeDoubP1[iGap][treeCali->thTCCutClusterStNum[11][0][0]]+=treeCali->thTCCutClusterStCharge[11][0][1];
          }
          else{
            c4YRNumClusDoub[iGap][treeCali->thTCCutClusterStNum[11][0][1]]++;
            c4YRSumChargeDoub0[iGap][treeCali->thTCCutClusterStNum[11][0][1]]+=treeCali->thTCCutClusterStCharge[11][0][1];
            c4YRSumChargeDoubM1[iGap][treeCali->thTCCutClusterStNum[11][0][1]]+=treeCali->thTCCutClusterStCharge[11][0][0];
          }
        }
        if(treeCali->thTCCutClusterStNum[11][0].size()==3){
          if(treeCali->thTCCutClusterStCharge[11][0][0]>treeCali->thTCCutClusterStCharge[11][0][1] && treeCali->thTCCutClusterStCharge[11][0][0]>treeCali->thTCCutClusterStCharge[11][0][2]){
            c4YRNumClusTrip[iGap][treeCali->thTCCutClusterStNum[11][0][0]]++;
            c4YRSumChargeTrip0[iGap][treeCali->thTCCutClusterStNum[11][0][0]]+=treeCali->thTCCutClusterStCharge[11][0][0];
            c4YRSumChargeTripP1[iGap][treeCali->thTCCutClusterStNum[11][0][0]]+=treeCali->thTCCutClusterStCharge[11][0][1];
            c4YRSumChargeTripP2[iGap][treeCali->thTCCutClusterStNum[11][0][0]]+=treeCali->thTCCutClusterStCharge[11][0][2];
          }
          if(treeCali->thTCCutClusterStCharge[11][0][1]>treeCali->thTCCutClusterStCharge[11][0][0] && treeCali->thTCCutClusterStCharge[11][0][1]>treeCali->thTCCutClusterStCharge[11][0][2]){
            c4YRNumClusTrip[iGap][treeCali->thTCCutClusterStNum[11][0][1]]++;
            c4YRSumChargeTrip0[iGap][treeCali->thTCCutClusterStNum[11][0][1]]+=treeCali->thTCCutClusterStCharge[11][0][1];
            c4YRSumChargeTripM1[iGap][treeCali->thTCCutClusterStNum[11][0][1]]+=treeCali->thTCCutClusterStCharge[11][0][0];
            c4YRSumChargeTripP1[iGap][treeCali->thTCCutClusterStNum[11][0][1]]+=treeCali->thTCCutClusterStCharge[11][0][2];
          }
          if(treeCali->thTCCutClusterStCharge[11][0][2]>treeCali->thTCCutClusterStCharge[11][0][0] && treeCali->thTCCutClusterStCharge[11][0][2]>treeCali->thTCCutClusterStCharge[11][0][1]){
            c4YRNumClusTrip[iGap][treeCali->thTCCutClusterStNum[11][0][2]]++;
            c4YRSumChargeTrip0[iGap][treeCali->thTCCutClusterStNum[11][0][2]]+=treeCali->thTCCutClusterStCharge[11][0][2];
            c4YRSumChargeTripM2[iGap][treeCali->thTCCutClusterStNum[11][0][2]]+=treeCali->thTCCutClusterStCharge[11][0][0];
            c4YRSumChargeTripM1[iGap][treeCali->thTCCutClusterStNum[11][0][2]]+=treeCali->thTCCutClusterStCharge[11][0][1];
          }
        }

        if(treeCali->thTCCutClusterStNum[9][0].size()==1){
          if(treeCali->thTCCutClusterStNum[9][0][0]==15 && cutFunc0313[iGap]->Eval(treeCali->thTCCutClusterTC[9][0]+treeCali->thTCCutClusterTC[11][0])>treeCali->thTCCutClusterTC[9][0]-treeCali->thTCCutClusterTC[11][0]) continue;
          else if((treeCali->thTCCutClusterStNum[9][0][0]==17 || treeCali->thTCCutClusterStNum[9][0][0]==20 || treeCali->thTCCutClusterStNum[9][0][0]==22) && cutFunc2429[iGap]->Eval(treeCali->thTCCutClusterTC[9][0]+treeCali->thTCCutClusterTC[11][0])>treeCali->thTCCutClusterTC[9][0]-treeCali->thTCCutClusterTC[11][0]) continue;
          else{
            c4XRNumClusSing[iGap][treeCali->thTCCutClusterStNum[9][0][0]]++;
            c4XRSumChargeSing[iGap][treeCali->thTCCutClusterStNum[9][0][0]]+=treeCali->thTCCutClusterStCharge[9][0][0];
          }
        }
        if(treeCali->thTCCutClusterStNum[9][0].size()==2){
          if(treeCali->thTCCutClusterStCharge[9][0][0]>treeCali->thTCCutClusterStCharge[9][0][1]){
            if((treeCali->thTCCutClusterStNum[9][0][0]==17 || treeCali->thTCCutClusterStNum[9][0][0]==22) && cutFunc2429[iGap]->Eval(treeCali->thTCCutClusterTC[9][0]+treeCali->thTCCutClusterTC[11][0])>treeCali->thTCCutClusterTC[9][0]-treeCali->thTCCutClusterTC[11][0]) continue;
            else{
              c4XRNumClusDoub[iGap][treeCali->thTCCutClusterStNum[9][0][0]]++;
              c4XRSumChargeDoub0[iGap][treeCali->thTCCutClusterStNum[9][0][0]]+=treeCali->thTCCutClusterStCharge[9][0][0];
              c4XRSumChargeDoubP1[iGap][treeCali->thTCCutClusterStNum[9][0][0]]+=treeCali->thTCCutClusterStCharge[9][0][1];
            }
          }
          else{
            if(treeCali->thTCCutClusterStNum[9][0][1]==15 && cutFunc0313[iGap]->Eval(treeCali->thTCCutClusterTC[9][0]+treeCali->thTCCutClusterTC[11][0])>treeCali->thTCCutClusterTC[9][0]-treeCali->thTCCutClusterTC[11][0]) continue;
            else if(treeCali->thTCCutClusterStNum[9][0][1]==20 && cutFunc2429[iGap]->Eval(treeCali->thTCCutClusterTC[9][0]+treeCali->thTCCutClusterTC[11][0])>treeCali->thTCCutClusterTC[9][0]-treeCali->thTCCutClusterTC[11][0]) continue;
            else{
              c4XRNumClusDoub[iGap][treeCali->thTCCutClusterStNum[9][0][1]]++;
              c4XRSumChargeDoub0[iGap][treeCali->thTCCutClusterStNum[9][0][1]]+=treeCali->thTCCutClusterStCharge[9][0][1];
              c4XRSumChargeDoubM1[iGap][treeCali->thTCCutClusterStNum[9][0][1]]+=treeCali->thTCCutClusterStCharge[9][0][0];
            }
          }
        }
        if(treeCali->thTCCutClusterStNum[9][0].size()==3){
          if(treeCali->thTCCutClusterStCharge[9][0][0]>treeCali->thTCCutClusterStCharge[9][0][1] && treeCali->thTCCutClusterStCharge[9][0][0]>treeCali->thTCCutClusterStCharge[9][0][2]){
            if((treeCali->thTCCutClusterStNum[9][0][0]==17 || treeCali->thTCCutClusterStNum[9][0][0]==22) && cutFunc2429[iGap]->Eval(treeCali->thTCCutClusterTC[9][0]+treeCali->thTCCutClusterTC[11][0])>treeCali->thTCCutClusterTC[9][0]-treeCali->thTCCutClusterTC[11][0]) continue;
            else{
              c4XRNumClusTrip[iGap][treeCali->thTCCutClusterStNum[9][0][0]]++;
              c4XRSumChargeTrip0[iGap][treeCali->thTCCutClusterStNum[9][0][0]]+=treeCali->thTCCutClusterStCharge[9][0][0];
              c4XRSumChargeTripP1[iGap][treeCali->thTCCutClusterStNum[9][0][0]]+=treeCali->thTCCutClusterStCharge[9][0][1];
              c4XRSumChargeTripP2[iGap][treeCali->thTCCutClusterStNum[9][0][0]]+=treeCali->thTCCutClusterStCharge[9][0][2];
            }
          }
          if(treeCali->thTCCutClusterStCharge[9][0][1]>treeCali->thTCCutClusterStCharge[9][0][0] && treeCali->thTCCutClusterStCharge[9][0][1]>treeCali->thTCCutClusterStCharge[9][0][2]){
            c4XRNumClusTrip[iGap][treeCali->thTCCutClusterStNum[9][0][1]]++;
            c4XRSumChargeTrip0[iGap][treeCali->thTCCutClusterStNum[9][0][1]]+=treeCali->thTCCutClusterStCharge[9][0][1];
            c4XRSumChargeTripM1[iGap][treeCali->thTCCutClusterStNum[9][0][1]]+=treeCali->thTCCutClusterStCharge[9][0][0];
            c4XRSumChargeTripP1[iGap][treeCali->thTCCutClusterStNum[9][0][1]]+=treeCali->thTCCutClusterStCharge[9][0][2];
          }
          if(treeCali->thTCCutClusterStCharge[9][0][2]>treeCali->thTCCutClusterStCharge[9][0][0] && treeCali->thTCCutClusterStCharge[9][0][2]>treeCali->thTCCutClusterStCharge[9][0][1]){
            if(treeCali->thTCCutClusterStNum[9][0][2]==15 && cutFunc0313[iGap]->Eval(treeCali->thTCCutClusterTC[9][0]+treeCali->thTCCutClusterTC[11][0])>treeCali->thTCCutClusterTC[9][0]-treeCali->thTCCutClusterTC[11][0]) continue;
            else if(treeCali->thTCCutClusterStNum[9][0][2]==20 && cutFunc2429[iGap]->Eval(treeCali->thTCCutClusterTC[9][0]+treeCali->thTCCutClusterTC[11][0])>treeCali->thTCCutClusterTC[9][0]-treeCali->thTCCutClusterTC[11][0]) continue;
            else{
              c4XRNumClusTrip[iGap][treeCali->thTCCutClusterStNum[9][0][2]]++;
              c4XRSumChargeTrip0[iGap][treeCali->thTCCutClusterStNum[9][0][2]]+=treeCali->thTCCutClusterStCharge[9][0][2];
              c4XRSumChargeTripM2[iGap][treeCali->thTCCutClusterStNum[9][0][2]]+=treeCali->thTCCutClusterStCharge[9][0][0];
              c4XRSumChargeTripM1[iGap][treeCali->thTCCutClusterStNum[9][0][2]]+=treeCali->thTCCutClusterStCharge[9][0][1];
            }
          }
        }
      }
    }
  } 

  if(gapNumTOF2!=7 && gapNumTOF2!=11 && treeCali->thTCCutNumCluster[8]==1 && treeCali->thTCCutNumCluster[10]==1 && gapNumTOF2!=-1){
    for(int iGap=0;iGap<6;iGap++){
      if(gapNumTOF2==iGap+7){
        if(treeCali->thTCCutClusterStNum[10][0].size()==1){
          c4YLNumClusSing[iGap][treeCali->thTCCutClusterStNum[10][0][0]]++;
          c4YLSumChargeSing[iGap][treeCali->thTCCutClusterStNum[10][0][0]]+=treeCali->thTCCutClusterStCharge[10][0][0];
        }
        if(treeCali->thTCCutClusterStNum[10][0].size()==2){
          if(treeCali->thTCCutClusterStCharge[10][0][0]>treeCali->thTCCutClusterStCharge[10][0][1]){
            c4YLNumClusDoub[iGap][treeCali->thTCCutClusterStNum[10][0][0]]++;
            c4YLSumChargeDoub0[iGap][treeCali->thTCCutClusterStNum[10][0][0]]+=treeCali->thTCCutClusterStCharge[10][0][0];
            c4YLSumChargeDoubP1[iGap][treeCali->thTCCutClusterStNum[10][0][0]]+=treeCali->thTCCutClusterStCharge[10][0][1];
          }
          else{
            c4YLNumClusDoub[iGap][treeCali->thTCCutClusterStNum[10][0][1]]++;
            c4YLSumChargeDoub0[iGap][treeCali->thTCCutClusterStNum[10][0][1]]+=treeCali->thTCCutClusterStCharge[10][0][1];
            c4YLSumChargeDoubM1[iGap][treeCali->thTCCutClusterStNum[10][0][1]]+=treeCali->thTCCutClusterStCharge[10][0][0];
          }
        }
        if(treeCali->thTCCutClusterStNum[10][0].size()==3){
          if(treeCali->thTCCutClusterStCharge[10][0][0]>treeCali->thTCCutClusterStCharge[10][0][1] && treeCali->thTCCutClusterStCharge[10][0][0]>treeCali->thTCCutClusterStCharge[10][0][2]){
            c4YLNumClusTrip[iGap][treeCali->thTCCutClusterStNum[10][0][0]]++;
            c4YLSumChargeTrip0[iGap][treeCali->thTCCutClusterStNum[10][0][0]]+=treeCali->thTCCutClusterStCharge[10][0][0];
            c4YLSumChargeTripP1[iGap][treeCali->thTCCutClusterStNum[10][0][0]]+=treeCali->thTCCutClusterStCharge[10][0][1];
            c4YLSumChargeTripP2[iGap][treeCali->thTCCutClusterStNum[10][0][0]]+=treeCali->thTCCutClusterStCharge[10][0][2];
          }
          if(treeCali->thTCCutClusterStCharge[10][0][1]>treeCali->thTCCutClusterStCharge[10][0][0] && treeCali->thTCCutClusterStCharge[10][0][1]>treeCali->thTCCutClusterStCharge[10][0][2]){
            c4YLNumClusTrip[iGap][treeCali->thTCCutClusterStNum[10][0][1]]++;
            c4YLSumChargeTrip0[iGap][treeCali->thTCCutClusterStNum[10][0][1]]+=treeCali->thTCCutClusterStCharge[10][0][1];
            c4YLSumChargeTripM1[iGap][treeCali->thTCCutClusterStNum[10][0][1]]+=treeCali->thTCCutClusterStCharge[10][0][0];
            c4YLSumChargeTripP1[iGap][treeCali->thTCCutClusterStNum[10][0][1]]+=treeCali->thTCCutClusterStCharge[10][0][2];
          }
          if(treeCali->thTCCutClusterStCharge[10][0][2]>treeCali->thTCCutClusterStCharge[10][0][0] && treeCali->thTCCutClusterStCharge[10][0][2]>treeCali->thTCCutClusterStCharge[10][0][1]){
            c4YLNumClusTrip[iGap][treeCali->thTCCutClusterStNum[10][0][2]]++;
            c4YLSumChargeTrip0[iGap][treeCali->thTCCutClusterStNum[10][0][2]]+=treeCali->thTCCutClusterStCharge[10][0][2];
            c4YLSumChargeTripM2[iGap][treeCali->thTCCutClusterStNum[10][0][2]]+=treeCali->thTCCutClusterStCharge[10][0][0];
            c4YLSumChargeTripM1[iGap][treeCali->thTCCutClusterStNum[10][0][2]]+=treeCali->thTCCutClusterStCharge[10][0][1];
          }
        }
        if(treeCali->thTCCutClusterStNum[8][0].size()==1){
          c4XLNumClusSing[iGap][treeCali->thTCCutClusterStNum[8][0][0]]++;
          c4XLSumChargeSing[iGap][treeCali->thTCCutClusterStNum[8][0][0]]+=treeCali->thTCCutClusterStCharge[8][0][0];
        }
        if(treeCali->thTCCutClusterStNum[8][0].size()==2){
          if(treeCali->thTCCutClusterStCharge[8][0][0]>treeCali->thTCCutClusterStCharge[8][0][1]){
            c4XLNumClusDoub[iGap][treeCali->thTCCutClusterStNum[8][0][0]]++;
            c4XLSumChargeDoub0[iGap][treeCali->thTCCutClusterStNum[8][0][0]]+=treeCali->thTCCutClusterStCharge[8][0][0];
            c4XLSumChargeDoubP1[iGap][treeCali->thTCCutClusterStNum[8][0][0]]+=treeCali->thTCCutClusterStCharge[8][0][1];
          }
          else{
            c4XLNumClusDoub[iGap][treeCali->thTCCutClusterStNum[8][0][1]]++;
            c4XLSumChargeDoub0[iGap][treeCali->thTCCutClusterStNum[8][0][1]]+=treeCali->thTCCutClusterStCharge[8][0][1];
            c4XLSumChargeDoubM1[iGap][treeCali->thTCCutClusterStNum[8][0][1]]+=treeCali->thTCCutClusterStCharge[8][0][0];
          }
        }
        if(treeCali->thTCCutClusterStNum[8][0].size()==3){
          if(treeCali->thTCCutClusterStCharge[8][0][0]>treeCali->thTCCutClusterStCharge[8][0][1] && treeCali->thTCCutClusterStCharge[8][0][0]>treeCali->thTCCutClusterStCharge[8][0][2]){
            c4XLNumClusTrip[iGap][treeCali->thTCCutClusterStNum[8][0][0]]++;
            c4XLSumChargeTrip0[iGap][treeCali->thTCCutClusterStNum[8][0][0]]+=treeCali->thTCCutClusterStCharge[8][0][0];
            c4XLSumChargeTripP1[iGap][treeCali->thTCCutClusterStNum[8][0][0]]+=treeCali->thTCCutClusterStCharge[8][0][1];
            c4XLSumChargeTripP2[iGap][treeCali->thTCCutClusterStNum[8][0][0]]+=treeCali->thTCCutClusterStCharge[8][0][2];
          }
          if(treeCali->thTCCutClusterStCharge[8][0][1]>treeCali->thTCCutClusterStCharge[8][0][0] && treeCali->thTCCutClusterStCharge[8][0][1]>treeCali->thTCCutClusterStCharge[8][0][2]){
            c4XLNumClusTrip[iGap][treeCali->thTCCutClusterStNum[8][0][1]]++;
            c4XLSumChargeTrip0[iGap][treeCali->thTCCutClusterStNum[8][0][1]]+=treeCali->thTCCutClusterStCharge[8][0][1];
            c4XLSumChargeTripM1[iGap][treeCali->thTCCutClusterStNum[8][0][1]]+=treeCali->thTCCutClusterStCharge[8][0][0];
            c4XLSumChargeTripP1[iGap][treeCali->thTCCutClusterStNum[8][0][1]]+=treeCali->thTCCutClusterStCharge[8][0][2];
          }
          if(treeCali->thTCCutClusterStCharge[8][0][2]>treeCali->thTCCutClusterStCharge[8][0][0] && treeCali->thTCCutClusterStCharge[8][0][2]>treeCali->thTCCutClusterStCharge[8][0][1]){
            c4XLNumClusTrip[iGap][treeCali->thTCCutClusterStNum[8][0][2]]++;
            c4XLSumChargeTrip0[iGap][treeCali->thTCCutClusterStNum[8][0][2]]+=treeCali->thTCCutClusterStCharge[8][0][2];
            c4XLSumChargeTripM2[iGap][treeCali->thTCCutClusterStNum[8][0][2]]+=treeCali->thTCCutClusterStCharge[8][0][0];
            c4XLSumChargeTripM1[iGap][treeCali->thTCCutClusterStNum[8][0][2]]+=treeCali->thTCCutClusterStCharge[8][0][1];
          }
        }
      }
    }
  }
  else if(gapNumTOF2==7 && treeCali->thTCCutNumClusterC4XLGap7Swap==1 && treeCali->thTCCutNumCluster[10]==1){
    if(treeCali->thTCCutClusterStNum[10][0].size()==1){
      c4YLNumClusSing[0][treeCali->thTCCutClusterStNum[10][0][0]]++;
      c4YLSumChargeSing[0][treeCali->thTCCutClusterStNum[10][0][0]]+=treeCali->thTCCutClusterStCharge[10][0][0];
      }
    if(treeCali->thTCCutClusterStNum[10][0].size()==2){
      if(treeCali->thTCCutClusterStCharge[10][0][0]>treeCali->thTCCutClusterStCharge[10][0][1]){
        c4YLNumClusDoub[0][treeCali->thTCCutClusterStNum[10][0][0]]++;
        c4YLSumChargeDoub0[0][treeCali->thTCCutClusterStNum[10][0][0]]+=treeCali->thTCCutClusterStCharge[10][0][0];
        c4YLSumChargeDoubP1[0][treeCali->thTCCutClusterStNum[10][0][0]]+=treeCali->thTCCutClusterStCharge[10][0][1];
      }
      else{
        c4YLNumClusDoub[0][treeCali->thTCCutClusterStNum[10][0][1]]++;
        c4YLSumChargeDoub0[0][treeCali->thTCCutClusterStNum[10][0][1]]+=treeCali->thTCCutClusterStCharge[10][0][1];
        c4YLSumChargeDoubM1[0][treeCali->thTCCutClusterStNum[10][0][1]]+=treeCali->thTCCutClusterStCharge[10][0][0];
      }
    }
    if(treeCali->thTCCutClusterStNum[10][0].size()==3){
      if(treeCali->thTCCutClusterStCharge[10][0][0]>treeCali->thTCCutClusterStCharge[10][0][1] && treeCali->thTCCutClusterStCharge[10][0][0]>treeCali->thTCCutClusterStCharge[10][0][2]){
        c4YLNumClusTrip[0][treeCali->thTCCutClusterStNum[10][0][0]]++;
        c4YLSumChargeTrip0[0][treeCali->thTCCutClusterStNum[10][0][0]]+=treeCali->thTCCutClusterStCharge[10][0][0];
        c4YLSumChargeTripP1[0][treeCali->thTCCutClusterStNum[10][0][0]]+=treeCali->thTCCutClusterStCharge[10][0][1];
        c4YLSumChargeTripP2[0][treeCali->thTCCutClusterStNum[10][0][0]]+=treeCali->thTCCutClusterStCharge[10][0][2];
      }
      if(treeCali->thTCCutClusterStCharge[10][0][1]>treeCali->thTCCutClusterStCharge[10][0][0] && treeCali->thTCCutClusterStCharge[10][0][1]>treeCali->thTCCutClusterStCharge[10][0][2]){
        c4YLNumClusTrip[0][treeCali->thTCCutClusterStNum[10][0][1]]++;
        c4YLSumChargeTrip0[0][treeCali->thTCCutClusterStNum[10][0][1]]+=treeCali->thTCCutClusterStCharge[10][0][1];
        c4YLSumChargeTripM1[0][treeCali->thTCCutClusterStNum[10][0][1]]+=treeCali->thTCCutClusterStCharge[10][0][0];
        c4YLSumChargeTripP1[0][treeCali->thTCCutClusterStNum[10][0][1]]+=treeCali->thTCCutClusterStCharge[10][0][2];
      }
      if(treeCali->thTCCutClusterStCharge[10][0][2]>treeCali->thTCCutClusterStCharge[10][0][0] && treeCali->thTCCutClusterStCharge[10][0][2]>treeCali->thTCCutClusterStCharge[10][0][1]){
        c4YLNumClusTrip[0][treeCali->thTCCutClusterStNum[10][0][2]]++;
        c4YLSumChargeTrip0[0][treeCali->thTCCutClusterStNum[10][0][2]]+=treeCali->thTCCutClusterStCharge[10][0][2];
        c4YLSumChargeTripM2[0][treeCali->thTCCutClusterStNum[10][0][2]]+=treeCali->thTCCutClusterStCharge[10][0][0];
        c4YLSumChargeTripM1[0][treeCali->thTCCutClusterStNum[10][0][2]]+=treeCali->thTCCutClusterStCharge[10][0][1];
      }
    }
    if(treeCali->thTCCutClusterStNumC4XLGap7Swap[0].size()==1){
      c4XLNumClusSing[0][treeCali->thTCCutClusterStNumC4XLGap7Swap[0][0]]++;
      c4XLSumChargeSing[0][treeCali->thTCCutClusterStNumC4XLGap7Swap[0][0]]+=treeCali->thTCCutClusterStChargeC4XLGap7Swap[0][0];
    }
    if(treeCali->thTCCutClusterStNumC4XLGap7Swap[0].size()==2){
      if(treeCali->thTCCutClusterStChargeC4XLGap7Swap[0][0]>treeCali->thTCCutClusterStChargeC4XLGap7Swap[0][1]){
        c4XLNumClusDoub[0][treeCali->thTCCutClusterStNumC4XLGap7Swap[0][0]]++;
        c4XLSumChargeDoub0[0][treeCali->thTCCutClusterStNumC4XLGap7Swap[0][0]]+=treeCali->thTCCutClusterStChargeC4XLGap7Swap[0][0];
        c4XLSumChargeDoubP1[0][treeCali->thTCCutClusterStNumC4XLGap7Swap[0][0]]+=treeCali->thTCCutClusterStChargeC4XLGap7Swap[0][1];
      }
      else{
        c4XLNumClusDoub[0][treeCali->thTCCutClusterStNumC4XLGap7Swap[0][1]]++;
        c4XLSumChargeDoub0[0][treeCali->thTCCutClusterStNumC4XLGap7Swap[0][1]]+=treeCali->thTCCutClusterStChargeC4XLGap7Swap[0][1];
        c4XLSumChargeDoubM1[0][treeCali->thTCCutClusterStNumC4XLGap7Swap[0][1]]+=treeCali->thTCCutClusterStChargeC4XLGap7Swap[0][0];
      }
    }
    if(treeCali->thTCCutClusterStNumC4XLGap7Swap[0].size()==3){
      if(treeCali->thTCCutClusterStChargeC4XLGap7Swap[0][0]>treeCali->thTCCutClusterStChargeC4XLGap7Swap[0][1] && treeCali->thTCCutClusterStChargeC4XLGap7Swap[0][0]>treeCali->thTCCutClusterStChargeC4XLGap7Swap[0][2]){
        c4XLNumClusTrip[0][treeCali->thTCCutClusterStNumC4XLGap7Swap[0][0]]++;
        c4XLSumChargeTrip0[0][treeCali->thTCCutClusterStNumC4XLGap7Swap[0][0]]+=treeCali->thTCCutClusterStChargeC4XLGap7Swap[0][0];
        c4XLSumChargeTripP1[0][treeCali->thTCCutClusterStNumC4XLGap7Swap[0][0]]+=treeCali->thTCCutClusterStChargeC4XLGap7Swap[0][1];
        c4XLSumChargeTripP2[0][treeCali->thTCCutClusterStNumC4XLGap7Swap[0][0]]+=treeCali->thTCCutClusterStChargeC4XLGap7Swap[0][2];
      }
      if(treeCali->thTCCutClusterStChargeC4XLGap7Swap[0][1]>treeCali->thTCCutClusterStChargeC4XLGap7Swap[0][0] && treeCali->thTCCutClusterStChargeC4XLGap7Swap[0][1]>treeCali->thTCCutClusterStChargeC4XLGap7Swap[0][2]){
        c4XLNumClusTrip[0][treeCali->thTCCutClusterStNumC4XLGap7Swap[0][1]]++;
        c4XLSumChargeTrip0[0][treeCali->thTCCutClusterStNumC4XLGap7Swap[0][1]]+=treeCali->thTCCutClusterStChargeC4XLGap7Swap[0][1];
        c4XLSumChargeTripM1[0][treeCali->thTCCutClusterStNumC4XLGap7Swap[0][1]]+=treeCali->thTCCutClusterStChargeC4XLGap7Swap[0][0];
        c4XLSumChargeTripP1[0][treeCali->thTCCutClusterStNumC4XLGap7Swap[0][1]]+=treeCali->thTCCutClusterStChargeC4XLGap7Swap[0][2];
      }
      if(treeCali->thTCCutClusterStChargeC4XLGap7Swap[0][2]>treeCali->thTCCutClusterStChargeC4XLGap7Swap[0][0] && treeCali->thTCCutClusterStChargeC4XLGap7Swap[0][2]>treeCali->thTCCutClusterStChargeC4XLGap7Swap[0][1]){
        c4XLNumClusTrip[0][treeCali->thTCCutClusterStNumC4XLGap7Swap[0][2]]++;
        c4XLSumChargeTrip0[0][treeCali->thTCCutClusterStNumC4XLGap7Swap[0][2]]+=treeCali->thTCCutClusterStChargeC4XLGap7Swap[0][2];
        c4XLSumChargeTripM2[0][treeCali->thTCCutClusterStNumC4XLGap7Swap[0][2]]+=treeCali->thTCCutClusterStChargeC4XLGap7Swap[0][0];
        c4XLSumChargeTripM1[0][treeCali->thTCCutClusterStNumC4XLGap7Swap[0][2]]+=treeCali->thTCCutClusterStChargeC4XLGap7Swap[0][1];
      }
    }
  }
  else if(gapNumTOF2==11 && treeCali->thTCCutNumClusterC4XLGap11Swap==1 && treeCali->thTCCutNumCluster[10]==1){
    if(treeCali->thTCCutClusterStNum[10][0].size()==1){
      c4YLNumClusSing[4][treeCali->thTCCutClusterStNum[10][0][0]]++;
      c4YLSumChargeSing[4][treeCali->thTCCutClusterStNum[10][0][0]]+=treeCali->thTCCutClusterStCharge[10][0][0];
    }
    if(treeCali->thTCCutClusterStNum[10][0].size()==2){
      if(treeCali->thTCCutClusterStCharge[10][0][0]>treeCali->thTCCutClusterStCharge[10][0][1]){
        c4YLNumClusDoub[4][treeCali->thTCCutClusterStNum[10][0][0]]++;
        c4YLSumChargeDoub0[4][treeCali->thTCCutClusterStNum[10][0][0]]+=treeCali->thTCCutClusterStCharge[10][0][0];
        c4YLSumChargeDoubP1[4][treeCali->thTCCutClusterStNum[10][0][0]]+=treeCali->thTCCutClusterStCharge[10][0][1];
      }
      else{
        c4YLNumClusDoub[4][treeCali->thTCCutClusterStNum[10][0][1]]++;
        c4YLSumChargeDoub0[4][treeCali->thTCCutClusterStNum[10][0][1]]+=treeCali->thTCCutClusterStCharge[10][0][1];
        c4YLSumChargeDoubM1[4][treeCali->thTCCutClusterStNum[10][0][1]]+=treeCali->thTCCutClusterStCharge[10][0][0];
      }
    }
    if(treeCali->thTCCutClusterStNum[10][0].size()==3){
      if(treeCali->thTCCutClusterStCharge[10][0][0]>treeCali->thTCCutClusterStCharge[10][0][1] && treeCali->thTCCutClusterStCharge[10][0][0]>treeCali->thTCCutClusterStCharge[10][0][2]){
        c4YLNumClusTrip[4][treeCali->thTCCutClusterStNum[10][0][0]]++;
        c4YLSumChargeTrip0[4][treeCali->thTCCutClusterStNum[10][0][0]]+=treeCali->thTCCutClusterStCharge[10][0][0];
        c4YLSumChargeTripP1[4][treeCali->thTCCutClusterStNum[10][0][0]]+=treeCali->thTCCutClusterStCharge[10][0][1];
        c4YLSumChargeTripP2[4][treeCali->thTCCutClusterStNum[10][0][0]]+=treeCali->thTCCutClusterStCharge[10][0][2];
      }
      if(treeCali->thTCCutClusterStCharge[10][0][1]>treeCali->thTCCutClusterStCharge[10][0][0] && treeCali->thTCCutClusterStCharge[10][0][1]>treeCali->thTCCutClusterStCharge[10][0][2]){
        c4YLNumClusTrip[4][treeCali->thTCCutClusterStNum[10][0][1]]++;
        c4YLSumChargeTrip0[4][treeCali->thTCCutClusterStNum[10][0][1]]+=treeCali->thTCCutClusterStCharge[10][0][1];
        c4YLSumChargeTripM1[4][treeCali->thTCCutClusterStNum[10][0][1]]+=treeCali->thTCCutClusterStCharge[10][0][0];
        c4YLSumChargeTripP1[4][treeCali->thTCCutClusterStNum[10][0][1]]+=treeCali->thTCCutClusterStCharge[10][0][2];
      }
      if(treeCali->thTCCutClusterStCharge[10][0][2]>treeCali->thTCCutClusterStCharge[10][0][0] && treeCali->thTCCutClusterStCharge[10][0][2]>treeCali->thTCCutClusterStCharge[10][0][1]){
        c4YLNumClusTrip[4][treeCali->thTCCutClusterStNum[10][0][2]]++;
        c4YLSumChargeTrip0[4][treeCali->thTCCutClusterStNum[10][0][2]]+=treeCali->thTCCutClusterStCharge[10][0][2];
        c4YLSumChargeTripM2[4][treeCali->thTCCutClusterStNum[10][0][2]]+=treeCali->thTCCutClusterStCharge[10][0][0];
        c4YLSumChargeTripM1[4][treeCali->thTCCutClusterStNum[10][0][2]]+=treeCali->thTCCutClusterStCharge[10][0][1];
      }
    }
    if(treeCali->thTCCutClusterStNumC4XLGap11Swap[0].size()==1){
      c4XLNumClusSing[4][treeCali->thTCCutClusterStNumC4XLGap11Swap[0][0]]++;
      c4XLSumChargeSing[4][treeCali->thTCCutClusterStNumC4XLGap11Swap[0][0]]+=treeCali->thTCCutClusterStChargeC4XLGap11Swap[0][0];
    }
    if(treeCali->thTCCutClusterStNumC4XLGap11Swap[0].size()==2){
      if(treeCali->thTCCutClusterStChargeC4XLGap11Swap[0][0]>treeCali->thTCCutClusterStChargeC4XLGap11Swap[0][1]){
        c4XLNumClusDoub[4][treeCali->thTCCutClusterStNumC4XLGap11Swap[0][0]]++;
        c4XLSumChargeDoub0[4][treeCali->thTCCutClusterStNumC4XLGap11Swap[0][0]]+=treeCali->thTCCutClusterStChargeC4XLGap11Swap[0][0];
        c4XLSumChargeDoubP1[4][treeCali->thTCCutClusterStNumC4XLGap11Swap[0][0]]+=treeCali->thTCCutClusterStChargeC4XLGap11Swap[0][1];
      }
      else{
        c4XLNumClusDoub[4][treeCali->thTCCutClusterStNumC4XLGap11Swap[0][1]]++;
        c4XLSumChargeDoub0[4][treeCali->thTCCutClusterStNumC4XLGap11Swap[0][1]]+=treeCali->thTCCutClusterStChargeC4XLGap11Swap[0][1];
        c4XLSumChargeDoubM1[4][treeCali->thTCCutClusterStNumC4XLGap11Swap[0][1]]+=treeCali->thTCCutClusterStChargeC4XLGap11Swap[0][0];
      }
    }
    if(treeCali->thTCCutClusterStNumC4XLGap11Swap[0].size()==3){
      if(treeCali->thTCCutClusterStChargeC4XLGap11Swap[0][0]>treeCali->thTCCutClusterStChargeC4XLGap11Swap[0][1] && treeCali->thTCCutClusterStChargeC4XLGap11Swap[0][0]>treeCali->thTCCutClusterStChargeC4XLGap11Swap[0][2]){
        c4XLNumClusTrip[4][treeCali->thTCCutClusterStNumC4XLGap11Swap[0][0]]++;
        c4XLSumChargeTrip0[4][treeCali->thTCCutClusterStNumC4XLGap11Swap[0][0]]+=treeCali->thTCCutClusterStChargeC4XLGap11Swap[0][0];
        c4XLSumChargeTripP1[4][treeCali->thTCCutClusterStNumC4XLGap11Swap[0][0]]+=treeCali->thTCCutClusterStChargeC4XLGap11Swap[0][1];
        c4XLSumChargeTripP2[4][treeCali->thTCCutClusterStNumC4XLGap11Swap[0][0]]+=treeCali->thTCCutClusterStChargeC4XLGap11Swap[0][2];
      }
      if(treeCali->thTCCutClusterStChargeC4XLGap11Swap[0][1]>treeCali->thTCCutClusterStChargeC4XLGap11Swap[0][0] && treeCali->thTCCutClusterStChargeC4XLGap11Swap[0][1]>treeCali->thTCCutClusterStChargeC4XLGap11Swap[0][2]){
        c4XLNumClusTrip[4][treeCali->thTCCutClusterStNumC4XLGap11Swap[0][1]]++;
        c4XLSumChargeTrip0[4][treeCali->thTCCutClusterStNumC4XLGap11Swap[0][1]]+=treeCali->thTCCutClusterStChargeC4XLGap11Swap[0][1];
        c4XLSumChargeTripM1[4][treeCali->thTCCutClusterStNumC4XLGap11Swap[0][1]]+=treeCali->thTCCutClusterStChargeC4XLGap11Swap[0][0];
        c4XLSumChargeTripP1[4][treeCali->thTCCutClusterStNumC4XLGap11Swap[0][1]]+=treeCali->thTCCutClusterStChargeC4XLGap11Swap[0][2];
      }
      if(treeCali->thTCCutClusterStChargeC4XLGap11Swap[0][2]>treeCali->thTCCutClusterStChargeC4XLGap11Swap[0][0] && treeCali->thTCCutClusterStChargeC4XLGap11Swap[0][2]>treeCali->thTCCutClusterStChargeC4XLGap11Swap[0][1]){
        c4XLNumClusTrip[4][treeCali->thTCCutClusterStNumC4XLGap11Swap[0][2]]++;
        c4XLSumChargeTrip0[4][treeCali->thTCCutClusterStNumC4XLGap11Swap[0][2]]+=treeCali->thTCCutClusterStChargeC4XLGap11Swap[0][2];
        c4XLSumChargeTripM2[4][treeCali->thTCCutClusterStNumC4XLGap11Swap[0][2]]+=treeCali->thTCCutClusterStChargeC4XLGap11Swap[0][0];
        c4XLSumChargeTripM1[4][treeCali->thTCCutClusterStNumC4XLGap11Swap[0][2]]+=treeCali->thTCCutClusterStChargeC4XLGap11Swap[0][1];
      }
    }
  }


  if(treeCali->thTCCutNumCluster[1]==1 && treeCali->thTCCutNumCluster[3]==1 && gapNumTOF2!=-1){
    for(int iGap=0;iGap<6;iGap++){
      if(gapNumTOF2==iGap+1){
        if(treeCali->thTCCutClusterStNum[3][0].size()==1){
          c2YRNumClusSing[iGap][treeCali->thTCCutClusterStNum[3][0][0]]++;
          c2YRSumChargeSing[iGap][treeCali->thTCCutClusterStNum[3][0][0]]+=treeCali->thTCCutClusterStCharge[3][0][0];
        }
        if(treeCali->thTCCutClusterStNum[3][0].size()==2){
          if(treeCali->thTCCutClusterStCharge[3][0][0]>treeCali->thTCCutClusterStCharge[3][0][1]){
            c2YRNumClusDoub[iGap][treeCali->thTCCutClusterStNum[3][0][0]]++;
            c2YRSumChargeDoub0[iGap][treeCali->thTCCutClusterStNum[3][0][0]]+=treeCali->thTCCutClusterStCharge[3][0][0];
            c2YRSumChargeDoubP1[iGap][treeCali->thTCCutClusterStNum[3][0][0]]+=treeCali->thTCCutClusterStCharge[3][0][1];
          }
          else{
            c2YRNumClusDoub[iGap][treeCali->thTCCutClusterStNum[3][0][1]]++;
            c2YRSumChargeDoub0[iGap][treeCali->thTCCutClusterStNum[3][0][1]]+=treeCali->thTCCutClusterStCharge[3][0][1];
            c2YRSumChargeDoubM1[iGap][treeCali->thTCCutClusterStNum[3][0][1]]+=treeCali->thTCCutClusterStCharge[3][0][0];
          }
        }
        if(treeCali->thTCCutClusterStNum[3][0].size()==3){
          if(treeCali->thTCCutClusterStCharge[3][0][0]>treeCali->thTCCutClusterStCharge[3][0][1] && treeCali->thTCCutClusterStCharge[3][0][0]>treeCali->thTCCutClusterStCharge[3][0][2]){
            c2YRNumClusTrip[iGap][treeCali->thTCCutClusterStNum[3][0][0]]++;
            c2YRSumChargeTrip0[iGap][treeCali->thTCCutClusterStNum[3][0][0]]+=treeCali->thTCCutClusterStCharge[3][0][0];
            c2YRSumChargeTripP1[iGap][treeCali->thTCCutClusterStNum[3][0][0]]+=treeCali->thTCCutClusterStCharge[3][0][1];
            c2YRSumChargeTripP2[iGap][treeCali->thTCCutClusterStNum[3][0][0]]+=treeCali->thTCCutClusterStCharge[3][0][2];
          }
          if(treeCali->thTCCutClusterStCharge[3][0][1]>treeCali->thTCCutClusterStCharge[3][0][0] && treeCali->thTCCutClusterStCharge[3][0][1]>treeCali->thTCCutClusterStCharge[3][0][2]){
            c2YRNumClusTrip[iGap][treeCali->thTCCutClusterStNum[3][0][1]]++;
            c2YRSumChargeTrip0[iGap][treeCali->thTCCutClusterStNum[3][0][1]]+=treeCali->thTCCutClusterStCharge[3][0][1];
            c2YRSumChargeTripM1[iGap][treeCali->thTCCutClusterStNum[3][0][1]]+=treeCali->thTCCutClusterStCharge[3][0][0];
            c2YRSumChargeTripP1[iGap][treeCali->thTCCutClusterStNum[3][0][1]]+=treeCali->thTCCutClusterStCharge[3][0][2];
          }
          if(treeCali->thTCCutClusterStCharge[3][0][2]>treeCali->thTCCutClusterStCharge[3][0][0] && treeCali->thTCCutClusterStCharge[3][0][2]>treeCali->thTCCutClusterStCharge[3][0][1]){
            c2YRNumClusTrip[iGap][treeCali->thTCCutClusterStNum[3][0][2]]++;
            c2YRSumChargeTrip0[iGap][treeCali->thTCCutClusterStNum[3][0][2]]+=treeCali->thTCCutClusterStCharge[3][0][2];
            c2YRSumChargeTripM2[iGap][treeCali->thTCCutClusterStNum[3][0][2]]+=treeCali->thTCCutClusterStCharge[3][0][0];
            c2YRSumChargeTripM1[iGap][treeCali->thTCCutClusterStNum[3][0][2]]+=treeCali->thTCCutClusterStCharge[3][0][1];
          }
        }

        if(treeCali->thTCCutClusterStNum[1][0].size()==1){
          c2XRNumClusSing[iGap][treeCali->thTCCutClusterStNum[1][0][0]]++;
          c2XRSumChargeSing[iGap][treeCali->thTCCutClusterStNum[1][0][0]]+=treeCali->thTCCutClusterStCharge[1][0][0];
        }
        if(treeCali->thTCCutClusterStNum[1][0].size()==2){
          if(treeCali->thTCCutClusterStCharge[1][0][0]>treeCali->thTCCutClusterStCharge[1][0][1]){
            c2XRNumClusDoub[iGap][treeCali->thTCCutClusterStNum[1][0][0]]++;
            c2XRSumChargeDoub0[iGap][treeCali->thTCCutClusterStNum[1][0][0]]+=treeCali->thTCCutClusterStCharge[1][0][0];
            c2XRSumChargeDoubP1[iGap][treeCali->thTCCutClusterStNum[1][0][0]]+=treeCali->thTCCutClusterStCharge[1][0][1];
          }
          else{
            c2XRNumClusDoub[iGap][treeCali->thTCCutClusterStNum[1][0][1]]++;
            c2XRSumChargeDoub0[iGap][treeCali->thTCCutClusterStNum[1][0][1]]+=treeCali->thTCCutClusterStCharge[1][0][1];
            c2XRSumChargeDoubM1[iGap][treeCali->thTCCutClusterStNum[1][0][1]]+=treeCali->thTCCutClusterStCharge[1][0][0];
          }
        }
        if(treeCali->thTCCutClusterStNum[1][0].size()==3){
          if(treeCali->thTCCutClusterStCharge[1][0][0]>treeCali->thTCCutClusterStCharge[1][0][1] && treeCali->thTCCutClusterStCharge[1][0][0]>treeCali->thTCCutClusterStCharge[1][0][2]){
            c2XRNumClusTrip[iGap][treeCali->thTCCutClusterStNum[1][0][0]]++;
            c2XRSumChargeTrip0[iGap][treeCali->thTCCutClusterStNum[1][0][0]]+=treeCali->thTCCutClusterStCharge[1][0][0];
            c2XRSumChargeTripP1[iGap][treeCali->thTCCutClusterStNum[1][0][0]]+=treeCali->thTCCutClusterStCharge[1][0][1];
            c2XRSumChargeTripP2[iGap][treeCali->thTCCutClusterStNum[1][0][0]]+=treeCali->thTCCutClusterStCharge[1][0][2];
          }
          if(treeCali->thTCCutClusterStCharge[1][0][1]>treeCali->thTCCutClusterStCharge[1][0][0] && treeCali->thTCCutClusterStCharge[1][0][1]>treeCali->thTCCutClusterStCharge[1][0][2]){
            c2XRNumClusTrip[iGap][treeCali->thTCCutClusterStNum[1][0][1]]++;
            c2XRSumChargeTrip0[iGap][treeCali->thTCCutClusterStNum[1][0][1]]+=treeCali->thTCCutClusterStCharge[1][0][1];
            c2XRSumChargeTripM1[iGap][treeCali->thTCCutClusterStNum[1][0][1]]+=treeCali->thTCCutClusterStCharge[1][0][0];
            c2XRSumChargeTripP1[iGap][treeCali->thTCCutClusterStNum[1][0][1]]+=treeCali->thTCCutClusterStCharge[1][0][2];
          }
          if(treeCali->thTCCutClusterStCharge[1][0][2]>treeCali->thTCCutClusterStCharge[1][0][0] && treeCali->thTCCutClusterStCharge[1][0][2]>treeCali->thTCCutClusterStCharge[1][0][1]){
            c2XRNumClusTrip[iGap][treeCali->thTCCutClusterStNum[1][0][2]]++;
            c2XRSumChargeTrip0[iGap][treeCali->thTCCutClusterStNum[1][0][2]]+=treeCali->thTCCutClusterStCharge[1][0][2];
            c2XRSumChargeTripM2[iGap][treeCali->thTCCutClusterStNum[1][0][2]]+=treeCali->thTCCutClusterStCharge[1][0][0];
            c2XRSumChargeTripM1[iGap][treeCali->thTCCutClusterStNum[1][0][2]]+=treeCali->thTCCutClusterStCharge[1][0][1];
          }
        }
      }
    }
  }

  if(treeCali->thTCCutNumCluster[0]==1 && treeCali->thTCCutNumCluster[2]==1 && gapNumTOF2!=-1){
    for(int iGap=0;iGap<6;iGap++){
      if(gapNumTOF2==iGap+7){
        if(treeCali->thTCCutClusterStNum[2][0].size()==1){
          c2YLNumClusSing[iGap][treeCali->thTCCutClusterStNum[2][0][0]]++;
          c2YLSumChargeSing[iGap][treeCali->thTCCutClusterStNum[2][0][0]]+=treeCali->thTCCutClusterStCharge[2][0][0];
        }
        if(treeCali->thTCCutClusterStNum[2][0].size()==2){
          if(treeCali->thTCCutClusterStCharge[2][0][0]>treeCali->thTCCutClusterStCharge[2][0][1]){
            c2YLNumClusDoub[iGap][treeCali->thTCCutClusterStNum[2][0][0]]++;
            c2YLSumChargeDoub0[iGap][treeCali->thTCCutClusterStNum[2][0][0]]+=treeCali->thTCCutClusterStCharge[2][0][0];
            c2YLSumChargeDoubP1[iGap][treeCali->thTCCutClusterStNum[2][0][0]]+=treeCali->thTCCutClusterStCharge[2][0][1];
          }
          else{
            c2YLNumClusDoub[iGap][treeCali->thTCCutClusterStNum[2][0][1]]++;
            c2YLSumChargeDoub0[iGap][treeCali->thTCCutClusterStNum[2][0][1]]+=treeCali->thTCCutClusterStCharge[2][0][1];
            c2YLSumChargeDoubM1[iGap][treeCali->thTCCutClusterStNum[2][0][1]]+=treeCali->thTCCutClusterStCharge[2][0][0];
          }
        }
        if(treeCali->thTCCutClusterStNum[2][0].size()==3){
          if(treeCali->thTCCutClusterStCharge[2][0][0]>treeCali->thTCCutClusterStCharge[2][0][1] && treeCali->thTCCutClusterStCharge[2][0][0]>treeCali->thTCCutClusterStCharge[2][0][2]){
            c2YLNumClusTrip[iGap][treeCali->thTCCutClusterStNum[2][0][0]]++;
            c2YLSumChargeTrip0[iGap][treeCali->thTCCutClusterStNum[2][0][0]]+=treeCali->thTCCutClusterStCharge[2][0][0];
            c2YLSumChargeTripP1[iGap][treeCali->thTCCutClusterStNum[2][0][0]]+=treeCali->thTCCutClusterStCharge[2][0][1];
            c2YLSumChargeTripP2[iGap][treeCali->thTCCutClusterStNum[2][0][0]]+=treeCali->thTCCutClusterStCharge[2][0][2];
          }
          if(treeCali->thTCCutClusterStCharge[2][0][1]>treeCali->thTCCutClusterStCharge[2][0][0] && treeCali->thTCCutClusterStCharge[2][0][1]>treeCali->thTCCutClusterStCharge[2][0][2]){
            c2YLNumClusTrip[iGap][treeCali->thTCCutClusterStNum[2][0][1]]++;
            c2YLSumChargeTrip0[iGap][treeCali->thTCCutClusterStNum[2][0][1]]+=treeCali->thTCCutClusterStCharge[2][0][1];
            c2YLSumChargeTripM1[iGap][treeCali->thTCCutClusterStNum[2][0][1]]+=treeCali->thTCCutClusterStCharge[2][0][0];
            c2YLSumChargeTripP1[iGap][treeCali->thTCCutClusterStNum[2][0][1]]+=treeCali->thTCCutClusterStCharge[2][0][2];
          }
          if(treeCali->thTCCutClusterStCharge[2][0][2]>treeCali->thTCCutClusterStCharge[2][0][0] && treeCali->thTCCutClusterStCharge[2][0][2]>treeCali->thTCCutClusterStCharge[2][0][1]){
            c2YLNumClusTrip[iGap][treeCali->thTCCutClusterStNum[2][0][2]]++;
            c2YLSumChargeTrip0[iGap][treeCali->thTCCutClusterStNum[2][0][2]]+=treeCali->thTCCutClusterStCharge[2][0][2];
            c2YLSumChargeTripM2[iGap][treeCali->thTCCutClusterStNum[2][0][2]]+=treeCali->thTCCutClusterStCharge[2][0][0];
            c2YLSumChargeTripM1[iGap][treeCali->thTCCutClusterStNum[2][0][2]]+=treeCali->thTCCutClusterStCharge[2][0][1];
          }
        }
        if(treeCali->thTCCutClusterStNum[0][0].size()==1){
          c2XLNumClusSing[iGap][treeCali->thTCCutClusterStNum[0][0][0]]++;
          c2XLSumChargeSing[iGap][treeCali->thTCCutClusterStNum[0][0][0]]+=treeCali->thTCCutClusterStCharge[0][0][0];
        }
        if(treeCali->thTCCutClusterStNum[0][0].size()==2){
          if(treeCali->thTCCutClusterStCharge[0][0][0]>treeCali->thTCCutClusterStCharge[0][0][1]){
            c2XLNumClusDoub[iGap][treeCali->thTCCutClusterStNum[0][0][0]]++;
            c2XLSumChargeDoub0[iGap][treeCali->thTCCutClusterStNum[0][0][0]]+=treeCali->thTCCutClusterStCharge[0][0][0];
            c2XLSumChargeDoubP1[iGap][treeCali->thTCCutClusterStNum[0][0][0]]+=treeCali->thTCCutClusterStCharge[0][0][1];
          }
          else{
            c2XLNumClusDoub[iGap][treeCali->thTCCutClusterStNum[0][0][1]]++;
            c2XLSumChargeDoub0[iGap][treeCali->thTCCutClusterStNum[0][0][1]]+=treeCali->thTCCutClusterStCharge[0][0][1];
            c2XLSumChargeDoubM1[iGap][treeCali->thTCCutClusterStNum[0][0][1]]+=treeCali->thTCCutClusterStCharge[0][0][0];
          }
        }
        if(treeCali->thTCCutClusterStNum[0][0].size()==3){
          if(treeCali->thTCCutClusterStCharge[0][0][0]>treeCali->thTCCutClusterStCharge[0][0][1] && treeCali->thTCCutClusterStCharge[0][0][0]>treeCali->thTCCutClusterStCharge[0][0][2]){
            c2XLNumClusTrip[iGap][treeCali->thTCCutClusterStNum[0][0][0]]++;
            c2XLSumChargeTrip0[iGap][treeCali->thTCCutClusterStNum[0][0][0]]+=treeCali->thTCCutClusterStCharge[0][0][0];
            c2XLSumChargeTripP1[iGap][treeCali->thTCCutClusterStNum[0][0][0]]+=treeCali->thTCCutClusterStCharge[0][0][1];
            c2XLSumChargeTripP2[iGap][treeCali->thTCCutClusterStNum[0][0][0]]+=treeCali->thTCCutClusterStCharge[0][0][2];
          }
          if(treeCali->thTCCutClusterStCharge[0][0][1]>treeCali->thTCCutClusterStCharge[0][0][0] && treeCali->thTCCutClusterStCharge[0][0][1]>treeCali->thTCCutClusterStCharge[0][0][2]){
            c2XLNumClusTrip[iGap][treeCali->thTCCutClusterStNum[0][0][1]]++;
            c2XLSumChargeTrip0[iGap][treeCali->thTCCutClusterStNum[0][0][1]]+=treeCali->thTCCutClusterStCharge[0][0][1];
            c2XLSumChargeTripM1[iGap][treeCali->thTCCutClusterStNum[0][0][1]]+=treeCali->thTCCutClusterStCharge[0][0][0];
            c2XLSumChargeTripP1[iGap][treeCali->thTCCutClusterStNum[0][0][1]]+=treeCali->thTCCutClusterStCharge[0][0][2];
          }
          if(treeCali->thTCCutClusterStCharge[0][0][2]>treeCali->thTCCutClusterStCharge[0][0][0] && treeCali->thTCCutClusterStCharge[0][0][2]>treeCali->thTCCutClusterStCharge[0][0][1]){
            c2XLNumClusTrip[iGap][treeCali->thTCCutClusterStNum[0][0][2]]++;
            c2XLSumChargeTrip0[iGap][treeCali->thTCCutClusterStNum[0][0][2]]+=treeCali->thTCCutClusterStCharge[0][0][2];
            c2XLSumChargeTripM2[iGap][treeCali->thTCCutClusterStNum[0][0][2]]+=treeCali->thTCCutClusterStCharge[0][0][0];
            c2XLSumChargeTripM1[iGap][treeCali->thTCCutClusterStNum[0][0][2]]+=treeCali->thTCCutClusterStCharge[0][0][1];
          }
        }
      }
    }
  }

  return 0;
}

static int statThreshold=3500;

static vector<int> c3YRAccpStNum[6],c3YRRejeStNum[6];
static double c3YRChargeAverLevel[6];
static int c3YRNumClusTot[6][c3YRNumSt];
static double c3YRaM2[6][c3YRNumSt]={{0}},c3YRaM1[6][c3YRNumSt]={{0}},c3YRa0[6][c3YRNumSt]={{0}},c3YRaP1[6][c3YRNumSt]={{0}},c3YRaP2[6][c3YRNumSt]={{0}};
static double c3YRb[6][c3YRNumSt][c3YRNumSt]={{{0}}},c3YRc[6][c3YRNumSt]={{0}};
static double c3YRGainWeightAccp[6][c3YRNumSt]={{0}},c3YRGainWeight[6][c3YRNumSt]={{0}};

static vector<int> c3XRAccpStNum[6],c3XRRejeStNum[6];
static double c3XRChargeAverLevel[6];
static int c3XRNumClusTot[6][c3XRNumSt];
static double c3XRaM2[6][c3XRNumSt]={{0}},c3XRaM1[6][c3XRNumSt]={{0}},c3XRa0[6][c3XRNumSt]={{0}},c3XRaP1[6][c3XRNumSt]={{0}},c3XRaP2[6][c3XRNumSt]={{0}};
static double c3XRb[6][c3XRNumSt][c3XRNumSt]={{{0}}},c3XRc[6][c3XRNumSt]={{0}};
static double c3XRGainWeightAccp[6][c3XRNumSt]={{0}},c3XRGainWeight[6][c3XRNumSt]={{0}};

static vector<int> c3YLAccpStNum[6],c3YLRejeStNum[6];
static double c3YLChargeAverLevel[6];
static int c3YLNumClusTot[6][c3YLNumSt];
static double c3YLaM2[6][c3YLNumSt]={{0}},c3YLaM1[6][c3YLNumSt]={{0}},c3YLa0[6][c3YLNumSt]={{0}},c3YLaP1[6][c3YLNumSt]={{0}},c3YLaP2[6][c3YLNumSt]={{0}};
static double c3YLb[6][c3YLNumSt][c3YLNumSt]={{{0}}},c3YLc[6][c3YLNumSt]={{0}};
static double c3YLGainWeightAccp[6][c3YLNumSt]={{0}},c3YLGainWeight[6][c3YLNumSt]={{0}};

static vector<int> c3XLAccpStNum[6],c3XLRejeStNum[6];
static double c3XLChargeAverLevel[6];
static int c3XLNumClusTot[6][c3XLNumSt];
static double c3XLaM2[6][c3XLNumSt]={{0}},c3XLaM1[6][c3XLNumSt]={{0}},c3XLa0[6][c3XLNumSt]={{0}},c3XLaP1[6][c3XLNumSt]={{0}},c3XLaP2[6][c3XLNumSt]={{0}};
static double c3XLb[6][c3XLNumSt][c3XLNumSt]={{{0}}},c3XLc[6][c3XLNumSt]={{0}};
static double c3XLGainWeightAccp[6][c3XLNumSt]={{0}},c3XLGainWeight[6][c3XLNumSt]={{0}};

static vector<int> c4YRAccpStNum[6],c4YRRejeStNum[6];
static double c4YRChargeAverLevel[6];
static int c4YRNumClusTot[6][c4YRNumSt];
static double c4YRaM2[6][c4YRNumSt]={{0}},c4YRaM1[6][c4YRNumSt]={{0}},c4YRa0[6][c4YRNumSt]={{0}},c4YRaP1[6][c4YRNumSt]={{0}},c4YRaP2[6][c4YRNumSt]={{0}};
static double c4YRb[6][c4YRNumSt][c4YRNumSt]={{{0}}},c4YRc[6][c4YRNumSt]={{0}};
static double c4YRGainWeightAccp[6][c4YRNumSt]={{0}},c4YRGainWeight[6][c4YRNumSt]={{0}};

static vector<int> c4XRAccpStNum[6],c4XRRejeStNum[6];
static double c4XRChargeAverLevel[6];
static int c4XRNumClusTot[6][c4XRNumSt];
static double c4XRaM2[6][c4XRNumSt]={{0}},c4XRaM1[6][c4XRNumSt]={{0}},c4XRa0[6][c4XRNumSt]={{0}},c4XRaP1[6][c4XRNumSt]={{0}},c4XRaP2[6][c4XRNumSt]={{0}};
static double c4XRb[6][c4XRNumSt][c4XRNumSt]={{{0}}},c4XRc[6][c4XRNumSt]={{0}};
static double c4XRGainWeightAccp[6][c4XRNumSt]={{0}},c4XRGainWeight[6][c4XRNumSt]={{0}};

static vector<int> c4YLAccpStNum[6],c4YLRejeStNum[6];
static double c4YLChargeAverLevel[6];
static int c4YLNumClusTot[6][c4YLNumSt];
static double c4YLaM2[6][c4YLNumSt]={{0}},c4YLaM1[6][c4YLNumSt]={{0}},c4YLa0[6][c4YLNumSt]={{0}},c4YLaP1[6][c4YLNumSt]={{0}},c4YLaP2[6][c4YLNumSt]={{0}};
static double c4YLb[6][c4YLNumSt][c4YLNumSt]={{{0}}},c4YLc[6][c4YLNumSt]={{0}};
static double c4YLGainWeightAccp[6][c4YLNumSt]={{0}},c4YLGainWeight[6][c4YLNumSt]={{0}};

static vector<int> c4XLAccpStNum[6],c4XLRejeStNum[6];
static double c4XLChargeAverLevel[6];
static int c4XLNumClusTot[6][c4XLNumSt];
static double c4XLaM2[6][c4XLNumSt]={{0}},c4XLaM1[6][c4XLNumSt]={{0}},c4XLa0[6][c4XLNumSt]={{0}},c4XLaP1[6][c4XLNumSt]={{0}},c4XLaP2[6][c4XLNumSt]={{0}};
static double c4XLb[6][c4XLNumSt][c4XLNumSt]={{{0}}},c4XLc[6][c4XLNumSt]={{0}};
static double c4XLGainWeightAccp[6][c4XLNumSt]={{0}},c4XLGainWeight[6][c4XLNumSt]={{0}};

static vector<int> c2YRAccpStNum[6],c2YRRejeStNum[6];
static double c2YRChargeAverLevel[6];
static int c2YRNumClusTot[6][c2YRNumSt];
static double c2YRaM2[6][c2YRNumSt]={{0}},c2YRaM1[6][c2YRNumSt]={{0}},c2YRa0[6][c2YRNumSt]={{0}},c2YRaP1[6][c2YRNumSt]={{0}},c2YRaP2[6][c2YRNumSt]={{0}};
static double c2YRb[6][c2YRNumSt][c2YRNumSt]={{{0}}},c2YRc[6][c2YRNumSt]={{0}};
static double c2YRGainWeightAccp[6][c2YRNumSt]={{0}},c2YRGainWeight[6][c2YRNumSt]={{0}};

static vector<int> c2XRAccpStNum[6],c2XRRejeStNum[6];
static double c2XRChargeAverLevel[6];
static int c2XRNumClusTot[6][c2XRNumSt];
static double c2XRaM2[6][c2XRNumSt]={{0}},c2XRaM1[6][c2XRNumSt]={{0}},c2XRa0[6][c2XRNumSt]={{0}},c2XRaP1[6][c2XRNumSt]={{0}},c2XRaP2[6][c2XRNumSt]={{0}};
static double c2XRb[6][c2XRNumSt][c2XRNumSt]={{{0}}},c2XRc[6][c2XRNumSt]={{0}};
static double c2XRGainWeightAccp[6][c2XRNumSt]={{0}},c2XRGainWeight[6][c2XRNumSt]={{0}};

static vector<int> c2YLAccpStNum[6],c2YLRejeStNum[6];
static double c2YLChargeAverLevel[6];
static int c2YLNumClusTot[6][c2YLNumSt];
static double c2YLaM2[6][c2YLNumSt]={{0}},c2YLaM1[6][c2YLNumSt]={{0}},c2YLa0[6][c2YLNumSt]={{0}},c2YLaP1[6][c2YLNumSt]={{0}},c2YLaP2[6][c2YLNumSt]={{0}};
static double c2YLb[6][c2YLNumSt][c2YLNumSt]={{{0}}},c2YLc[6][c2YLNumSt]={{0}};
static double c2YLGainWeightAccp[6][c2YLNumSt]={{0}},c2YLGainWeight[6][c2YLNumSt]={{0}};

static vector<int> c2XLAccpStNum[6],c2XLRejeStNum[6];
static double c2XLChargeAverLevel[6];
static int c2XLNumClusTot[6][c2XLNumSt];
static double c2XLaM2[6][c2XLNumSt]={{0}},c2XLaM1[6][c2XLNumSt]={{0}},c2XLa0[6][c2XLNumSt]={{0}},c2XLaP1[6][c2XLNumSt]={{0}},c2XLaP2[6][c2XLNumSt]={{0}};
static double c2XLb[6][c2XLNumSt][c2XLNumSt]={{{0}}},c2XLc[6][c2XLNumSt]={{0}};
static double c2XLGainWeightAccp[6][c2XLNumSt]={{0}},c2XLGainWeight[6][c2XLNumSt]={{0}};


static double chargeTot;
static int clusTot;
Long_t Det_MWPC::done_amplifierGainCalibration(){
  for(int iGap=0;iGap<6;iGap++){
    for(int iSt=0;iSt<c3YRNumSt;iSt++){
      c3YRNumClusTot[iGap][iSt]=c3YRNumClusSing[iGap][iSt]+c3YRNumClusDoub[iGap][iSt]+c3YRNumClusTrip[iGap][iSt];
      if(c3YRNumClusTot[iGap][iSt]>statThreshold){
        c3YRaM2[iGap][iSt]=c3YRSumChargeTripM2[iGap][iSt]/c3YRNumClusTot[iGap][iSt];
        c3YRaM1[iGap][iSt]=(c3YRSumChargeDoubM1[iGap][iSt]+c3YRSumChargeTripM1[iGap][iSt])/c3YRNumClusTot[iGap][iSt];
        c3YRa0[iGap][iSt]=(c3YRSumChargeSing[iGap][iSt]+c3YRSumChargeDoub0[iGap][iSt]+c3YRSumChargeTrip0[iGap][iSt])/c3YRNumClusTot[iGap][iSt];
        c3YRaP1[iGap][iSt]=(c3YRSumChargeDoubP1[iGap][iSt]+c3YRSumChargeTripP1[iGap][iSt])/c3YRNumClusTot[iGap][iSt];
        c3YRaP2[iGap][iSt]=c3YRSumChargeTripP2[iGap][iSt]/c3YRNumClusTot[iGap][iSt];
        c3YRAccpStNum[iGap].push_back(iSt); 
      }
      else c3YRRejeStNum[iGap].push_back(iSt);
    }
  }
  for(int iGap=0;iGap<6;iGap++){
    if(c3YRAccpStNum[iGap].size()<=4){
      cout<<"Error: Number of accepted strips with satisfication of statistical threshold is equal to or less than 4 for gap "<<iGap+1<<" of c3YR!"<<endl;
      cout<<"Please reset statistical threshold."<<endl;
      return 0;
    }
    chargeTot=0;
    clusTot=0;
    for(UInt_t i=0;i<c3YRAccpStNum[iGap].size();i++){
      chargeTot+=(c3YRaM2[iGap][c3YRAccpStNum[iGap][i]]+c3YRaM1[iGap][c3YRAccpStNum[iGap][i]]+c3YRa0[iGap][c3YRAccpStNum[iGap][i]]+c3YRaP1[iGap][c3YRAccpStNum[iGap][i]]+c3YRaP2[iGap][c3YRAccpStNum[iGap][i]])*c3YRNumClusTot[iGap][c3YRAccpStNum[iGap][i]];
      clusTot+=c3YRNumClusTot[iGap][c3YRAccpStNum[iGap][i]];
    }
    c3YRChargeAverLevel[iGap]=chargeTot/clusTot;
  }    
  for(int iGap=0;iGap<6;iGap++){
    c3YRb[iGap][0][0]=c3YRa0[iGap][0];
    c3YRb[iGap][0][1]=c3YRaP1[iGap][0];
    c3YRb[iGap][0][2]=c3YRaP2[iGap][0];
    c3YRc[iGap][0]=c3YRChargeAverLevel[iGap];
    c3YRb[iGap][1][0]=c3YRaM1[iGap][1];
    c3YRb[iGap][1][1]=c3YRa0[iGap][1];
    c3YRb[iGap][1][2]=c3YRaP1[iGap][1];
    c3YRb[iGap][1][3]=c3YRaP2[iGap][1];
    c3YRc[iGap][1]=c3YRChargeAverLevel[iGap];
    for(UInt_t i=2;i<=c3YRNumSt-3;i++){
      c3YRb[iGap][i][i-2]=c3YRaM2[iGap][i];
      c3YRb[iGap][i][i-1]=c3YRaM1[iGap][i];
      c3YRb[iGap][i][i]=c3YRa0[iGap][i];
      c3YRb[iGap][i][i+1]=c3YRaP1[iGap][i];
      c3YRb[iGap][i][i+2]=c3YRaP2[iGap][i];
      c3YRc[iGap][i]=c3YRChargeAverLevel[iGap];
    }
    c3YRb[iGap][c3YRNumSt-2][c3YRNumSt-4]=c3YRaM2[iGap][c3YRNumSt-2];
    c3YRb[iGap][c3YRNumSt-2][c3YRNumSt-3]=c3YRaM1[iGap][c3YRNumSt-2];
    c3YRb[iGap][c3YRNumSt-2][c3YRNumSt-2]=c3YRa0[iGap][c3YRNumSt-2];
    c3YRb[iGap][c3YRNumSt-2][c3YRNumSt-1]=c3YRaP1[iGap][c3YRNumSt-2];
    c3YRc[iGap][c3YRNumSt-2]=c3YRChargeAverLevel[iGap];
    c3YRb[iGap][c3YRNumSt-1][c3YRNumSt-3]=c3YRaM2[iGap][c3YRNumSt-1];
    c3YRb[iGap][c3YRNumSt-1][c3YRNumSt-2]=c3YRaM1[iGap][c3YRNumSt-1];
    c3YRb[iGap][c3YRNumSt-1][c3YRNumSt-1]=c3YRa0[iGap][c3YRNumSt-1];
    c3YRc[iGap][c3YRNumSt-1]=c3YRChargeAverLevel[iGap];
  }
  for(int iGap=0;iGap<6;iGap++){
    TMatrix c3YRB(c3YRAccpStNum[iGap].size(),c3YRAccpStNum[iGap].size()), c3YRW(c3YRAccpStNum[iGap].size(),c3YRAccpStNum[iGap].size());
    for(UInt_t i=0;i<c3YRAccpStNum[iGap].size();i++){
      for(UInt_t j=0;j<c3YRAccpStNum[iGap].size();j++){
        c3YRB(i,j)=c3YRb[iGap][c3YRAccpStNum[iGap][i]][c3YRAccpStNum[iGap][j]];
        for(UInt_t k=0;k<c3YRRejeStNum[iGap].size();k++){
          c3YRB(i,j)+=c3YRb[iGap][c3YRAccpStNum[iGap][i]][c3YRRejeStNum[iGap][k]]/(c3YRNumSt-c3YRRejeStNum[iGap].size());
        }
      }
    }
    
    for(UInt_t j=0;j<c3YRAccpStNum[iGap].size();j++){
      c3YRW=c3YRB;
      for(UInt_t i=0;i<c3YRAccpStNum[iGap].size();i++){
        c3YRW(i,j)=c3YRc[iGap][i];
      }
      c3YRGainWeightAccp[iGap][c3YRAccpStNum[iGap][j]]=c3YRW.Determinant()/c3YRB.Determinant();
    }
  }
  double c3YRGainWeightAccpSum[6]={0},c3YRGainWeightAccpAver[6];
  for(int iGap=0;iGap<6;iGap++){
    for(UInt_t i=0;i<c3YRAccpStNum[iGap].size();i++){
      c3YRGainWeightAccpSum[iGap]+=c3YRGainWeightAccp[iGap][c3YRAccpStNum[iGap][i]]; 
    }
    if(c3YRAccpStNum[iGap].size()!=0) c3YRGainWeightAccpAver[iGap]=c3YRGainWeightAccpSum[iGap]/c3YRAccpStNum[iGap].size();
    else c3YRGainWeightAccpAver[iGap]=0;
  }
  for(int iGap=0;iGap<6;iGap++){
    if(c3YRAccpStNum[iGap][0]!=0){
      for(int j=0;j<c3YRAccpStNum[iGap][0];j++){
        c3YRGainWeight[iGap][j]=c3YRGainWeightAccpAver[iGap];
      }
    }
    for(UInt_t i=0;i<c3YRAccpStNum[iGap].size()-1;i++){
      c3YRGainWeight[iGap][c3YRAccpStNum[iGap][i]]=c3YRGainWeightAccp[iGap][c3YRAccpStNum[iGap][i]];
      c3YRGainWeight[iGap][c3YRAccpStNum[iGap][i+1]]=c3YRGainWeightAccp[iGap][c3YRAccpStNum[iGap][i+1]];
      for(int j=c3YRAccpStNum[iGap][i]+1;j<c3YRAccpStNum[iGap][i+1];j++){
        c3YRGainWeight[iGap][j]=c3YRGainWeightAccpAver[iGap];
      }
    }
    if(c3YRAccpStNum[iGap][c3YRAccpStNum[iGap].size()-1]!=c3YRNumSt-1){
      for(int j=c3YRAccpStNum[iGap][c3YRAccpStNum[iGap].size()-1]+1;j<c3YRNumSt;j++){
        c3YRGainWeight[iGap][j]=c3YRGainWeightAccpAver[iGap];
      }
    }
  }    
  ofstream outPutNumClusterC3YR("number-of-cluster-C3YR");
  for(int iGap=0;iGap<6;iGap++){
    for(int iSt=0;iSt<c3YRNumSt;iSt++){
      outPutNumClusterC3YR<<c3YRNumClusTot[iGap][iSt]<<"  ";
    }
    outPutNumClusterC3YR<<endl;
  }
  outPutNumClusterC3YR.close();

  ofstream outPutGainC3YR("gain-weight-C3YR");
  for(int iGap=0;iGap<6;iGap++){
    for(int i=0;i<c3YRNumSt;i++){
      outPutGainC3YR<<c3YRGainWeight[iGap][i]<<"  ";
    }
    outPutGainC3YR<<endl;
  }
  outPutGainC3YR.close();

  for(int iGap=0;iGap<6;iGap++){
    for(int iSt=0;iSt<c3XRNumSt;iSt++){
      c3XRNumClusTot[iGap][iSt]=c3XRNumClusSing[iGap][iSt]+c3XRNumClusDoub[iGap][iSt]+c3XRNumClusTrip[iGap][iSt];
      if(c3XRNumClusTot[iGap][iSt]>statThreshold){
        c3XRaM2[iGap][iSt]=c3XRSumChargeTripM2[iGap][iSt]/c3XRNumClusTot[iGap][iSt];
        c3XRaM1[iGap][iSt]=(c3XRSumChargeDoubM1[iGap][iSt]+c3XRSumChargeTripM1[iGap][iSt])/c3XRNumClusTot[iGap][iSt];
        c3XRa0[iGap][iSt]=(c3XRSumChargeSing[iGap][iSt]+c3XRSumChargeDoub0[iGap][iSt]+c3XRSumChargeTrip0[iGap][iSt])/c3XRNumClusTot[iGap][iSt];
        c3XRaP1[iGap][iSt]=(c3XRSumChargeDoubP1[iGap][iSt]+c3XRSumChargeTripP1[iGap][iSt])/c3XRNumClusTot[iGap][iSt];
        c3XRaP2[iGap][iSt]=c3XRSumChargeTripP2[iGap][iSt]/c3XRNumClusTot[iGap][iSt];
        c3XRAccpStNum[iGap].push_back(iSt);
      }
      else c3XRRejeStNum[iGap].push_back(iSt);
    }
  }
  for(int iGap=0;iGap<6;iGap++){
    if(c3XRAccpStNum[iGap].size()<=4){
      cout<<"Error: Number of accepted strips with satisfication of statistical threshold is equal to or less than 4 for gap "<<iGap+1<<" of c3XR!"<<endl;
      cout<<"Please reset statistical threshold."<<endl;
      return 0;
    }
    chargeTot=0;
    clusTot=0;
    for(UInt_t i=0;i<c3XRAccpStNum[iGap].size();i++){
      chargeTot+=(c3XRaM2[iGap][c3XRAccpStNum[iGap][i]]+c3XRaM1[iGap][c3XRAccpStNum[iGap][i]]+c3XRa0[iGap][c3XRAccpStNum[iGap][i]]+c3XRaP1[iGap][c3XRAccpStNum[iGap][i]]+c3XRaP2[iGap][c3XRAccpStNum[iGap][i]])*c3XRNumClusTot[iGap][c3XRAccpStNum[iGap][i]];
      clusTot+=c3XRNumClusTot[iGap][c3XRAccpStNum[iGap][i]];
    }
    c3XRChargeAverLevel[iGap]=chargeTot/clusTot;
  }
  for(int iGap=0;iGap<6;iGap++){
    c3XRb[iGap][0][0]=c3XRa0[iGap][0];
    c3XRb[iGap][0][1]=c3XRaP1[iGap][0];
    c3XRb[iGap][0][2]=c3XRaP2[iGap][0];
    c3XRc[iGap][0]=c3XRChargeAverLevel[iGap];
    c3XRb[iGap][1][0]=c3XRaM1[iGap][1];
    c3XRb[iGap][1][1]=c3XRa0[iGap][1];
    c3XRb[iGap][1][2]=c3XRaP1[iGap][1];
    c3XRb[iGap][1][3]=c3XRaP2[iGap][1];
    c3XRc[iGap][1]=c3XRChargeAverLevel[iGap];
    for(UInt_t i=2;i<=c3XRNumSt-3;i++){
      c3XRb[iGap][i][i-2]=c3XRaM2[iGap][i];
      c3XRb[iGap][i][i-1]=c3XRaM1[iGap][i];
      c3XRb[iGap][i][i]=c3XRa0[iGap][i];
      c3XRb[iGap][i][i+1]=c3XRaP1[iGap][i];
      c3XRb[iGap][i][i+2]=c3XRaP2[iGap][i];
      c3XRc[iGap][i]=c3XRChargeAverLevel[iGap];
    }
    c3XRb[iGap][c3XRNumSt-2][c3XRNumSt-4]=c3XRaM2[iGap][c3XRNumSt-2];
    c3XRb[iGap][c3XRNumSt-2][c3XRNumSt-3]=c3XRaM1[iGap][c3XRNumSt-2];
    c3XRb[iGap][c3XRNumSt-2][c3XRNumSt-2]=c3XRa0[iGap][c3XRNumSt-2];
    c3XRb[iGap][c3XRNumSt-2][c3XRNumSt-1]=c3XRaP1[iGap][c3XRNumSt-2];
    c3XRc[iGap][c3XRNumSt-2]=c3XRChargeAverLevel[iGap];
    c3XRb[iGap][c3XRNumSt-1][c3XRNumSt-3]=c3XRaM2[iGap][c3XRNumSt-1];
    c3XRb[iGap][c3XRNumSt-1][c3XRNumSt-2]=c3XRaM1[iGap][c3XRNumSt-1];
    c3XRb[iGap][c3XRNumSt-1][c3XRNumSt-1]=c3XRa0[iGap][c3XRNumSt-1];
    c3XRc[iGap][c3XRNumSt-1]=c3XRChargeAverLevel[iGap];
  }

  vector<int> c3XRRejeStNumCable[6][4];
  for(int iGap=0;iGap<6;iGap++){
    for(UInt_t i=0;i<c3XRRejeStNum[iGap].size();i++){
      for(int j=0;j<4;j++){
        if(c3XRRejeStNum[iGap][i]>=16*j && c3XRRejeStNum[iGap][i]<16*(j+1)) c3XRRejeStNumCable[iGap][j].push_back(c3XRRejeStNum[iGap][i]);
      }
    }
  }
  for(int iGap=0;iGap<6;iGap++){
    TMatrix c3XRB(c3XRAccpStNum[iGap].size(),c3XRAccpStNum[iGap].size()), c3XRW(c3XRAccpStNum[iGap].size(),c3XRAccpStNum[iGap].size());
    for(UInt_t i=0;i<c3XRAccpStNum[iGap].size();i++){
      for(UInt_t j=0;j<c3XRAccpStNum[iGap].size();j++){
        c3XRB(i,j)=c3XRb[iGap][c3XRAccpStNum[iGap][i]][c3XRAccpStNum[iGap][j]];
        for(int k=0;k<4;k++){
          for(UInt_t l=0;l<c3XRRejeStNumCable[iGap][k].size();l++){
            if(c3XRAccpStNum[iGap][j]>=16*k && c3XRAccpStNum[iGap][j]<16*(k+1)) c3XRB(i,j)+=c3XRb[iGap][c3XRAccpStNum[iGap][i]][c3XRRejeStNumCable[iGap][k][l]]/(16-c3XRRejeStNumCable[iGap][k].size()); 
          }
        }
      }
    }
    for(UInt_t j=0;j<c3XRAccpStNum[iGap].size();j++){
      c3XRW=c3XRB;
      for(UInt_t i=0;i<c3XRAccpStNum[iGap].size();i++){
        c3XRW(i,j)=c3XRc[iGap][i];
      }
      c3XRGainWeightAccp[iGap][c3XRAccpStNum[iGap][j]]=c3XRW.Determinant()/c3XRB.Determinant();
    }
  }
  double c3XRGainWeightAccpSum[6][4]={{0}},c3XRGainWeightAccpAver[6][4];
  int c3XRAccpNumStCable[6][4]={{0}};
  for(int iGap=0;iGap<6;iGap++){
    for(UInt_t i=0;i<c3XRAccpStNum[iGap].size();i++){
      for(int j=0;j<4;j++){
        if(c3XRAccpStNum[iGap][i]>=16*j && c3XRAccpStNum[iGap][i]<16*(j+1)){
          c3XRGainWeightAccpSum[iGap][j]+=c3XRGainWeightAccp[iGap][c3XRAccpStNum[iGap][i]];
          c3XRAccpNumStCable[iGap][j]++;
          break;
        }
      }
    }
    for(int j=0;j<4;j++){
      if(c3XRAccpNumStCable[iGap][j]!=0) c3XRGainWeightAccpAver[iGap][j]=c3XRGainWeightAccpSum[iGap][j]/c3XRAccpNumStCable[iGap][j];
      else c3XRGainWeightAccpAver[iGap][j]=1;
    }
  }
  for(int iGap=0;iGap<6;iGap++){
    if(c3XRAccpStNum[iGap][0]!=0){
      for(int j=0;j<c3XRAccpStNum[iGap][0];j++){
        for(int k=0;k<4;k++){
          if(j>=16*k && j<16*(k+1)){
            c3XRGainWeight[iGap][j]=c3XRGainWeightAccpAver[iGap][k];
            break;
          }
        }
      }
    }
    for(UInt_t i=0;i<c3XRAccpStNum[iGap].size()-1;i++){
      c3XRGainWeight[iGap][c3XRAccpStNum[iGap][i]]=c3XRGainWeightAccp[iGap][c3XRAccpStNum[iGap][i]];
      c3XRGainWeight[iGap][c3XRAccpStNum[iGap][i+1]]=c3XRGainWeightAccp[iGap][c3XRAccpStNum[iGap][i+1]];
      for(int j=c3XRAccpStNum[iGap][i]+1;j<c3XRAccpStNum[iGap][i+1];j++){
        for(int k=0;k<4;k++){
          if(j>=16*k && j<16*(k+1)){
            c3XRGainWeight[iGap][j]=c3XRGainWeightAccpAver[iGap][k];
            break;
          }
        }
      }
    }
    if(c3XRAccpStNum[iGap][c3XRAccpStNum[iGap].size()-1]!=c3XRNumSt-1){
      for(int j=c3XRAccpStNum[iGap][c3XRAccpStNum[iGap].size()-1]+1;j<c3XRNumSt;j++){
        for(int k=0;k<4;k++){
          if(j>=16*k && j<16*(k+1)){
            c3XRGainWeight[iGap][j]=c3XRGainWeightAccpAver[iGap][k];
            break;
          }
        }
      }
    }
  }
  ofstream outPutNumClusterC3XR("number-of-cluster-C3XR");
  for(int iGap=0;iGap<6;iGap++){
    for(int iSt=0;iSt<c3XRNumSt;iSt++){
      outPutNumClusterC3XR<<c3XRNumClusTot[iGap][iSt]<<"  ";
    }
    outPutNumClusterC3XR<<endl;
  }
  outPutNumClusterC3XR.close();

  ofstream outPutGainC3XR("gain-weight-C3XR");
  for(int iGap=0;iGap<6;iGap++){
    for(int i=0;i<c3XRNumSt;i++){
      outPutGainC3XR<<c3XRGainWeight[iGap][i]<<"  ";
    }
    outPutGainC3XR<<endl;
  }
  outPutGainC3XR.close();

  ofstream outPutNormXYC3R("normalized-factor-C3R");
  for(int iGap=0;iGap<6;iGap++){
    outPutNormXYC3R<<c3XRChargeAverLevel[iGap]/c3YRChargeAverLevel[iGap]<<"  ";
  }
  outPutNormXYC3R.close();

  for(int iGap=0;iGap<6;iGap++){
    for(int iSt=0;iSt<c3YLNumSt;iSt++){
      c3YLNumClusTot[iGap][iSt]=c3YLNumClusSing[iGap][iSt]+c3YLNumClusDoub[iGap][iSt]+c3YLNumClusTrip[iGap][iSt];
      if(c3YLNumClusTot[iGap][iSt]>statThreshold){
        c3YLaM2[iGap][iSt]=c3YLSumChargeTripM2[iGap][iSt]/c3YLNumClusTot[iGap][iSt];
        c3YLaM1[iGap][iSt]=(c3YLSumChargeDoubM1[iGap][iSt]+c3YLSumChargeTripM1[iGap][iSt])/c3YLNumClusTot[iGap][iSt];
        c3YLa0[iGap][iSt]=(c3YLSumChargeSing[iGap][iSt]+c3YLSumChargeDoub0[iGap][iSt]+c3YLSumChargeTrip0[iGap][iSt])/c3YLNumClusTot[iGap][iSt];
        c3YLaP1[iGap][iSt]=(c3YLSumChargeDoubP1[iGap][iSt]+c3YLSumChargeTripP1[iGap][iSt])/c3YLNumClusTot[iGap][iSt];
        c3YLaP2[iGap][iSt]=c3YLSumChargeTripP2[iGap][iSt]/c3YLNumClusTot[iGap][iSt];
        c3YLAccpStNum[iGap].push_back(iSt);
      }
      else c3YLRejeStNum[iGap].push_back(iSt);
    }
  }
  for(int iGap=0;iGap<6;iGap++){
    if(c3YLAccpStNum[iGap].size()<=4){
      cout<<"Error: Number of accepted strips with satisfication of statistical threshold is equal to or less than 4 for gap "<<iGap+1<<" of c3YL!"<<endl;
      cout<<"Please reset statistical threshold."<<endl;
      return 0;
    }
    chargeTot=0;
    clusTot=0;
    for(UInt_t i=0;i<c3YLAccpStNum[iGap].size();i++){
      chargeTot+=(c3YLaM2[iGap][c3YLAccpStNum[iGap][i]]+c3YLaM1[iGap][c3YLAccpStNum[iGap][i]]+c3YLa0[iGap][c3YLAccpStNum[iGap][i]]+c3YLaP1[iGap][c3YLAccpStNum[iGap][i]]+c3YLaP2[iGap][c3YLAccpStNum[iGap][i]])*c3YLNumClusTot[iGap][c3YLAccpStNum[iGap][i]];
      clusTot+=c3YLNumClusTot[iGap][c3YLAccpStNum[iGap][i]];
    }
    c3YLChargeAverLevel[iGap]=chargeTot/clusTot;
  }
  for(int iGap=0;iGap<6;iGap++){
    c3YLb[iGap][0][0]=c3YLa0[iGap][0];
    c3YLb[iGap][0][1]=c3YLaP1[iGap][0];
    c3YLb[iGap][0][2]=c3YLaP2[iGap][0];
    c3YLc[iGap][0]=c3YLChargeAverLevel[iGap];
    c3YLb[iGap][1][0]=c3YLaM1[iGap][1];
    c3YLb[iGap][1][1]=c3YLa0[iGap][1];
    c3YLb[iGap][1][2]=c3YLaP1[iGap][1];
    c3YLb[iGap][1][3]=c3YLaP2[iGap][1];
    c3YLc[iGap][1]=c3YLChargeAverLevel[iGap];
    for(UInt_t i=2;i<=c3YLNumSt-3;i++){
      c3YLb[iGap][i][i-2]=c3YLaM2[iGap][i];
      c3YLb[iGap][i][i-1]=c3YLaM1[iGap][i];
      c3YLb[iGap][i][i]=c3YLa0[iGap][i];
      c3YLb[iGap][i][i+1]=c3YLaP1[iGap][i];
      c3YLb[iGap][i][i+2]=c3YLaP2[iGap][i];
      c3YLc[iGap][i]=c3YLChargeAverLevel[iGap];
    }
    c3YLb[iGap][c3YLNumSt-2][c3YLNumSt-4]=c3YLaM2[iGap][c3YLNumSt-2];
    c3YLb[iGap][c3YLNumSt-2][c3YLNumSt-3]=c3YLaM1[iGap][c3YLNumSt-2];
    c3YLb[iGap][c3YLNumSt-2][c3YLNumSt-2]=c3YLa0[iGap][c3YLNumSt-2];
    c3YLb[iGap][c3YLNumSt-2][c3YLNumSt-1]=c3YLaP1[iGap][c3YLNumSt-2];
    c3YLc[iGap][c3YLNumSt-2]=c3YLChargeAverLevel[iGap];
    c3YLb[iGap][c3YLNumSt-1][c3YLNumSt-3]=c3YLaM2[iGap][c3YLNumSt-1];
    c3YLb[iGap][c3YLNumSt-1][c3YLNumSt-2]=c3YLaM1[iGap][c3YLNumSt-1];
    c3YLb[iGap][c3YLNumSt-1][c3YLNumSt-1]=c3YLa0[iGap][c3YLNumSt-1];
    c3YLc[iGap][c3YLNumSt-1]=c3YLChargeAverLevel[iGap];
  }
  for(int iGap=0;iGap<6;iGap++){
    TMatrix c3YLB(c3YLAccpStNum[iGap].size(),c3YLAccpStNum[iGap].size()), c3YLW(c3YLAccpStNum[iGap].size(),c3YLAccpStNum[iGap].size());
    for(UInt_t i=0;i<c3YLAccpStNum[iGap].size();i++){
      for(UInt_t j=0;j<c3YLAccpStNum[iGap].size();j++){
        c3YLB(i,j)=c3YLb[iGap][c3YLAccpStNum[iGap][i]][c3YLAccpStNum[iGap][j]];
        for(UInt_t k=0;k<c3YLRejeStNum[iGap].size();k++){
          c3YLB(i,j)+=c3YLb[iGap][c3YLAccpStNum[iGap][i]][c3YLRejeStNum[iGap][k]]/(c3YLNumSt-c3YLRejeStNum[iGap].size());
        }
      }
    }

    for(UInt_t j=0;j<c3YLAccpStNum[iGap].size();j++){
      c3YLW=c3YLB;
      for(UInt_t i=0;i<c3YLAccpStNum[iGap].size();i++){
        c3YLW(i,j)=c3YLc[iGap][i];
      }
      c3YLGainWeightAccp[iGap][c3YLAccpStNum[iGap][j]]=c3YLW.Determinant()/c3YLB.Determinant();
    }
  }
  double c3YLGainWeightAccpSum[6]={0},c3YLGainWeightAccpAver[6];
  for(int iGap=0;iGap<6;iGap++){
    for(UInt_t i=0;i<c3YLAccpStNum[iGap].size();i++){
      c3YLGainWeightAccpSum[iGap]+=c3YLGainWeightAccp[iGap][c3YLAccpStNum[iGap][i]];
    }
    if(c3YLAccpStNum[iGap].size()!=0) c3YLGainWeightAccpAver[iGap]=c3YLGainWeightAccpSum[iGap]/c3YLAccpStNum[iGap].size();
    else c3YLGainWeightAccpAver[iGap]=0;
  }
  for(int iGap=0;iGap<6;iGap++){
    if(c3YLAccpStNum[iGap][0]!=0){
      for(int j=0;j<c3YLAccpStNum[iGap][0];j++){
        c3YLGainWeight[iGap][j]=c3YLGainWeightAccpAver[iGap];
      }
    }
    for(UInt_t i=0;i<c3YLAccpStNum[iGap].size()-1;i++){
      c3YLGainWeight[iGap][c3YLAccpStNum[iGap][i]]=c3YLGainWeightAccp[iGap][c3YLAccpStNum[iGap][i]];
      c3YLGainWeight[iGap][c3YLAccpStNum[iGap][i+1]]=c3YLGainWeightAccp[iGap][c3YLAccpStNum[iGap][i+1]];
      for(int j=c3YLAccpStNum[iGap][i]+1;j<c3YLAccpStNum[iGap][i+1];j++){
        c3YLGainWeight[iGap][j]=c3YLGainWeightAccpAver[iGap];
      }
    }
    if(c3YLAccpStNum[iGap][c3YLAccpStNum[iGap].size()-1]!=c3YLNumSt-1){
      for(int j=c3YLAccpStNum[iGap][c3YLAccpStNum[iGap].size()-1]+1;j<c3YLNumSt;j++){
        c3YLGainWeight[iGap][j]=c3YLGainWeightAccpAver[iGap];
      }
    }
  }
  ofstream outPutNumClusterC3YL("number-of-cluster-C3YL");
  for(int iGap=0;iGap<6;iGap++){
    for(int iSt=0;iSt<c3YLNumSt;iSt++){
      outPutNumClusterC3YL<<c3YLNumClusTot[iGap][iSt]<<"  ";
    }
    outPutNumClusterC3YL<<endl;
  }
  outPutNumClusterC3YL.close();

  ofstream outPutGainC3YL("gain-weight-C3YL");
  for(int iGap=0;iGap<6;iGap++){
    for(int i=0;i<c3YLNumSt;i++){
      outPutGainC3YL<<c3YLGainWeight[iGap][i]<<"  ";
    }
    outPutGainC3YL<<endl;
  }
  outPutGainC3YL.close();

  for(int iGap=0;iGap<6;iGap++){
    for(int iSt=0;iSt<c3XLNumSt;iSt++){
      c3XLNumClusTot[iGap][iSt]=c3XLNumClusSing[iGap][iSt]+c3XLNumClusDoub[iGap][iSt]+c3XLNumClusTrip[iGap][iSt];
      if(c3XLNumClusTot[iGap][iSt]>statThreshold){
        c3XLaM2[iGap][iSt]=c3XLSumChargeTripM2[iGap][iSt]/c3XLNumClusTot[iGap][iSt];
        c3XLaM1[iGap][iSt]=(c3XLSumChargeDoubM1[iGap][iSt]+c3XLSumChargeTripM1[iGap][iSt])/c3XLNumClusTot[iGap][iSt];
        c3XLa0[iGap][iSt]=(c3XLSumChargeSing[iGap][iSt]+c3XLSumChargeDoub0[iGap][iSt]+c3XLSumChargeTrip0[iGap][iSt])/c3XLNumClusTot[iGap][iSt];
        c3XLaP1[iGap][iSt]=(c3XLSumChargeDoubP1[iGap][iSt]+c3XLSumChargeTripP1[iGap][iSt])/c3XLNumClusTot[iGap][iSt];
        c3XLaP2[iGap][iSt]=c3XLSumChargeTripP2[iGap][iSt]/c3XLNumClusTot[iGap][iSt];
        c3XLAccpStNum[iGap].push_back(iSt);
      }
      else c3XLRejeStNum[iGap].push_back(iSt);
    }
  }
  for(int iGap=0;iGap<6;iGap++){
    if(c3XLAccpStNum[iGap].size()<=4){
      cout<<"Error: Number of accepted strips with satisfication of statistical threshold is equal to or less than 4 for gap "<<iGap+1<<" of c3XL!"<<endl;
      cout<<"Please reset statistical threshold."<<endl;
      return 0;
    }
    chargeTot=0;
    clusTot=0;
    for(UInt_t i=0;i<c3XLAccpStNum[iGap].size();i++){
      chargeTot+=(c3XLaM2[iGap][c3XLAccpStNum[iGap][i]]+c3XLaM1[iGap][c3XLAccpStNum[iGap][i]]+c3XLa0[iGap][c3XLAccpStNum[iGap][i]]+c3XLaP1[iGap][c3XLAccpStNum[iGap][i]]+c3XLaP2[iGap][c3XLAccpStNum[iGap][i]])*c3XLNumClusTot[iGap][c3XLAccpStNum[iGap][i]];
      clusTot+=c3XLNumClusTot[iGap][c3XLAccpStNum[iGap][i]];
    }
    c3XLChargeAverLevel[iGap]=chargeTot/clusTot;
  }
  for(int iGap=0;iGap<6;iGap++){
    c3XLb[iGap][0][0]=c3XLa0[iGap][0];
    c3XLb[iGap][0][1]=c3XLaP1[iGap][0];
    c3XLb[iGap][0][2]=c3XLaP2[iGap][0];
    c3XLc[iGap][0]=c3XLChargeAverLevel[iGap];
    c3XLb[iGap][1][0]=c3XLaM1[iGap][1];
    c3XLb[iGap][1][1]=c3XLa0[iGap][1];
    c3XLb[iGap][1][2]=c3XLaP1[iGap][1];
    c3XLb[iGap][1][3]=c3XLaP2[iGap][1];
    c3XLc[iGap][1]=c3XLChargeAverLevel[iGap];
    for(UInt_t i=2;i<=c3XLNumSt-3;i++){
      c3XLb[iGap][i][i-2]=c3XLaM2[iGap][i];
      c3XLb[iGap][i][i-1]=c3XLaM1[iGap][i];
      c3XLb[iGap][i][i]=c3XLa0[iGap][i];
      c3XLb[iGap][i][i+1]=c3XLaP1[iGap][i];
      c3XLb[iGap][i][i+2]=c3XLaP2[iGap][i];
      c3XLc[iGap][i]=c3XLChargeAverLevel[iGap];
    }
    c3XLb[iGap][c3XLNumSt-2][c3XLNumSt-4]=c3XLaM2[iGap][c3XLNumSt-2];
    c3XLb[iGap][c3XLNumSt-2][c3XLNumSt-3]=c3XLaM1[iGap][c3XLNumSt-2];
    c3XLb[iGap][c3XLNumSt-2][c3XLNumSt-2]=c3XLa0[iGap][c3XLNumSt-2];
    c3XLb[iGap][c3XLNumSt-2][c3XLNumSt-1]=c3XLaP1[iGap][c3XLNumSt-2];
    c3XLc[iGap][c3XLNumSt-2]=c3XLChargeAverLevel[iGap];
    c3XLb[iGap][c3XLNumSt-1][c3XLNumSt-3]=c3XLaM2[iGap][c3XLNumSt-1];
    c3XLb[iGap][c3XLNumSt-1][c3XLNumSt-2]=c3XLaM1[iGap][c3XLNumSt-1];
    c3XLb[iGap][c3XLNumSt-1][c3XLNumSt-1]=c3XLa0[iGap][c3XLNumSt-1];
    c3XLc[iGap][c3XLNumSt-1]=c3XLChargeAverLevel[iGap];
  }

  vector<int> c3XLRejeStNumCable[6][4];
  for(int iGap=0;iGap<6;iGap++){
    for(UInt_t i=0;i<c3XLRejeStNum[iGap].size();i++){
      for(int j=0;j<4;j++){
        if(c3XLRejeStNum[iGap][i]>=16*j && c3XLRejeStNum[iGap][i]<16*(j+1)) c3XLRejeStNumCable[iGap][j].push_back(c3XLRejeStNum[iGap][i]);
      }
    }
  }
  for(int iGap=0;iGap<6;iGap++){
    TMatrix c3XLB(c3XLAccpStNum[iGap].size(),c3XLAccpStNum[iGap].size()), c3XLW(c3XLAccpStNum[iGap].size(),c3XLAccpStNum[iGap].size());
    for(UInt_t i=0;i<c3XLAccpStNum[iGap].size();i++){
      for(UInt_t j=0;j<c3XLAccpStNum[iGap].size();j++){
        c3XLB(i,j)=c3XLb[iGap][c3XLAccpStNum[iGap][i]][c3XLAccpStNum[iGap][j]];
        for(int k=0;k<4;k++){
          for(UInt_t l=0;l<c3XLRejeStNumCable[iGap][k].size();l++){
            if(c3XLAccpStNum[iGap][j]>=16*k && c3XLAccpStNum[iGap][j]<16*(k+1)) c3XLB(i,j)+=c3XLb[iGap][c3XLAccpStNum[iGap][i]][c3XLRejeStNumCable[iGap][k][l]]/(16-c3XLRejeStNumCable[iGap][k].size());
          }
        }
      }
    }
    for(UInt_t j=0;j<c3XLAccpStNum[iGap].size();j++){
      c3XLW=c3XLB;
      for(UInt_t i=0;i<c3XLAccpStNum[iGap].size();i++){
        c3XLW(i,j)=c3XLc[iGap][i];
      }
      c3XLGainWeightAccp[iGap][c3XLAccpStNum[iGap][j]]=c3XLW.Determinant()/c3XLB.Determinant();
    }
  }
  double c3XLGainWeightAccpSum[6][4]={{0}},c3XLGainWeightAccpAver[6][4];
  int c3XLAccpNumStCable[6][4]={{0}};
  for(int iGap=0;iGap<6;iGap++){
    for(UInt_t i=0;i<c3XLAccpStNum[iGap].size();i++){
      for(int j=0;j<4;j++){
        if(c3XLAccpStNum[iGap][i]>=16*j && c3XLAccpStNum[iGap][i]<16*(j+1)){
          c3XLGainWeightAccpSum[iGap][j]+=c3XLGainWeightAccp[iGap][c3XLAccpStNum[iGap][i]];
          c3XLAccpNumStCable[iGap][j]++;
          break;
        }
      }
    }
    for(int j=0;j<4;j++){
      if(c3XLAccpNumStCable[iGap][j]!=0) c3XLGainWeightAccpAver[iGap][j]=c3XLGainWeightAccpSum[iGap][j]/c3XLAccpNumStCable[iGap][j];
      else c3XLGainWeightAccpAver[iGap][j]=1;
    }
  }
  for(int iGap=0;iGap<6;iGap++){
    if(c3XLAccpStNum[iGap][0]!=0){
      for(int j=0;j<c3XLAccpStNum[iGap][0];j++){
        for(int k=0;k<4;k++){
          if(j>=16*k && j<16*(k+1)){
            c3XLGainWeight[iGap][j]=c3XLGainWeightAccpAver[iGap][k];
            break;
          }
        }
      }
    }
    for(UInt_t i=0;i<c3XLAccpStNum[iGap].size()-1;i++){
      c3XLGainWeight[iGap][c3XLAccpStNum[iGap][i]]=c3XLGainWeightAccp[iGap][c3XLAccpStNum[iGap][i]];
      c3XLGainWeight[iGap][c3XLAccpStNum[iGap][i+1]]=c3XLGainWeightAccp[iGap][c3XLAccpStNum[iGap][i+1]];
      for(int j=c3XLAccpStNum[iGap][i]+1;j<c3XLAccpStNum[iGap][i+1];j++){
        for(int k=0;k<4;k++){
          if(j>=16*k && j<16*(k+1)){
            c3XLGainWeight[iGap][j]=c3XLGainWeightAccpAver[iGap][k];
            break;
          }
        }
      }
    }
    if(c3XLAccpStNum[iGap][c3XLAccpStNum[iGap].size()-1]!=c3XLNumSt-1){
      for(int j=c3XLAccpStNum[iGap][c3XLAccpStNum[iGap].size()-1]+1;j<c3XLNumSt;j++){
        for(int k=0;k<4;k++){
          if(j>=16*k && j<16*(k+1)){
            c3XLGainWeight[iGap][j]=c3XLGainWeightAccpAver[iGap][k];
            break;
          }
        }
      }
    }
  }
  ofstream outPutNumClusterC3XL("number-of-cluster-C3XL");
  for(int iGap=0;iGap<6;iGap++){
    for(int iSt=0;iSt<c3XLNumSt;iSt++){
      outPutNumClusterC3XL<<c3XLNumClusTot[iGap][iSt]<<"  ";
    }
    outPutNumClusterC3XL<<endl;
  }
  outPutNumClusterC3XL.close();

  ofstream outPutGainC3XL("gain-weight-C3XL");
  for(int iGap=0;iGap<6;iGap++){
    for(int i=0;i<c3XLNumSt;i++){
      outPutGainC3XL<<c3XLGainWeight[iGap][i]<<"  ";
    }
    outPutGainC3XL<<endl;
  }
  outPutGainC3XL.close();
  ofstream outPutNormXYC3L("normalized-factor-C3L");
  for(int iGap=0;iGap<6;iGap++){
    outPutNormXYC3L<<c3XLChargeAverLevel[iGap]/c3YLChargeAverLevel[iGap]<<"  ";
  }
  outPutNormXYC3L.close();

  for(int iGap=0;iGap<6;iGap++){
    for(int iSt=0;iSt<c4YRNumSt;iSt++){
      c4YRNumClusTot[iGap][iSt]=c4YRNumClusSing[iGap][iSt]+c4YRNumClusDoub[iGap][iSt]+c4YRNumClusTrip[iGap][iSt];
      if(c4YRNumClusTot[iGap][iSt]>statThreshold){
        c4YRaM2[iGap][iSt]=c4YRSumChargeTripM2[iGap][iSt]/c4YRNumClusTot[iGap][iSt];
        c4YRaM1[iGap][iSt]=(c4YRSumChargeDoubM1[iGap][iSt]+c4YRSumChargeTripM1[iGap][iSt])/c4YRNumClusTot[iGap][iSt];
        c4YRa0[iGap][iSt]=(c4YRSumChargeSing[iGap][iSt]+c4YRSumChargeDoub0[iGap][iSt]+c4YRSumChargeTrip0[iGap][iSt])/c4YRNumClusTot[iGap][iSt];
        c4YRaP1[iGap][iSt]=(c4YRSumChargeDoubP1[iGap][iSt]+c4YRSumChargeTripP1[iGap][iSt])/c4YRNumClusTot[iGap][iSt];
        c4YRaP2[iGap][iSt]=c4YRSumChargeTripP2[iGap][iSt]/c4YRNumClusTot[iGap][iSt];
        c4YRAccpStNum[iGap].push_back(iSt);
      }
      else c4YRRejeStNum[iGap].push_back(iSt);
    }
  }
  for(int iGap=0;iGap<6;iGap++){
    if(c4YRAccpStNum[iGap].size()<=4){
      cout<<"Error: Number of accepted strips with satisfication of statistical threshold is equal to or less than 4 for gap "<<iGap+1<<" of c4YR!"<<endl;
      cout<<"Please reset statistical threshold."<<endl;
      return 0;
    }
    chargeTot=0;
    clusTot=0;
    for(UInt_t i=0;i<c4YRAccpStNum[iGap].size();i++){
      chargeTot+=(c4YRaM2[iGap][c4YRAccpStNum[iGap][i]]+c4YRaM1[iGap][c4YRAccpStNum[iGap][i]]+c4YRa0[iGap][c4YRAccpStNum[iGap][i]]+c4YRaP1[iGap][c4YRAccpStNum[iGap][i]]+c4YRaP2[iGap][c4YRAccpStNum[iGap][i]])*c4YRNumClusTot[iGap][c4YRAccpStNum[iGap][i]];
      clusTot+=c4YRNumClusTot[iGap][c4YRAccpStNum[iGap][i]];
    }
    c4YRChargeAverLevel[iGap]=chargeTot/clusTot;
  }
  for(int iGap=0;iGap<6;iGap++){
    c4YRb[iGap][0][0]=c4YRa0[iGap][0];
    c4YRb[iGap][0][1]=c4YRaP1[iGap][0];
    c4YRb[iGap][0][2]=c4YRaP2[iGap][0];
    c4YRc[iGap][0]=c4YRChargeAverLevel[iGap];
    c4YRb[iGap][1][0]=c4YRaM1[iGap][1];
    c4YRb[iGap][1][1]=c4YRa0[iGap][1];
    c4YRb[iGap][1][2]=c4YRaP1[iGap][1];
    c4YRb[iGap][1][3]=c4YRaP2[iGap][1];
    c4YRc[iGap][1]=c4YRChargeAverLevel[iGap];
    for(UInt_t i=2;i<=c4YRNumSt-3;i++){
      c4YRb[iGap][i][i-2]=c4YRaM2[iGap][i];
      c4YRb[iGap][i][i-1]=c4YRaM1[iGap][i];
      c4YRb[iGap][i][i]=c4YRa0[iGap][i];
      c4YRb[iGap][i][i+1]=c4YRaP1[iGap][i];
      c4YRb[iGap][i][i+2]=c4YRaP2[iGap][i];
      c4YRc[iGap][i]=c4YRChargeAverLevel[iGap];
    }
    c4YRb[iGap][c4YRNumSt-2][c4YRNumSt-4]=c4YRaM2[iGap][c4YRNumSt-2];
    c4YRb[iGap][c4YRNumSt-2][c4YRNumSt-3]=c4YRaM1[iGap][c4YRNumSt-2];
    c4YRb[iGap][c4YRNumSt-2][c4YRNumSt-2]=c4YRa0[iGap][c4YRNumSt-2];
    c4YRb[iGap][c4YRNumSt-2][c4YRNumSt-1]=c4YRaP1[iGap][c4YRNumSt-2];
    c4YRc[iGap][c4YRNumSt-2]=c4YRChargeAverLevel[iGap];
    c4YRb[iGap][c4YRNumSt-1][c4YRNumSt-3]=c4YRaM2[iGap][c4YRNumSt-1];
    c4YRb[iGap][c4YRNumSt-1][c4YRNumSt-2]=c4YRaM1[iGap][c4YRNumSt-1];
    c4YRb[iGap][c4YRNumSt-1][c4YRNumSt-1]=c4YRa0[iGap][c4YRNumSt-1];
    c4YRc[iGap][c4YRNumSt-1]=c4YRChargeAverLevel[iGap];
  }
  for(int iGap=0;iGap<6;iGap++){
    TMatrix c4YRB(c4YRAccpStNum[iGap].size(),c4YRAccpStNum[iGap].size()), c4YRW(c4YRAccpStNum[iGap].size(),c4YRAccpStNum[iGap].size());
    for(UInt_t i=0;i<c4YRAccpStNum[iGap].size();i++){
      for(UInt_t j=0;j<c4YRAccpStNum[iGap].size();j++){
        c4YRB(i,j)=c4YRb[iGap][c4YRAccpStNum[iGap][i]][c4YRAccpStNum[iGap][j]];
        for(UInt_t k=0;k<c4YRRejeStNum[iGap].size();k++){
          c4YRB(i,j)+=c4YRb[iGap][c4YRAccpStNum[iGap][i]][c4YRRejeStNum[iGap][k]]/(c4YRNumSt-c4YRRejeStNum[iGap].size());
        }
      }
    }

    for(UInt_t j=0;j<c4YRAccpStNum[iGap].size();j++){
      c4YRW=c4YRB;
      for(UInt_t i=0;i<c4YRAccpStNum[iGap].size();i++){
        c4YRW(i,j)=c4YRc[iGap][i];
      }
      c4YRGainWeightAccp[iGap][c4YRAccpStNum[iGap][j]]=c4YRW.Determinant()/c4YRB.Determinant();
    }
  }
  double c4YRGainWeightAccpSum[6]={0},c4YRGainWeightAccpAver[6];
  for(int iGap=0;iGap<6;iGap++){
    for(UInt_t i=0;i<c4YRAccpStNum[iGap].size();i++){
      c4YRGainWeightAccpSum[iGap]+=c4YRGainWeightAccp[iGap][c4YRAccpStNum[iGap][i]];
    }
    if(c4YRAccpStNum[iGap].size()!=0) c4YRGainWeightAccpAver[iGap]=c4YRGainWeightAccpSum[iGap]/c4YRAccpStNum[iGap].size();
    else c4YRGainWeightAccpAver[iGap]=0;
  }
  for(int iGap=0;iGap<6;iGap++){
    if(c4YRAccpStNum[iGap][0]!=0){
      for(int j=0;j<c4YRAccpStNum[iGap][0];j++){
        c4YRGainWeight[iGap][j]=c4YRGainWeightAccpAver[iGap];
      }
    }
    for(UInt_t i=0;i<c4YRAccpStNum[iGap].size()-1;i++){
      c4YRGainWeight[iGap][c4YRAccpStNum[iGap][i]]=c4YRGainWeightAccp[iGap][c4YRAccpStNum[iGap][i]];
      c4YRGainWeight[iGap][c4YRAccpStNum[iGap][i+1]]=c4YRGainWeightAccp[iGap][c4YRAccpStNum[iGap][i+1]];
      for(int j=c4YRAccpStNum[iGap][i]+1;j<c4YRAccpStNum[iGap][i+1];j++){
        c4YRGainWeight[iGap][j]=c4YRGainWeightAccpAver[iGap];
      }
    }
    if(c4YRAccpStNum[iGap][c4YRAccpStNum[iGap].size()-1]!=c4YRNumSt-1){
      for(int j=c4YRAccpStNum[iGap][c4YRAccpStNum[iGap].size()-1]+1;j<c4YRNumSt;j++){
        c4YRGainWeight[iGap][j]=c4YRGainWeightAccpAver[iGap];
      }
    }
  }
  ofstream outPutNumClusterC4YR("number-of-cluster-C4YR");
  for(int iGap=0;iGap<6;iGap++){
    for(int iSt=0;iSt<c4YRNumSt;iSt++){
      outPutNumClusterC4YR<<c4YRNumClusTot[iGap][iSt]<<"  ";
    }
    outPutNumClusterC4YR<<endl;
  }
  outPutNumClusterC4YR.close();

  ofstream outPutGainC4YR("gain-weight-C4YR");
  for(int iGap=0;iGap<6;iGap++){
    for(int i=0;i<c4YRNumSt;i++){
      outPutGainC4YR<<c4YRGainWeight[iGap][i]<<"  ";
    }
    outPutGainC4YR<<endl;
  }
  outPutGainC4YR.close();

  for(int iGap=0;iGap<6;iGap++){
    for(int iSt=0;iSt<c4XRNumSt;iSt++){
      c4XRNumClusTot[iGap][iSt]=c4XRNumClusSing[iGap][iSt]+c4XRNumClusDoub[iGap][iSt]+c4XRNumClusTrip[iGap][iSt];
      if(c4XRNumClusTot[iGap][iSt]>statThreshold){
        c4XRaM2[iGap][iSt]=c4XRSumChargeTripM2[iGap][iSt]/c4XRNumClusTot[iGap][iSt];
        c4XRaM1[iGap][iSt]=(c4XRSumChargeDoubM1[iGap][iSt]+c4XRSumChargeTripM1[iGap][iSt])/c4XRNumClusTot[iGap][iSt];
        c4XRa0[iGap][iSt]=(c4XRSumChargeSing[iGap][iSt]+c4XRSumChargeDoub0[iGap][iSt]+c4XRSumChargeTrip0[iGap][iSt])/c4XRNumClusTot[iGap][iSt];
        c4XRaP1[iGap][iSt]=(c4XRSumChargeDoubP1[iGap][iSt]+c4XRSumChargeTripP1[iGap][iSt])/c4XRNumClusTot[iGap][iSt];
        c4XRaP2[iGap][iSt]=c4XRSumChargeTripP2[iGap][iSt]/c4XRNumClusTot[iGap][iSt];
        c4XRAccpStNum[iGap].push_back(iSt);
      }
      else c4XRRejeStNum[iGap].push_back(iSt);
    }
  }
  for(int iGap=0;iGap<6;iGap++){
    if(c4XRAccpStNum[iGap].size()<=4){
      cout<<"Error: Number of accepted strips with satisfication of statistical threshold is equal to or less than 4 for gap "<<iGap+1<<" of c4XR!"<<endl;
      cout<<"Please reset statistical threshold."<<endl;
      return 0;
    }
    chargeTot=0;
    clusTot=0;
    for(UInt_t i=0;i<c4XRAccpStNum[iGap].size();i++){
      chargeTot+=(c4XRaM2[iGap][c4XRAccpStNum[iGap][i]]+c4XRaM1[iGap][c4XRAccpStNum[iGap][i]]+c4XRa0[iGap][c4XRAccpStNum[iGap][i]]+c4XRaP1[iGap][c4XRAccpStNum[iGap][i]]+c4XRaP2[iGap][c4XRAccpStNum[iGap][i]])*c4XRNumClusTot[iGap][c4XRAccpStNum[iGap][i]];
      clusTot+=c4XRNumClusTot[iGap][c4XRAccpStNum[iGap][i]];
    }
    c4XRChargeAverLevel[iGap]=chargeTot/clusTot;
  }
  for(int iGap=0;iGap<6;iGap++){
    c4XRb[iGap][0][0]=c4XRa0[iGap][0];
    c4XRb[iGap][0][1]=c4XRaP1[iGap][0];
    c4XRb[iGap][0][2]=c4XRaP2[iGap][0];
    c4XRc[iGap][0]=c4XRChargeAverLevel[iGap];
    c4XRb[iGap][1][0]=c4XRaM1[iGap][1];
    c4XRb[iGap][1][1]=c4XRa0[iGap][1];
    c4XRb[iGap][1][2]=c4XRaP1[iGap][1];
    c4XRb[iGap][1][3]=c4XRaP2[iGap][1];
    c4XRc[iGap][1]=c4XRChargeAverLevel[iGap];
    for(UInt_t i=2;i<=c4XRNumSt-3;i++){
      c4XRb[iGap][i][i-2]=c4XRaM2[iGap][i];
      c4XRb[iGap][i][i-1]=c4XRaM1[iGap][i];
      c4XRb[iGap][i][i]=c4XRa0[iGap][i];
      c4XRb[iGap][i][i+1]=c4XRaP1[iGap][i];
      c4XRb[iGap][i][i+2]=c4XRaP2[iGap][i];
      c4XRc[iGap][i]=c4XRChargeAverLevel[iGap];
    }
    c4XRb[iGap][c4XRNumSt-2][c4XRNumSt-4]=c4XRaM2[iGap][c4XRNumSt-2];
    c4XRb[iGap][c4XRNumSt-2][c4XRNumSt-3]=c4XRaM1[iGap][c4XRNumSt-2];
    c4XRb[iGap][c4XRNumSt-2][c4XRNumSt-2]=c4XRa0[iGap][c4XRNumSt-2];
    c4XRb[iGap][c4XRNumSt-2][c4XRNumSt-1]=c4XRaP1[iGap][c4XRNumSt-2];
    c4XRc[iGap][c4XRNumSt-2]=c4XRChargeAverLevel[iGap];
    c4XRb[iGap][c4XRNumSt-1][c4XRNumSt-3]=c4XRaM2[iGap][c4XRNumSt-1];
    c4XRb[iGap][c4XRNumSt-1][c4XRNumSt-2]=c4XRaM1[iGap][c4XRNumSt-1];
    c4XRb[iGap][c4XRNumSt-1][c4XRNumSt-1]=c4XRa0[iGap][c4XRNumSt-1];
    c4XRc[iGap][c4XRNumSt-1]=c4XRChargeAverLevel[iGap];
  }

  vector<int> c4XRRejeStNumCable[6][5];
  for(int iGap=0;iGap<6;iGap++){
    for(UInt_t i=0;i<c4XRRejeStNum[iGap].size();i++){
      for(int j=0;j<5;j++){
        if(j<4){
          if(c4XRRejeStNum[iGap][i]>=16*j && c4XRRejeStNum[iGap][i]<16*(j+1)) c4XRRejeStNumCable[iGap][j].push_back(c4XRRejeStNum[iGap][i]);
        }
        else{
          if(c4XRRejeStNum[iGap][i]>=64 && c4XRRejeStNum[iGap][i]<72) c4XRRejeStNumCable[iGap][j].push_back(c4XRRejeStNum[iGap][i]);
        }
      }
    }
  }
  for(int iGap=0;iGap<6;iGap++){
    TMatrix c4XRB(c4XRAccpStNum[iGap].size(),c4XRAccpStNum[iGap].size()), c4XRW(c4XRAccpStNum[iGap].size(),c4XRAccpStNum[iGap].size());
    for(UInt_t i=0;i<c4XRAccpStNum[iGap].size();i++){
      for(UInt_t j=0;j<c4XRAccpStNum[iGap].size();j++){
        c4XRB(i,j)=c4XRb[iGap][c4XRAccpStNum[iGap][i]][c4XRAccpStNum[iGap][j]];
        for(int k=0;k<5;k++){
          if(k<4){
            for(UInt_t l=0;l<c4XRRejeStNumCable[iGap][k].size();l++){
              if(c4XRAccpStNum[iGap][j]>=16*k && c4XRAccpStNum[iGap][j]<16*(k+1)) c4XRB(i,j)+=c4XRb[iGap][c4XRAccpStNum[iGap][i]][c4XRRejeStNumCable[iGap][k][l]]/(16-c4XRRejeStNumCable[iGap][k].size());
            }
          }
          else{
            for(UInt_t l=0;l<c4XRRejeStNumCable[iGap][k].size();l++){
              if(c4XRAccpStNum[iGap][j]>=64 && c4XRAccpStNum[iGap][j]<72) c4XRB(i,j)+=c4XRb[iGap][c4XRAccpStNum[iGap][i]][c4XRRejeStNumCable[iGap][k][l]]/(8-c4XRRejeStNumCable[iGap][k].size());
            }
          }
        }
      }
    }
    for(UInt_t j=0;j<c4XRAccpStNum[iGap].size();j++){
      c4XRW=c4XRB;
      for(UInt_t i=0;i<c4XRAccpStNum[iGap].size();i++){
        c4XRW(i,j)=c4XRc[iGap][i];
      }
      c4XRGainWeightAccp[iGap][c4XRAccpStNum[iGap][j]]=c4XRW.Determinant()/c4XRB.Determinant();
    }
  }
  double c4XRGainWeightAccpSum[6][5]={{0}},c4XRGainWeightAccpAver[6][5];
  int c4XRAccpNumStCable[6][5]={{0}};
  for(int iGap=0;iGap<6;iGap++){
    for(UInt_t i=0;i<c4XRAccpStNum[iGap].size();i++){
      for(int j=0;j<5;j++){
        if(j<4){
          if(c4XRAccpStNum[iGap][i]>=16*j && c4XRAccpStNum[iGap][i]<16*(j+1)){
            c4XRGainWeightAccpSum[iGap][j]+=c4XRGainWeightAccp[iGap][c4XRAccpStNum[iGap][i]];
            c4XRAccpNumStCable[iGap][j]++;
            break;
          }
        }
        else{
          if(c4XRAccpStNum[iGap][i]>=64 && c4XRAccpStNum[iGap][i]<72){
            c4XRGainWeightAccpSum[iGap][j]+=c4XRGainWeightAccp[iGap][c4XRAccpStNum[iGap][i]];
            c4XRAccpNumStCable[iGap][j]++;
            break;
          }
        }
      }
    }
    for(int j=0;j<5;j++){
      if(c4XRAccpNumStCable[iGap][j]!=0) c4XRGainWeightAccpAver[iGap][j]=c4XRGainWeightAccpSum[iGap][j]/c4XRAccpNumStCable[iGap][j];
      else c4XRGainWeightAccpAver[iGap][j]=1;
    }
  }
  for(int iGap=0;iGap<6;iGap++){
    if(c4XRAccpStNum[iGap][0]!=0){
      for(int j=0;j<c4XRAccpStNum[iGap][0];j++){
        for(int k=0;k<5;k++){
          if(k<4){
            if(j>=16*k && j<16*(k+1)){
              c4XRGainWeight[iGap][j]=c4XRGainWeightAccpAver[iGap][k];
              break;
            }
          }
          else{
            if(j>=64 && j<72){
              c4XRGainWeight[iGap][j]=c4XRGainWeightAccpAver[iGap][k];
              break;
            }
          }
        }
      }
    }
    for(UInt_t i=0;i<c4XRAccpStNum[iGap].size()-1;i++){
      c4XRGainWeight[iGap][c4XRAccpStNum[iGap][i]]=c4XRGainWeightAccp[iGap][c4XRAccpStNum[iGap][i]];
      c4XRGainWeight[iGap][c4XRAccpStNum[iGap][i+1]]=c4XRGainWeightAccp[iGap][c4XRAccpStNum[iGap][i+1]];
      for(int j=c4XRAccpStNum[iGap][i]+1;j<c4XRAccpStNum[iGap][i+1];j++){
        for(int k=0;k<5;k++){
          if(k<4){
            if(j>=16*k && j<16*(k+1)){
              c4XRGainWeight[iGap][j]=c4XRGainWeightAccpAver[iGap][k];
              break;
            }
          }
          else{
            if(j>=64 && j<72){
              c4XRGainWeight[iGap][j]=c4XRGainWeightAccpAver[iGap][k];
              break;
            }
          }
        }
      }
    }
    if(c4XRAccpStNum[iGap][c4XRAccpStNum[iGap].size()-1]!=c4XRNumSt-1){
      for(int j=c4XRAccpStNum[iGap][c4XRAccpStNum[iGap].size()-1]+1;j<c4XRNumSt;j++){
        for(int k=0;k<5;k++){
          if(k<4){
            if(j>=16*k && j<16*(k+1)){
              c4XRGainWeight[iGap][j]=c4XRGainWeightAccpAver[iGap][k];
              break;
            }
          }
          else{
            if(j>=64 && j<72){
              c4XRGainWeight[iGap][j]=c4XRGainWeightAccpAver[iGap][k];
              break;
            }
          }
        }
      }
    }
  }
  ofstream outPutNumClusterC4XR("number-of-cluster-C4XR");
  for(int iGap=0;iGap<6;iGap++){
    for(int iSt=0;iSt<c4XRNumSt;iSt++){
      outPutNumClusterC4XR<<c4XRNumClusTot[iGap][iSt]<<"  ";
    }
    outPutNumClusterC4XR<<endl;
  }
  outPutNumClusterC4XR.close();

  ofstream outPutGainC4XR("gain-weight-C4XR");
  for(int iGap=0;iGap<6;iGap++){
    for(int i=0;i<c4XRNumSt;i++){
      outPutGainC4XR<<c4XRGainWeight[iGap][i]<<"  ";
    }
    outPutGainC4XR<<endl;
  }
  outPutGainC4XR.close();

  ofstream outPutNormXYC4R("normalized-factor-C4R");
  for(int iGap=0;iGap<6;iGap++){
    outPutNormXYC4R<<c4XRChargeAverLevel[iGap]/c4YRChargeAverLevel[iGap]<<"  ";
  }
  outPutNormXYC4R.close();

  for(int iGap=0;iGap<6;iGap++){
    for(int iSt=0;iSt<c4YLNumSt;iSt++){
      c4YLNumClusTot[iGap][iSt]=c4YLNumClusSing[iGap][iSt]+c4YLNumClusDoub[iGap][iSt]+c4YLNumClusTrip[iGap][iSt];
      if(c4YLNumClusTot[iGap][iSt]>statThreshold){
        c4YLaM2[iGap][iSt]=c4YLSumChargeTripM2[iGap][iSt]/c4YLNumClusTot[iGap][iSt];
        c4YLaM1[iGap][iSt]=(c4YLSumChargeDoubM1[iGap][iSt]+c4YLSumChargeTripM1[iGap][iSt])/c4YLNumClusTot[iGap][iSt];
        c4YLa0[iGap][iSt]=(c4YLSumChargeSing[iGap][iSt]+c4YLSumChargeDoub0[iGap][iSt]+c4YLSumChargeTrip0[iGap][iSt])/c4YLNumClusTot[iGap][iSt];
        c4YLaP1[iGap][iSt]=(c4YLSumChargeDoubP1[iGap][iSt]+c4YLSumChargeTripP1[iGap][iSt])/c4YLNumClusTot[iGap][iSt];
        c4YLaP2[iGap][iSt]=c4YLSumChargeTripP2[iGap][iSt]/c4YLNumClusTot[iGap][iSt];
        c4YLAccpStNum[iGap].push_back(iSt);
      }
      else c4YLRejeStNum[iGap].push_back(iSt);
    }
  }
  for(int iGap=0;iGap<6;iGap++){
    if(c4YLAccpStNum[iGap].size()<=4){
      cout<<"Error: Number of accepted strips with satisfication of statistical threshold is equal to or less than 4 for gap "<<iGap+1<<" of c4YL!"<<endl;
      cout<<"Please reset statistical threshold."<<endl;
      return 0;
    }
    chargeTot=0;
    clusTot=0;
    for(UInt_t i=0;i<c4YLAccpStNum[iGap].size();i++){
      chargeTot+=(c4YLaM2[iGap][c4YLAccpStNum[iGap][i]]+c4YLaM1[iGap][c4YLAccpStNum[iGap][i]]+c4YLa0[iGap][c4YLAccpStNum[iGap][i]]+c4YLaP1[iGap][c4YLAccpStNum[iGap][i]]+c4YLaP2[iGap][c4YLAccpStNum[iGap][i]])*c4YLNumClusTot[iGap][c4YLAccpStNum[iGap][i]];
      clusTot+=c4YLNumClusTot[iGap][c4YLAccpStNum[iGap][i]];
    }
    c4YLChargeAverLevel[iGap]=chargeTot/clusTot;
  }
  for(int iGap=0;iGap<6;iGap++){
    c4YLb[iGap][0][0]=c4YLa0[iGap][0];
    c4YLb[iGap][0][1]=c4YLaP1[iGap][0];
    c4YLb[iGap][0][2]=c4YLaP2[iGap][0];
    c4YLc[iGap][0]=c4YLChargeAverLevel[iGap];
    c4YLb[iGap][1][0]=c4YLaM1[iGap][1];
    c4YLb[iGap][1][1]=c4YLa0[iGap][1];
    c4YLb[iGap][1][2]=c4YLaP1[iGap][1];
    c4YLb[iGap][1][3]=c4YLaP2[iGap][1];
    c4YLc[iGap][1]=c4YLChargeAverLevel[iGap];
    for(UInt_t i=2;i<=c4YLNumSt-3;i++){
      c4YLb[iGap][i][i-2]=c4YLaM2[iGap][i];
      c4YLb[iGap][i][i-1]=c4YLaM1[iGap][i];
      c4YLb[iGap][i][i]=c4YLa0[iGap][i];
      c4YLb[iGap][i][i+1]=c4YLaP1[iGap][i];
      c4YLb[iGap][i][i+2]=c4YLaP2[iGap][i];
      c4YLc[iGap][i]=c4YLChargeAverLevel[iGap];
    }
    c4YLb[iGap][c4YLNumSt-2][c4YLNumSt-4]=c4YLaM2[iGap][c4YLNumSt-2];
    c4YLb[iGap][c4YLNumSt-2][c4YLNumSt-3]=c4YLaM1[iGap][c4YLNumSt-2];
    c4YLb[iGap][c4YLNumSt-2][c4YLNumSt-2]=c4YLa0[iGap][c4YLNumSt-2];
    c4YLb[iGap][c4YLNumSt-2][c4YLNumSt-1]=c4YLaP1[iGap][c4YLNumSt-2];
    c4YLc[iGap][c4YLNumSt-2]=c4YLChargeAverLevel[iGap];
    c4YLb[iGap][c4YLNumSt-1][c4YLNumSt-3]=c4YLaM2[iGap][c4YLNumSt-1];
    c4YLb[iGap][c4YLNumSt-1][c4YLNumSt-2]=c4YLaM1[iGap][c4YLNumSt-1];
    c4YLb[iGap][c4YLNumSt-1][c4YLNumSt-1]=c4YLa0[iGap][c4YLNumSt-1];
    c4YLc[iGap][c4YLNumSt-1]=c4YLChargeAverLevel[iGap];
  }
  for(int iGap=0;iGap<6;iGap++){
    TMatrix c4YLB(c4YLAccpStNum[iGap].size(),c4YLAccpStNum[iGap].size()), c4YLW(c4YLAccpStNum[iGap].size(),c4YLAccpStNum[iGap].size());
    for(UInt_t i=0;i<c4YLAccpStNum[iGap].size();i++){
      for(UInt_t j=0;j<c4YLAccpStNum[iGap].size();j++){
        c4YLB(i,j)=c4YLb[iGap][c4YLAccpStNum[iGap][i]][c4YLAccpStNum[iGap][j]];
        for(UInt_t k=0;k<c4YLRejeStNum[iGap].size();k++){
          c4YLB(i,j)+=c4YLb[iGap][c4YLAccpStNum[iGap][i]][c4YLRejeStNum[iGap][k]]/(c4YLNumSt-c4YLRejeStNum[iGap].size());
        }
      }
    }

    for(UInt_t j=0;j<c4YLAccpStNum[iGap].size();j++){
      c4YLW=c4YLB;
      for(UInt_t i=0;i<c4YLAccpStNum[iGap].size();i++){
        c4YLW(i,j)=c4YLc[iGap][i];
      }
      c4YLGainWeightAccp[iGap][c4YLAccpStNum[iGap][j]]=c4YLW.Determinant()/c4YLB.Determinant();
    }
  }
  double c4YLGainWeightAccpSum[6]={0},c4YLGainWeightAccpAver[6];
  for(int iGap=0;iGap<6;iGap++){
    for(UInt_t i=0;i<c4YLAccpStNum[iGap].size();i++){
      c4YLGainWeightAccpSum[iGap]+=c4YLGainWeightAccp[iGap][c4YLAccpStNum[iGap][i]];
    }
    if(c4YLAccpStNum[iGap].size()!=0) c4YLGainWeightAccpAver[iGap]=c4YLGainWeightAccpSum[iGap]/c4YLAccpStNum[iGap].size();
    else c4YLGainWeightAccpAver[iGap]=0;
  }
  for(int iGap=0;iGap<6;iGap++){
    if(c4YLAccpStNum[iGap][0]!=0){
      for(int j=0;j<c4YLAccpStNum[iGap][0];j++){
        c4YLGainWeight[iGap][j]=c4YLGainWeightAccpAver[iGap];
      }
    }
    for(UInt_t i=0;i<c4YLAccpStNum[iGap].size()-1;i++){
      c4YLGainWeight[iGap][c4YLAccpStNum[iGap][i]]=c4YLGainWeightAccp[iGap][c4YLAccpStNum[iGap][i]];
      c4YLGainWeight[iGap][c4YLAccpStNum[iGap][i+1]]=c4YLGainWeightAccp[iGap][c4YLAccpStNum[iGap][i+1]];
      for(int j=c4YLAccpStNum[iGap][i]+1;j<c4YLAccpStNum[iGap][i+1];j++){
        c4YLGainWeight[iGap][j]=c4YLGainWeightAccpAver[iGap];
      }
    }
    if(c4YLAccpStNum[iGap][c4YLAccpStNum[iGap].size()-1]!=c4YLNumSt-1){
      for(int j=c4YLAccpStNum[iGap][c4YLAccpStNum[iGap].size()-1]+1;j<c4YLNumSt;j++){
        c4YLGainWeight[iGap][j]=c4YLGainWeightAccpAver[iGap];
      }
    }
  }
  ofstream outPutNumClusterC4YL("number-of-cluster-C4YL");
  for(int iGap=0;iGap<6;iGap++){
    for(int iSt=0;iSt<c4YLNumSt;iSt++){
      outPutNumClusterC4YL<<c4YLNumClusTot[iGap][iSt]<<"  ";
    }
    outPutNumClusterC4YL<<endl;
  }
  outPutNumClusterC4YL.close();

  ofstream outPutGainC4YL("gain-weight-C4YL");
  for(int iGap=0;iGap<6;iGap++){
    for(int i=0;i<c4YLNumSt;i++){
      outPutGainC4YL<<c4YLGainWeight[iGap][i]<<"  ";
    }
    outPutGainC4YL<<endl;
  }
  outPutGainC4YL.close();

  for(int iGap=0;iGap<6;iGap++){
    for(int iSt=0;iSt<c4XLNumSt;iSt++){
      c4XLNumClusTot[iGap][iSt]=c4XLNumClusSing[iGap][iSt]+c4XLNumClusDoub[iGap][iSt]+c4XLNumClusTrip[iGap][iSt];
      if(c4XLNumClusTot[iGap][iSt]>statThreshold){
        c4XLaM2[iGap][iSt]=c4XLSumChargeTripM2[iGap][iSt]/c4XLNumClusTot[iGap][iSt];
        c4XLaM1[iGap][iSt]=(c4XLSumChargeDoubM1[iGap][iSt]+c4XLSumChargeTripM1[iGap][iSt])/c4XLNumClusTot[iGap][iSt];
        c4XLa0[iGap][iSt]=(c4XLSumChargeSing[iGap][iSt]+c4XLSumChargeDoub0[iGap][iSt]+c4XLSumChargeTrip0[iGap][iSt])/c4XLNumClusTot[iGap][iSt];
        c4XLaP1[iGap][iSt]=(c4XLSumChargeDoubP1[iGap][iSt]+c4XLSumChargeTripP1[iGap][iSt])/c4XLNumClusTot[iGap][iSt];
        c4XLaP2[iGap][iSt]=c4XLSumChargeTripP2[iGap][iSt]/c4XLNumClusTot[iGap][iSt];
        c4XLAccpStNum[iGap].push_back(iSt);
      }
      else c4XLRejeStNum[iGap].push_back(iSt);
    }
  }
  for(int iGap=0;iGap<6;iGap++){
    if(c4XLAccpStNum[iGap].size()<=4){
      cout<<"Error: Number of accepted strips with satisfication of statistical threshold is equal to or less than 4 for gap "<<iGap+1<<" of c4XL!"<<endl;
      cout<<"Please reset statistical threshold."<<endl;
      return 0;
    }
    chargeTot=0;
    clusTot=0;
    for(UInt_t i=0;i<c4XLAccpStNum[iGap].size();i++){
      chargeTot+=(c4XLaM2[iGap][c4XLAccpStNum[iGap][i]]+c4XLaM1[iGap][c4XLAccpStNum[iGap][i]]+c4XLa0[iGap][c4XLAccpStNum[iGap][i]]+c4XLaP1[iGap][c4XLAccpStNum[iGap][i]]+c4XLaP2[iGap][c4XLAccpStNum[iGap][i]])*c4XLNumClusTot[iGap][c4XLAccpStNum[iGap][i]];
      clusTot+=c4XLNumClusTot[iGap][c4XLAccpStNum[iGap][i]];
    }
    c4XLChargeAverLevel[iGap]=chargeTot/clusTot;
  }
  for(int iGap=0;iGap<6;iGap++){
    c4XLb[iGap][0][0]=c4XLa0[iGap][0];
    c4XLb[iGap][0][1]=c4XLaP1[iGap][0];
    c4XLb[iGap][0][2]=c4XLaP2[iGap][0];
    c4XLc[iGap][0]=c4XLChargeAverLevel[iGap];
    c4XLb[iGap][1][0]=c4XLaM1[iGap][1];
    c4XLb[iGap][1][1]=c4XLa0[iGap][1];
    c4XLb[iGap][1][2]=c4XLaP1[iGap][1];
    c4XLb[iGap][1][3]=c4XLaP2[iGap][1];
    c4XLc[iGap][1]=c4XLChargeAverLevel[iGap];
    for(UInt_t i=2;i<=c4XLNumSt-3;i++){
      c4XLb[iGap][i][i-2]=c4XLaM2[iGap][i];
      c4XLb[iGap][i][i-1]=c4XLaM1[iGap][i];
      c4XLb[iGap][i][i]=c4XLa0[iGap][i];
      c4XLb[iGap][i][i+1]=c4XLaP1[iGap][i];
      c4XLb[iGap][i][i+2]=c4XLaP2[iGap][i];
      c4XLc[iGap][i]=c4XLChargeAverLevel[iGap];
    }
    c4XLb[iGap][c4XLNumSt-2][c4XLNumSt-4]=c4XLaM2[iGap][c4XLNumSt-2];
    c4XLb[iGap][c4XLNumSt-2][c4XLNumSt-3]=c4XLaM1[iGap][c4XLNumSt-2];
    c4XLb[iGap][c4XLNumSt-2][c4XLNumSt-2]=c4XLa0[iGap][c4XLNumSt-2];
    c4XLb[iGap][c4XLNumSt-2][c4XLNumSt-1]=c4XLaP1[iGap][c4XLNumSt-2];
    c4XLc[iGap][c4XLNumSt-2]=c4XLChargeAverLevel[iGap];
    c4XLb[iGap][c4XLNumSt-1][c4XLNumSt-3]=c4XLaM2[iGap][c4XLNumSt-1];
    c4XLb[iGap][c4XLNumSt-1][c4XLNumSt-2]=c4XLaM1[iGap][c4XLNumSt-1];
    c4XLb[iGap][c4XLNumSt-1][c4XLNumSt-1]=c4XLa0[iGap][c4XLNumSt-1];
    c4XLc[iGap][c4XLNumSt-1]=c4XLChargeAverLevel[iGap];
  }

  vector<int> c4XLRejeStNumCable[6][5];
  for(int iGap=0;iGap<6;iGap++){
    for(UInt_t i=0;i<c4XLRejeStNum[iGap].size();i++){
      for(int j=0;j<5;j++){
        if(j<4){
          if(c4XLRejeStNum[iGap][i]>=16*j && c4XLRejeStNum[iGap][i]<16*(j+1)) c4XLRejeStNumCable[iGap][j].push_back(c4XLRejeStNum[iGap][i]);
        }
        else{
          if(c4XLRejeStNum[iGap][i]>=64 && c4XLRejeStNum[iGap][i]<72) c4XLRejeStNumCable[iGap][j].push_back(c4XLRejeStNum[iGap][i]);
        }
      }
    }
  }
  for(int iGap=0;iGap<6;iGap++){
    TMatrix c4XLB(c4XLAccpStNum[iGap].size(),c4XLAccpStNum[iGap].size()), c4XLW(c4XLAccpStNum[iGap].size(),c4XLAccpStNum[iGap].size());
    for(UInt_t i=0;i<c4XLAccpStNum[iGap].size();i++){
      for(UInt_t j=0;j<c4XLAccpStNum[iGap].size();j++){
        c4XLB(i,j)=c4XLb[iGap][c4XLAccpStNum[iGap][i]][c4XLAccpStNum[iGap][j]];
        for(int k=0;k<5;k++){
          if(k<4){
            for(UInt_t l=0;l<c4XLRejeStNumCable[iGap][k].size();l++){
              if(c4XLAccpStNum[iGap][j]>=16*k && c4XLAccpStNum[iGap][j]<16*(k+1)) c4XLB(i,j)+=c4XLb[iGap][c4XLAccpStNum[iGap][i]][c4XLRejeStNumCable[iGap][k][l]]/(16-c4XLRejeStNumCable[iGap][k].size());
            }
          }
          else{
            for(UInt_t l=0;l<c4XLRejeStNumCable[iGap][k].size();l++){
              if(c4XLAccpStNum[iGap][j]>=64 && c4XLAccpStNum[iGap][j]<72) c4XLB(i,j)+=c4XLb[iGap][c4XLAccpStNum[iGap][i]][c4XLRejeStNumCable[iGap][k][l]]/(8-c4XLRejeStNumCable[iGap][k].size());
            }
          }
        }
      }
    }
    for(UInt_t j=0;j<c4XLAccpStNum[iGap].size();j++){
      c4XLW=c4XLB;
      for(UInt_t i=0;i<c4XLAccpStNum[iGap].size();i++){
        c4XLW(i,j)=c4XLc[iGap][i];
      }
      c4XLGainWeightAccp[iGap][c4XLAccpStNum[iGap][j]]=c4XLW.Determinant()/c4XLB.Determinant();
    }
  }
  double c4XLGainWeightAccpSum[6][5]={{0}},c4XLGainWeightAccpAver[6][5];
  int c4XLAccpNumStCable[6][5]={{0}};
  for(int iGap=0;iGap<6;iGap++){
    for(UInt_t i=0;i<c4XLAccpStNum[iGap].size();i++){
      for(int j=0;j<5;j++){
        if(j<4){
          if(c4XLAccpStNum[iGap][i]>=16*j && c4XLAccpStNum[iGap][i]<16*(j+1)){
            c4XLGainWeightAccpSum[iGap][j]+=c4XLGainWeightAccp[iGap][c4XLAccpStNum[iGap][i]];
            c4XLAccpNumStCable[iGap][j]++;
            break;
          }
        }
        else{
          if(c4XLAccpStNum[iGap][i]>=64 && c4XLAccpStNum[iGap][i]<72){
            c4XLGainWeightAccpSum[iGap][j]+=c4XLGainWeightAccp[iGap][c4XLAccpStNum[iGap][i]];
            c4XLAccpNumStCable[iGap][j]++;
            break;
          }
        }
      }
    }
    for(int j=0;j<5;j++){
      if(c4XLAccpNumStCable[iGap][j]!=0) c4XLGainWeightAccpAver[iGap][j]=c4XLGainWeightAccpSum[iGap][j]/c4XLAccpNumStCable[iGap][j];
      else c4XLGainWeightAccpAver[iGap][j]=1;
    }
  }
  for(int iGap=0;iGap<6;iGap++){
    if(c4XLAccpStNum[iGap][0]!=0){
      for(int j=0;j<c4XLAccpStNum[iGap][0];j++){
        for(int k=0;k<5;k++){
          if(k<4){
            if(j>=16*k && j<16*(k+1)){
              c4XLGainWeight[iGap][j]=c4XLGainWeightAccpAver[iGap][k];
              break;
            }
          }
          else{
            if(j>=64 && j<72){
              c4XLGainWeight[iGap][j]=c4XLGainWeightAccpAver[iGap][k];
              break;
            }
          }
        }
      }
    }
    for(UInt_t i=0;i<c4XLAccpStNum[iGap].size()-1;i++){
      c4XLGainWeight[iGap][c4XLAccpStNum[iGap][i]]=c4XLGainWeightAccp[iGap][c4XLAccpStNum[iGap][i]];
      c4XLGainWeight[iGap][c4XLAccpStNum[iGap][i+1]]=c4XLGainWeightAccp[iGap][c4XLAccpStNum[iGap][i+1]];
      for(int j=c4XLAccpStNum[iGap][i]+1;j<c4XLAccpStNum[iGap][i+1];j++){
        for(int k=0;k<5;k++){
          if(k<4){
            if(j>=16*k && j<16*(k+1)){
              c4XLGainWeight[iGap][j]=c4XLGainWeightAccpAver[iGap][k];
              break;
            }
          }
          else{
            if(j>=64 && j<72){
              c4XLGainWeight[iGap][j]=c4XLGainWeightAccpAver[iGap][k];
              break;
            }
          }
        }
      }
    }
    if(c4XLAccpStNum[iGap][c4XLAccpStNum[iGap].size()-1]!=c4XLNumSt-1){
      for(int j=c4XLAccpStNum[iGap][c4XLAccpStNum[iGap].size()-1]+1;j<c4XLNumSt;j++){
        for(int k=0;k<5;k++){
          if(k<4){
            if(j>=16*k && j<16*(k+1)){
              c4XLGainWeight[iGap][j]=c4XLGainWeightAccpAver[iGap][k];
              break;
            }
          }
          else{
            if(j>=64 && j<72){
              c4XLGainWeight[iGap][j]=c4XLGainWeightAccpAver[iGap][k];
              break;
            }
          }
        }
      }
    }
  }
  ofstream outPutNumClusterC4XL("number-of-cluster-C4XL");
  for(int iGap=0;iGap<6;iGap++){
    for(int iSt=0;iSt<c4XLNumSt;iSt++){
      outPutNumClusterC4XL<<c4XLNumClusTot[iGap][iSt]<<"  ";
    }
    outPutNumClusterC4XL<<endl;
  }
  outPutNumClusterC4XL.close();

  ofstream outPutGainC4XL("gain-weight-C4XL");
  for(int iGap=0;iGap<6;iGap++){
    for(int i=0;i<c4XLNumSt;i++){
      outPutGainC4XL<<c4XLGainWeight[iGap][i]<<"  ";
    }
    outPutGainC4XL<<endl;
  }
  outPutGainC4XL.close();
  ofstream outPutNormXYC4L("normalized-factor-C4L");
  for(int iGap=0;iGap<6;iGap++){
    outPutNormXYC4L<<c4XLChargeAverLevel[iGap]/c4YLChargeAverLevel[iGap]<<"  ";
  }
  outPutNormXYC4L.close();

  for(int iGap=0;iGap<6;iGap++){
    for(int iSt=0;iSt<c2YRNumSt;iSt++){
      c2YRNumClusTot[iGap][iSt]=c2YRNumClusSing[iGap][iSt]+c2YRNumClusDoub[iGap][iSt]+c2YRNumClusTrip[iGap][iSt];
      if(c2YRNumClusTot[iGap][iSt]>statThreshold){
        c2YRaM2[iGap][iSt]=c2YRSumChargeTripM2[iGap][iSt]/c2YRNumClusTot[iGap][iSt];
        c2YRaM1[iGap][iSt]=(c2YRSumChargeDoubM1[iGap][iSt]+c2YRSumChargeTripM1[iGap][iSt])/c2YRNumClusTot[iGap][iSt];
        c2YRa0[iGap][iSt]=(c2YRSumChargeSing[iGap][iSt]+c2YRSumChargeDoub0[iGap][iSt]+c2YRSumChargeTrip0[iGap][iSt])/c2YRNumClusTot[iGap][iSt];
        c2YRaP1[iGap][iSt]=(c2YRSumChargeDoubP1[iGap][iSt]+c2YRSumChargeTripP1[iGap][iSt])/c2YRNumClusTot[iGap][iSt];
        c2YRaP2[iGap][iSt]=c2YRSumChargeTripP2[iGap][iSt]/c2YRNumClusTot[iGap][iSt];
        c2YRAccpStNum[iGap].push_back(iSt);
      }
      else c2YRRejeStNum[iGap].push_back(iSt);
    }
  }
  for(int iGap=0;iGap<6;iGap++){
    if(c2YRAccpStNum[iGap].size()<=4){
      cout<<"Error: Number of accepted strips with satisfication of statistical threshold is equal to or less than 4 for gap "<<iGap+1<<" of c2YR!"<<endl;
      cout<<"Please reset statistical threshold."<<endl;
      return 0;
    }
    chargeTot=0;
    clusTot=0;
    for(UInt_t i=0;i<c2YRAccpStNum[iGap].size();i++){
      chargeTot+=(c2YRaM2[iGap][c2YRAccpStNum[iGap][i]]+c2YRaM1[iGap][c2YRAccpStNum[iGap][i]]+c2YRa0[iGap][c2YRAccpStNum[iGap][i]]+c2YRaP1[iGap][c2YRAccpStNum[iGap][i]]+c2YRaP2[iGap][c2YRAccpStNum[iGap][i]])*c2YRNumClusTot[iGap][c2YRAccpStNum[iGap][i]];
      clusTot+=c2YRNumClusTot[iGap][c2YRAccpStNum[iGap][i]];
    }
    c2YRChargeAverLevel[iGap]=chargeTot/clusTot;
  }
  for(int iGap=0;iGap<6;iGap++){
    c2YRb[iGap][0][0]=c2YRa0[iGap][0];
    c2YRb[iGap][0][1]=c2YRaP1[iGap][0];
    c2YRb[iGap][0][2]=c2YRaP2[iGap][0];
    c2YRc[iGap][0]=c2YRChargeAverLevel[iGap];
    c2YRb[iGap][1][0]=c2YRaM1[iGap][1];
    c2YRb[iGap][1][1]=c2YRa0[iGap][1];
    c2YRb[iGap][1][2]=c2YRaP1[iGap][1];
    c2YRb[iGap][1][3]=c2YRaP2[iGap][1];
    c2YRc[iGap][1]=c2YRChargeAverLevel[iGap];
    for(UInt_t i=2;i<=c2YRNumSt-3;i++){
      c2YRb[iGap][i][i-2]=c2YRaM2[iGap][i];
      c2YRb[iGap][i][i-1]=c2YRaM1[iGap][i];
      c2YRb[iGap][i][i]=c2YRa0[iGap][i];
      c2YRb[iGap][i][i+1]=c2YRaP1[iGap][i];
      c2YRb[iGap][i][i+2]=c2YRaP2[iGap][i];
      c2YRc[iGap][i]=c2YRChargeAverLevel[iGap];
    }
    c2YRb[iGap][c2YRNumSt-2][c2YRNumSt-4]=c2YRaM2[iGap][c2YRNumSt-2];
    c2YRb[iGap][c2YRNumSt-2][c2YRNumSt-3]=c2YRaM1[iGap][c2YRNumSt-2];
    c2YRb[iGap][c2YRNumSt-2][c2YRNumSt-2]=c2YRa0[iGap][c2YRNumSt-2];
    c2YRb[iGap][c2YRNumSt-2][c2YRNumSt-1]=c2YRaP1[iGap][c2YRNumSt-2];
    c2YRc[iGap][c2YRNumSt-2]=c2YRChargeAverLevel[iGap];
    c2YRb[iGap][c2YRNumSt-1][c2YRNumSt-3]=c2YRaM2[iGap][c2YRNumSt-1];
    c2YRb[iGap][c2YRNumSt-1][c2YRNumSt-2]=c2YRaM1[iGap][c2YRNumSt-1];
    c2YRb[iGap][c2YRNumSt-1][c2YRNumSt-1]=c2YRa0[iGap][c2YRNumSt-1];
    c2YRc[iGap][c2YRNumSt-1]=c2YRChargeAverLevel[iGap];
  }
  for(int iGap=0;iGap<6;iGap++){
    TMatrix c2YRB(c2YRAccpStNum[iGap].size(),c2YRAccpStNum[iGap].size()), c2YRW(c2YRAccpStNum[iGap].size(),c2YRAccpStNum[iGap].size());
    for(UInt_t i=0;i<c2YRAccpStNum[iGap].size();i++){
      for(UInt_t j=0;j<c2YRAccpStNum[iGap].size();j++){
        c2YRB(i,j)=c2YRb[iGap][c2YRAccpStNum[iGap][i]][c2YRAccpStNum[iGap][j]];
        for(UInt_t k=0;k<c2YRRejeStNum[iGap].size();k++){
          c2YRB(i,j)+=c2YRb[iGap][c2YRAccpStNum[iGap][i]][c2YRRejeStNum[iGap][k]]/(c2YRNumSt-c2YRRejeStNum[iGap].size());
        }
      }
    }

    for(UInt_t j=0;j<c2YRAccpStNum[iGap].size();j++){
      c2YRW=c2YRB;
      for(UInt_t i=0;i<c2YRAccpStNum[iGap].size();i++){
        c2YRW(i,j)=c2YRc[iGap][i];
      }
      c2YRGainWeightAccp[iGap][c2YRAccpStNum[iGap][j]]=c2YRW.Determinant()/c2YRB.Determinant();
    }
  }
  double c2YRGainWeightAccpSum[6]={0},c2YRGainWeightAccpAver[6];
  for(int iGap=0;iGap<6;iGap++){
    for(UInt_t i=0;i<c2YRAccpStNum[iGap].size();i++){
      c2YRGainWeightAccpSum[iGap]+=c2YRGainWeightAccp[iGap][c2YRAccpStNum[iGap][i]];
    }
    if(c2YRAccpStNum[iGap].size()!=0) c2YRGainWeightAccpAver[iGap]=c2YRGainWeightAccpSum[iGap]/c2YRAccpStNum[iGap].size();
    else c2YRGainWeightAccpAver[iGap]=0;
  }
  for(int iGap=0;iGap<6;iGap++){
    if(c2YRAccpStNum[iGap][0]!=0){
      for(int j=0;j<c2YRAccpStNum[iGap][0];j++){
        c2YRGainWeight[iGap][j]=c2YRGainWeightAccpAver[iGap];
      }
    }
    for(UInt_t i=0;i<c2YRAccpStNum[iGap].size()-1;i++){
      c2YRGainWeight[iGap][c2YRAccpStNum[iGap][i]]=c2YRGainWeightAccp[iGap][c2YRAccpStNum[iGap][i]];
      c2YRGainWeight[iGap][c2YRAccpStNum[iGap][i+1]]=c2YRGainWeightAccp[iGap][c2YRAccpStNum[iGap][i+1]];
      for(int j=c2YRAccpStNum[iGap][i]+1;j<c2YRAccpStNum[iGap][i+1];j++){
        c2YRGainWeight[iGap][j]=c2YRGainWeightAccpAver[iGap];
      }
    }
    if(c2YRAccpStNum[iGap][c2YRAccpStNum[iGap].size()-1]!=c2YRNumSt-1){
      for(int j=c2YRAccpStNum[iGap][c2YRAccpStNum[iGap].size()-1]+1;j<c2YRNumSt;j++){
        c2YRGainWeight[iGap][j]=c2YRGainWeightAccpAver[iGap];
      }
    }
  }
  ofstream outPutNumClusterC2YR("number-of-cluster-C2YR");
  for(int iGap=0;iGap<6;iGap++){
    for(int iSt=0;iSt<c2YRNumSt;iSt++){
      outPutNumClusterC2YR<<c2YRNumClusTot[iGap][iSt]<<"  ";
    }
    outPutNumClusterC2YR<<endl;
  }
  outPutNumClusterC2YR.close();

  ofstream outPutGainC2YR("gain-weight-C2YR");
  for(int iGap=0;iGap<6;iGap++){
    for(int i=0;i<c2YRNumSt;i++){
      outPutGainC2YR<<c2YRGainWeight[iGap][i]<<"  ";
    }
    outPutGainC2YR<<endl;
  }
  outPutGainC2YR.close();

  for(int iGap=0;iGap<6;iGap++){
    for(int iSt=0;iSt<c2XRNumSt;iSt++){
      c2XRNumClusTot[iGap][iSt]=c2XRNumClusSing[iGap][iSt]+c2XRNumClusDoub[iGap][iSt]+c2XRNumClusTrip[iGap][iSt];
      if(c2XRNumClusTot[iGap][iSt]>statThreshold){
        c2XRaM2[iGap][iSt]=c2XRSumChargeTripM2[iGap][iSt]/c2XRNumClusTot[iGap][iSt];
        c2XRaM1[iGap][iSt]=(c2XRSumChargeDoubM1[iGap][iSt]+c2XRSumChargeTripM1[iGap][iSt])/c2XRNumClusTot[iGap][iSt];
        c2XRa0[iGap][iSt]=(c2XRSumChargeSing[iGap][iSt]+c2XRSumChargeDoub0[iGap][iSt]+c2XRSumChargeTrip0[iGap][iSt])/c2XRNumClusTot[iGap][iSt];
        c2XRaP1[iGap][iSt]=(c2XRSumChargeDoubP1[iGap][iSt]+c2XRSumChargeTripP1[iGap][iSt])/c2XRNumClusTot[iGap][iSt];
        c2XRaP2[iGap][iSt]=c2XRSumChargeTripP2[iGap][iSt]/c2XRNumClusTot[iGap][iSt];
        c2XRAccpStNum[iGap].push_back(iSt);
      }
      else c2XRRejeStNum[iGap].push_back(iSt);
    }
  }
  for(int iGap=0;iGap<6;iGap++){
    if(c2XRAccpStNum[iGap].size()<=4){
      cout<<"Error: Number of accepted strips with satisfication of statistical threshold is equal to or less than 4 for gap "<<iGap+1<<" of c2XR!"<<endl;
      cout<<"Please reset statistical threshold."<<endl;
      return 0;
    }
    chargeTot=0;
    clusTot=0;
    for(UInt_t i=0;i<c2XRAccpStNum[iGap].size();i++){
      chargeTot+=(c2XRaM2[iGap][c2XRAccpStNum[iGap][i]]+c2XRaM1[iGap][c2XRAccpStNum[iGap][i]]+c2XRa0[iGap][c2XRAccpStNum[iGap][i]]+c2XRaP1[iGap][c2XRAccpStNum[iGap][i]]+c2XRaP2[iGap][c2XRAccpStNum[iGap][i]])*c2XRNumClusTot[iGap][c2XRAccpStNum[iGap][i]];
      clusTot+=c2XRNumClusTot[iGap][c2XRAccpStNum[iGap][i]];
    }
    c2XRChargeAverLevel[iGap]=chargeTot/clusTot;
  }
  for(int iGap=0;iGap<6;iGap++){
    c2XRb[iGap][0][0]=c2XRa0[iGap][0];
    c2XRb[iGap][0][1]=c2XRaP1[iGap][0];
    c2XRb[iGap][0][2]=c2XRaP2[iGap][0];
    c2XRc[iGap][0]=c2XRChargeAverLevel[iGap];
    c2XRb[iGap][1][0]=c2XRaM1[iGap][1];
    c2XRb[iGap][1][1]=c2XRa0[iGap][1];
    c2XRb[iGap][1][2]=c2XRaP1[iGap][1];
    c2XRb[iGap][1][3]=c2XRaP2[iGap][1];
    c2XRc[iGap][1]=c2XRChargeAverLevel[iGap];
    for(UInt_t i=2;i<=c2XRNumSt-3;i++){
      c2XRb[iGap][i][i-2]=c2XRaM2[iGap][i];
      c2XRb[iGap][i][i-1]=c2XRaM1[iGap][i];
      c2XRb[iGap][i][i]=c2XRa0[iGap][i];
      c2XRb[iGap][i][i+1]=c2XRaP1[iGap][i];
      c2XRb[iGap][i][i+2]=c2XRaP2[iGap][i];
      c2XRc[iGap][i]=c2XRChargeAverLevel[iGap];
    }
    c2XRb[iGap][c2XRNumSt-2][c2XRNumSt-4]=c2XRaM2[iGap][c2XRNumSt-2];
    c2XRb[iGap][c2XRNumSt-2][c2XRNumSt-3]=c2XRaM1[iGap][c2XRNumSt-2];
    c2XRb[iGap][c2XRNumSt-2][c2XRNumSt-2]=c2XRa0[iGap][c2XRNumSt-2];
    c2XRb[iGap][c2XRNumSt-2][c2XRNumSt-1]=c2XRaP1[iGap][c2XRNumSt-2];
    c2XRc[iGap][c2XRNumSt-2]=c2XRChargeAverLevel[iGap];
    c2XRb[iGap][c2XRNumSt-1][c2XRNumSt-3]=c2XRaM2[iGap][c2XRNumSt-1];
    c2XRb[iGap][c2XRNumSt-1][c2XRNumSt-2]=c2XRaM1[iGap][c2XRNumSt-1];
    c2XRb[iGap][c2XRNumSt-1][c2XRNumSt-1]=c2XRa0[iGap][c2XRNumSt-1];
    c2XRc[iGap][c2XRNumSt-1]=c2XRChargeAverLevel[iGap];
  }

  vector<int> c2XRRejeStNumCable[6][4];
  for(int iGap=0;iGap<6;iGap++){
    for(UInt_t i=0;i<c2XRRejeStNum[iGap].size();i++){
      for(int j=0;j<4;j++){
        if(C2XStTrans[c2XRRejeStNum[iGap][i]]>=14*j && C2XStTrans[c2XRRejeStNum[iGap][i]]<14*(j+1)) c2XRRejeStNumCable[iGap][j].push_back(c2XRRejeStNum[iGap][i]);
      }
    }
  }
  for(int iGap=0;iGap<6;iGap++){
    TMatrix c2XRB(c2XRAccpStNum[iGap].size(),c2XRAccpStNum[iGap].size()), c2XRW(c2XRAccpStNum[iGap].size(),c2XRAccpStNum[iGap].size());
    for(UInt_t i=0;i<c2XRAccpStNum[iGap].size();i++){
      for(UInt_t j=0;j<c2XRAccpStNum[iGap].size();j++){
        c2XRB(i,j)=c2XRb[iGap][c2XRAccpStNum[iGap][i]][c2XRAccpStNum[iGap][j]];
        for(int k=0;k<4;k++){
          for(UInt_t l=0;l<c2XRRejeStNumCable[iGap][k].size();l++){
            if(C2XStTrans[c2XRAccpStNum[iGap][j]]>=14*k && C2XStTrans[c2XRAccpStNum[iGap][j]]<14*(k+1)) c2XRB(i,j)+=c2XRb[iGap][c2XRAccpStNum[iGap][i]][c2XRRejeStNumCable[iGap][k][l]]/(14-c2XRRejeStNumCable[iGap][k].size());
          }
        }
      }
    }
    for(UInt_t j=0;j<c2XRAccpStNum[iGap].size();j++){
      c2XRW=c2XRB;
      for(UInt_t i=0;i<c2XRAccpStNum[iGap].size();i++){
        c2XRW(i,j)=c2XRc[iGap][i];
      }
      c2XRGainWeightAccp[iGap][c2XRAccpStNum[iGap][j]]=c2XRW.Determinant()/c2XRB.Determinant();
    }
  }
  double c2XRGainWeightAccpSum[6][4]={{0}},c2XRGainWeightAccpAver[6][4];
  int c2XRAccpNumStCable[6][4]={{0}};
  for(int iGap=0;iGap<6;iGap++){
    for(UInt_t i=0;i<c2XRAccpStNum[iGap].size();i++){
      for(int j=0;j<4;j++){
        if(C2XStTrans[c2XRAccpStNum[iGap][i]]>=14*j && C2XStTrans[c2XRAccpStNum[iGap][i]]<14*(j+1)){
          c2XRGainWeightAccpSum[iGap][j]+=c2XRGainWeightAccp[iGap][c2XRAccpStNum[iGap][i]];
          c2XRAccpNumStCable[iGap][j]++;
          break;
        }
      }
    }
    for(int j=0;j<4;j++){
      if(c2XRAccpNumStCable[iGap][j]!=0) c2XRGainWeightAccpAver[iGap][j]=c2XRGainWeightAccpSum[iGap][j]/c2XRAccpNumStCable[iGap][j];
      else c2XRGainWeightAccpAver[iGap][j]=1;
    }
  }
  for(int iGap=0;iGap<6;iGap++){
    if(c2XRAccpStNum[iGap][0]!=0){
      for(int j=0;j<c2XRAccpStNum[iGap][0];j++){
        for(int k=0;k<4;k++){
          if(C2XStTrans[j]>=14*k && C2XStTrans[j]<14*(k+1)){
            c2XRGainWeight[iGap][j]=c2XRGainWeightAccpAver[iGap][k];
            break;
          }
        }
      }
    }
    for(UInt_t i=0;i<c2XRAccpStNum[iGap].size()-1;i++){
      c2XRGainWeight[iGap][c2XRAccpStNum[iGap][i]]=c2XRGainWeightAccp[iGap][c2XRAccpStNum[iGap][i]];
      c2XRGainWeight[iGap][c2XRAccpStNum[iGap][i+1]]=c2XRGainWeightAccp[iGap][c2XRAccpStNum[iGap][i+1]];
      for(int j=c2XRAccpStNum[iGap][i]+1;j<c2XRAccpStNum[iGap][i+1];j++){
        for(int k=0;k<4;k++){
          if(C2XStTrans[j]>=14*k && C2XStTrans[j]<14*(k+1)){
            c2XRGainWeight[iGap][j]=c2XRGainWeightAccpAver[iGap][k];
            break;
          }
        }
      }
    }
    if(c2XRAccpStNum[iGap][c2XRAccpStNum[iGap].size()-1]!=c2XRNumSt-1){
      for(int j=c2XRAccpStNum[iGap][c2XRAccpStNum[iGap].size()-1]+1;j<c2XRNumSt;j++){
        for(int k=0;k<4;k++){
          if(C2XStTrans[j]>=14*k && C2XStTrans[j]<14*(k+1)){
            c2XRGainWeight[iGap][j]=c2XRGainWeightAccpAver[iGap][k];
            break;
          }
        }
      }
    }
  }
  ofstream outPutNumClusterC2XR("number-of-cluster-C2XR");
  for(int iGap=0;iGap<6;iGap++){
    for(int iSt=0;iSt<c2XRNumSt;iSt++){
      outPutNumClusterC2XR<<c2XRNumClusTot[iGap][iSt]<<"  ";
    }
    outPutNumClusterC2XR<<endl;
  }
  outPutNumClusterC2XR.close();

  ofstream outPutGainC2XR("gain-weight-C2XR");
  for(int iGap=0;iGap<6;iGap++){
    for(int i=0;i<c2XRNumSt;i++){
      outPutGainC2XR<<c2XRGainWeight[iGap][i]<<"  ";
    }
    outPutGainC2XR<<endl;
  }
  outPutGainC2XR.close();

  ofstream outPutNormXYC2R("normalized-factor-C2R");
  for(int iGap=0;iGap<6;iGap++){
    outPutNormXYC2R<<c2XRChargeAverLevel[iGap]/c2YRChargeAverLevel[iGap]<<"  ";
  }
  outPutNormXYC2R.close();

  for(int iGap=0;iGap<6;iGap++){
    for(int iSt=0;iSt<c2YLNumSt;iSt++){
      c2YLNumClusTot[iGap][iSt]=c2YLNumClusSing[iGap][iSt]+c2YLNumClusDoub[iGap][iSt]+c2YLNumClusTrip[iGap][iSt];
      if(c2YLNumClusTot[iGap][iSt]>statThreshold){
        c2YLaM2[iGap][iSt]=c2YLSumChargeTripM2[iGap][iSt]/c2YLNumClusTot[iGap][iSt];
        c2YLaM1[iGap][iSt]=(c2YLSumChargeDoubM1[iGap][iSt]+c2YLSumChargeTripM1[iGap][iSt])/c2YLNumClusTot[iGap][iSt];
        c2YLa0[iGap][iSt]=(c2YLSumChargeSing[iGap][iSt]+c2YLSumChargeDoub0[iGap][iSt]+c2YLSumChargeTrip0[iGap][iSt])/c2YLNumClusTot[iGap][iSt];
        c2YLaP1[iGap][iSt]=(c2YLSumChargeDoubP1[iGap][iSt]+c2YLSumChargeTripP1[iGap][iSt])/c2YLNumClusTot[iGap][iSt];
        c2YLaP2[iGap][iSt]=c2YLSumChargeTripP2[iGap][iSt]/c2YLNumClusTot[iGap][iSt];
        c2YLAccpStNum[iGap].push_back(iSt);
      }
      else c2YLRejeStNum[iGap].push_back(iSt);
    }
  }
  for(int iGap=0;iGap<6;iGap++){
    if(c2YLAccpStNum[iGap].size()<=4){
      cout<<"Error: Number of accepted strips with satisfication of statistical threshold is equal to or less than 4 for gap "<<iGap+1<<" of c2YL!"<<endl;
      cout<<"Please reset statistical threshold."<<endl;
      return 0;
    }
    chargeTot=0;
    clusTot=0;
    for(UInt_t i=0;i<c2YLAccpStNum[iGap].size();i++){
      chargeTot+=(c2YLaM2[iGap][c2YLAccpStNum[iGap][i]]+c2YLaM1[iGap][c2YLAccpStNum[iGap][i]]+c2YLa0[iGap][c2YLAccpStNum[iGap][i]]+c2YLaP1[iGap][c2YLAccpStNum[iGap][i]]+c2YLaP2[iGap][c2YLAccpStNum[iGap][i]])*c2YLNumClusTot[iGap][c2YLAccpStNum[iGap][i]];
      clusTot+=c2YLNumClusTot[iGap][c2YLAccpStNum[iGap][i]];
    }
    c2YLChargeAverLevel[iGap]=chargeTot/clusTot;
  }
  for(int iGap=0;iGap<6;iGap++){
    c2YLb[iGap][0][0]=c2YLa0[iGap][0];
    c2YLb[iGap][0][1]=c2YLaP1[iGap][0];
    c2YLb[iGap][0][2]=c2YLaP2[iGap][0];
    c2YLc[iGap][0]=c2YLChargeAverLevel[iGap];
    c2YLb[iGap][1][0]=c2YLaM1[iGap][1];
    c2YLb[iGap][1][1]=c2YLa0[iGap][1];
    c2YLb[iGap][1][2]=c2YLaP1[iGap][1];
    c2YLb[iGap][1][3]=c2YLaP2[iGap][1];
    c2YLc[iGap][1]=c2YLChargeAverLevel[iGap];
    for(UInt_t i=2;i<=c2YLNumSt-3;i++){
      c2YLb[iGap][i][i-2]=c2YLaM2[iGap][i];
      c2YLb[iGap][i][i-1]=c2YLaM1[iGap][i];
      c2YLb[iGap][i][i]=c2YLa0[iGap][i];
      c2YLb[iGap][i][i+1]=c2YLaP1[iGap][i];
      c2YLb[iGap][i][i+2]=c2YLaP2[iGap][i];
      c2YLc[iGap][i]=c2YLChargeAverLevel[iGap];
    }
    c2YLb[iGap][c2YLNumSt-2][c2YLNumSt-4]=c2YLaM2[iGap][c2YLNumSt-2];
    c2YLb[iGap][c2YLNumSt-2][c2YLNumSt-3]=c2YLaM1[iGap][c2YLNumSt-2];
    c2YLb[iGap][c2YLNumSt-2][c2YLNumSt-2]=c2YLa0[iGap][c2YLNumSt-2];
    c2YLb[iGap][c2YLNumSt-2][c2YLNumSt-1]=c2YLaP1[iGap][c2YLNumSt-2];
    c2YLc[iGap][c2YLNumSt-2]=c2YLChargeAverLevel[iGap];
    c2YLb[iGap][c2YLNumSt-1][c2YLNumSt-3]=c2YLaM2[iGap][c2YLNumSt-1];
    c2YLb[iGap][c2YLNumSt-1][c2YLNumSt-2]=c2YLaM1[iGap][c2YLNumSt-1];
    c2YLb[iGap][c2YLNumSt-1][c2YLNumSt-1]=c2YLa0[iGap][c2YLNumSt-1];
    c2YLc[iGap][c2YLNumSt-1]=c2YLChargeAverLevel[iGap];
  }
  for(int iGap=0;iGap<6;iGap++){
    TMatrix c2YLB(c2YLAccpStNum[iGap].size(),c2YLAccpStNum[iGap].size()), c2YLW(c2YLAccpStNum[iGap].size(),c2YLAccpStNum[iGap].size());
    for(UInt_t i=0;i<c2YLAccpStNum[iGap].size();i++){
      for(UInt_t j=0;j<c2YLAccpStNum[iGap].size();j++){
        c2YLB(i,j)=c2YLb[iGap][c2YLAccpStNum[iGap][i]][c2YLAccpStNum[iGap][j]];
        for(UInt_t k=0;k<c2YLRejeStNum[iGap].size();k++){
          c2YLB(i,j)+=c2YLb[iGap][c2YLAccpStNum[iGap][i]][c2YLRejeStNum[iGap][k]]/(c2YLNumSt-c2YLRejeStNum[iGap].size());
        }
      }
    }

    for(UInt_t j=0;j<c2YLAccpStNum[iGap].size();j++){
      c2YLW=c2YLB;
      for(UInt_t i=0;i<c2YLAccpStNum[iGap].size();i++){
        c2YLW(i,j)=c2YLc[iGap][i];
      }
      c2YLGainWeightAccp[iGap][c2YLAccpStNum[iGap][j]]=c2YLW.Determinant()/c2YLB.Determinant();
    }
  }
  double c2YLGainWeightAccpSum[6]={0},c2YLGainWeightAccpAver[6];
  for(int iGap=0;iGap<6;iGap++){
    for(UInt_t i=0;i<c2YLAccpStNum[iGap].size();i++){
      c2YLGainWeightAccpSum[iGap]+=c2YLGainWeightAccp[iGap][c2YLAccpStNum[iGap][i]];
    }
    if(c2YLAccpStNum[iGap].size()!=0) c2YLGainWeightAccpAver[iGap]=c2YLGainWeightAccpSum[iGap]/c2YLAccpStNum[iGap].size();
    else c2YLGainWeightAccpAver[iGap]=0;
  }
  for(int iGap=0;iGap<6;iGap++){
    if(c2YLAccpStNum[iGap][0]!=0){
      for(int j=0;j<c2YLAccpStNum[iGap][0];j++){
        c2YLGainWeight[iGap][j]=c2YLGainWeightAccpAver[iGap];
      }
    }
    for(UInt_t i=0;i<c2YLAccpStNum[iGap].size()-1;i++){
      c2YLGainWeight[iGap][c2YLAccpStNum[iGap][i]]=c2YLGainWeightAccp[iGap][c2YLAccpStNum[iGap][i]];
      c2YLGainWeight[iGap][c2YLAccpStNum[iGap][i+1]]=c2YLGainWeightAccp[iGap][c2YLAccpStNum[iGap][i+1]];
      for(int j=c2YLAccpStNum[iGap][i]+1;j<c2YLAccpStNum[iGap][i+1];j++){
        c2YLGainWeight[iGap][j]=c2YLGainWeightAccpAver[iGap];
      }
    }
    if(c2YLAccpStNum[iGap][c2YLAccpStNum[iGap].size()-1]!=c2YLNumSt-1){
      for(int j=c2YLAccpStNum[iGap][c2YLAccpStNum[iGap].size()-1]+1;j<c2YLNumSt;j++){
        c2YLGainWeight[iGap][j]=c2YLGainWeightAccpAver[iGap];
      }
    }
  }
  ofstream outPutNumClusterC2YL("number-of-cluster-C2YL");
  for(int iGap=0;iGap<6;iGap++){
    for(int iSt=0;iSt<c2YLNumSt;iSt++){
      outPutNumClusterC2YL<<c2YLNumClusTot[iGap][iSt]<<"  ";
    }
    outPutNumClusterC2YL<<endl;
  }
  outPutNumClusterC2YL.close();

  ofstream outPutGainC2YL("gain-weight-C2YL");
  for(int iGap=0;iGap<6;iGap++){
    for(int i=0;i<c2YLNumSt;i++){
      outPutGainC2YL<<c2YLGainWeight[iGap][i]<<"  ";
    }
    outPutGainC2YL<<endl;
  }
  outPutGainC2YL.close();

  for(int iGap=0;iGap<6;iGap++){
    for(int iSt=0;iSt<c2XLNumSt;iSt++){
      c2XLNumClusTot[iGap][iSt]=c2XLNumClusSing[iGap][iSt]+c2XLNumClusDoub[iGap][iSt]+c2XLNumClusTrip[iGap][iSt];
      if(c2XLNumClusTot[iGap][iSt]>statThreshold){
        c2XLaM2[iGap][iSt]=c2XLSumChargeTripM2[iGap][iSt]/c2XLNumClusTot[iGap][iSt];
        c2XLaM1[iGap][iSt]=(c2XLSumChargeDoubM1[iGap][iSt]+c2XLSumChargeTripM1[iGap][iSt])/c2XLNumClusTot[iGap][iSt];
        c2XLa0[iGap][iSt]=(c2XLSumChargeSing[iGap][iSt]+c2XLSumChargeDoub0[iGap][iSt]+c2XLSumChargeTrip0[iGap][iSt])/c2XLNumClusTot[iGap][iSt];
        c2XLaP1[iGap][iSt]=(c2XLSumChargeDoubP1[iGap][iSt]+c2XLSumChargeTripP1[iGap][iSt])/c2XLNumClusTot[iGap][iSt];
        c2XLaP2[iGap][iSt]=c2XLSumChargeTripP2[iGap][iSt]/c2XLNumClusTot[iGap][iSt];
        c2XLAccpStNum[iGap].push_back(iSt);
      }
      else c2XLRejeStNum[iGap].push_back(iSt);
    }
  }
  for(int iGap=0;iGap<6;iGap++){
    if(c2XLAccpStNum[iGap].size()<=4){
      cout<<"Error: Number of accepted strips with satisfication of statistical threshold is equal to or less than 4 for gap "<<iGap+1<<" of c2XL!"<<endl;
      cout<<"Please reset statistical threshold."<<endl;
      return 0;
    }
    chargeTot=0;
    clusTot=0;
    for(UInt_t i=0;i<c2XLAccpStNum[iGap].size();i++){
      chargeTot+=(c2XLaM2[iGap][c2XLAccpStNum[iGap][i]]+c2XLaM1[iGap][c2XLAccpStNum[iGap][i]]+c2XLa0[iGap][c2XLAccpStNum[iGap][i]]+c2XLaP1[iGap][c2XLAccpStNum[iGap][i]]+c2XLaP2[iGap][c2XLAccpStNum[iGap][i]])*c2XLNumClusTot[iGap][c2XLAccpStNum[iGap][i]];
      clusTot+=c2XLNumClusTot[iGap][c2XLAccpStNum[iGap][i]];
    }
    c2XLChargeAverLevel[iGap]=chargeTot/clusTot;
  }
  for(int iGap=0;iGap<6;iGap++){
    c2XLb[iGap][0][0]=c2XLa0[iGap][0];
    c2XLb[iGap][0][1]=c2XLaP1[iGap][0];
    c2XLb[iGap][0][2]=c2XLaP2[iGap][0];
    c2XLc[iGap][0]=c2XLChargeAverLevel[iGap];
    c2XLb[iGap][1][0]=c2XLaM1[iGap][1];
    c2XLb[iGap][1][1]=c2XLa0[iGap][1];
    c2XLb[iGap][1][2]=c2XLaP1[iGap][1];
    c2XLb[iGap][1][3]=c2XLaP2[iGap][1];
    c2XLc[iGap][1]=c2XLChargeAverLevel[iGap];
    for(UInt_t i=2;i<=c2XLNumSt-3;i++){
      c2XLb[iGap][i][i-2]=c2XLaM2[iGap][i];
      c2XLb[iGap][i][i-1]=c2XLaM1[iGap][i];
      c2XLb[iGap][i][i]=c2XLa0[iGap][i];
      c2XLb[iGap][i][i+1]=c2XLaP1[iGap][i];
      c2XLb[iGap][i][i+2]=c2XLaP2[iGap][i];
      c2XLc[iGap][i]=c2XLChargeAverLevel[iGap];
    }
    c2XLb[iGap][c2XLNumSt-2][c2XLNumSt-4]=c2XLaM2[iGap][c2XLNumSt-2];
    c2XLb[iGap][c2XLNumSt-2][c2XLNumSt-3]=c2XLaM1[iGap][c2XLNumSt-2];
    c2XLb[iGap][c2XLNumSt-2][c2XLNumSt-2]=c2XLa0[iGap][c2XLNumSt-2];
    c2XLb[iGap][c2XLNumSt-2][c2XLNumSt-1]=c2XLaP1[iGap][c2XLNumSt-2];
    c2XLc[iGap][c2XLNumSt-2]=c2XLChargeAverLevel[iGap];
    c2XLb[iGap][c2XLNumSt-1][c2XLNumSt-3]=c2XLaM2[iGap][c2XLNumSt-1];
    c2XLb[iGap][c2XLNumSt-1][c2XLNumSt-2]=c2XLaM1[iGap][c2XLNumSt-1];
    c2XLb[iGap][c2XLNumSt-1][c2XLNumSt-1]=c2XLa0[iGap][c2XLNumSt-1];
    c2XLc[iGap][c2XLNumSt-1]=c2XLChargeAverLevel[iGap];
  }

  vector<int> c2XLRejeStNumCable[6][4];
  for(int iGap=0;iGap<6;iGap++){
    for(UInt_t i=0;i<c2XLRejeStNum[iGap].size();i++){
      for(int j=0;j<4;j++){
        if(C2XStTrans[c2XLRejeStNum[iGap][i]]>=14*j && C2XStTrans[c2XLRejeStNum[iGap][i]]<14*(j+1)) c2XLRejeStNumCable[iGap][j].push_back(c2XLRejeStNum[iGap][i]);
      }
    }
  }
  for(int iGap=0;iGap<6;iGap++){
    TMatrix c2XLB(c2XLAccpStNum[iGap].size(),c2XLAccpStNum[iGap].size()), c2XLW(c2XLAccpStNum[iGap].size(),c2XLAccpStNum[iGap].size());
    for(UInt_t i=0;i<c2XLAccpStNum[iGap].size();i++){
      for(UInt_t j=0;j<c2XLAccpStNum[iGap].size();j++){
        c2XLB(i,j)=c2XLb[iGap][c2XLAccpStNum[iGap][i]][c2XLAccpStNum[iGap][j]];
        for(int k=0;k<4;k++){
          for(UInt_t l=0;l<c2XLRejeStNumCable[iGap][k].size();l++){
            if(C2XStTrans[c2XLAccpStNum[iGap][j]]>=14*k && C2XStTrans[c2XLAccpStNum[iGap][j]]<14*(k+1)) c2XLB(i,j)+=c2XLb[iGap][c2XLAccpStNum[iGap][i]][c2XLRejeStNumCable[iGap][k][l]]/(14-c2XLRejeStNumCable[iGap][k].size());
          }
        }
      }
    }
    for(UInt_t j=0;j<c2XLAccpStNum[iGap].size();j++){
      c2XLW=c2XLB;
      for(UInt_t i=0;i<c2XLAccpStNum[iGap].size();i++){
        c2XLW(i,j)=c2XLc[iGap][i];
      }
      c2XLGainWeightAccp[iGap][c2XLAccpStNum[iGap][j]]=c2XLW.Determinant()/c2XLB.Determinant();
    }
  }
  double c2XLGainWeightAccpSum[6][4]={{0}},c2XLGainWeightAccpAver[6][4];
  int c2XLAccpNumStCable[6][4]={{0}};
  for(int iGap=0;iGap<6;iGap++){
    for(UInt_t i=0;i<c2XLAccpStNum[iGap].size();i++){
      for(int j=0;j<4;j++){
        if(C2XStTrans[c2XLAccpStNum[iGap][i]]>=14*j && C2XStTrans[c2XLAccpStNum[iGap][i]]<14*(j+1)){
          c2XLGainWeightAccpSum[iGap][j]+=c2XLGainWeightAccp[iGap][c2XLAccpStNum[iGap][i]];
          c2XLAccpNumStCable[iGap][j]++;
          break;
        }
      }
    }
    for(int j=0;j<4;j++){
      if(c2XLAccpNumStCable[iGap][j]!=0) c2XLGainWeightAccpAver[iGap][j]=c2XLGainWeightAccpSum[iGap][j]/c2XLAccpNumStCable[iGap][j];
      else c2XLGainWeightAccpAver[iGap][j]=1;
    }
  }
  for(int iGap=0;iGap<6;iGap++){
    if(c2XLAccpStNum[iGap][0]!=0){
      for(int j=0;j<c2XLAccpStNum[iGap][0];j++){
        for(int k=0;k<4;k++){
          if(C2XStTrans[j]>=14*k && C2XStTrans[j]<14*(k+1)){
            c2XLGainWeight[iGap][j]=c2XLGainWeightAccpAver[iGap][k];
            break;
          }
        }
      }
    }
    for(UInt_t i=0;i<c2XLAccpStNum[iGap].size()-1;i++){
      c2XLGainWeight[iGap][c2XLAccpStNum[iGap][i]]=c2XLGainWeightAccp[iGap][c2XLAccpStNum[iGap][i]];
      c2XLGainWeight[iGap][c2XLAccpStNum[iGap][i+1]]=c2XLGainWeightAccp[iGap][c2XLAccpStNum[iGap][i+1]];
      for(int j=c2XLAccpStNum[iGap][i]+1;j<c2XLAccpStNum[iGap][i+1];j++){
        for(int k=0;k<4;k++){
          if(C2XStTrans[j]>=14*k && C2XStTrans[j]<14*(k+1)){
            c2XLGainWeight[iGap][j]=c2XLGainWeightAccpAver[iGap][k];
            break;
          }
        }
      }
    }
    if(c2XLAccpStNum[iGap][c2XLAccpStNum[iGap].size()-1]!=c2XLNumSt-1){
      for(int j=c2XLAccpStNum[iGap][c2XLAccpStNum[iGap].size()-1]+1;j<c2XLNumSt;j++){
        for(int k=0;k<4;k++){
          if(C2XStTrans[j]>=14*k && C2XStTrans[j]<14*(k+1)){
            c2XLGainWeight[iGap][j]=c2XLGainWeightAccpAver[iGap][k];
            break;
          }
        }
      }
    }
  }
  ofstream outPutNumClusterC2XL("number-of-cluster-C2XL");
  for(int iGap=0;iGap<6;iGap++){
    for(int iSt=0;iSt<c2XLNumSt;iSt++){
      outPutNumClusterC2XL<<c2XLNumClusTot[iGap][iSt]<<"  ";
    }
    outPutNumClusterC2XL<<endl;
  }
  outPutNumClusterC2XL.close();

  ofstream outPutGainC2XL("gain-weight-C2XL");
  for(int iGap=0;iGap<6;iGap++){
    for(int i=0;i<c2XLNumSt;i++){
      outPutGainC2XL<<c2XLGainWeight[iGap][i]<<"  ";
    }
    outPutGainC2XL<<endl;
  }
  outPutGainC2XL.close();
  ofstream outPutNormXYC2L("normalized-factor-C2L");
  for(int iGap=0;iGap<6;iGap++){
    outPutNormXYC2L<<c2XLChargeAverLevel[iGap]/c2YLChargeAverLevel[iGap]<<"  ";
  }
  outPutNormXYC2L.close();
  return 0;
};

