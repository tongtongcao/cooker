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
#include "MwpcGainCali.h"
using namespace std;

static string inputFilePath(getenv("mwpc_inputFiles"));
static string inputFile=inputFilePath+"/gapNumTable/run3990-4059";
static ifstream input(inputFile.c_str());
static int runNum, eventNum, gapNumTOF1, gapNumTOF2;

TH2D *histADCVsStBefore[12][6];
TH2D *histADCVsStAfter[12][6];
TH1D *histPosBefore[12][6];
TH1D *histPosAfter[12][6];
TH2D *histDiffVsSumBefore[6][6];
TH2D *histDiffVsSumAfter[6][6];

Long_t Det_MWPC::histos_testGain(){
  ostringstream sstr_name,sstr_title;
  for(int i=0;i<12;i++){
    for(int j=0;j<6;j++){
      sstr_name<<"ADC vs strip number for gap "<<(!(i%2))*6+1+j<<" of C"<<i/4+2<<(char)((i%4)/2+88)<<(char)((i%4%2)*6+76)<<" before gain calibration";
      sstr_title<<"Gap "<<(!(i%2))*6+1+j<<";Strip Number;ADC";
      histADCVsStBefore[i][j]=new TH2D(sstr_name.str().c_str(),sstr_title.str().c_str(),((i%4)/2*(-5)+(i/4)*(!((i%4)/2))+7)*8,-0.5,((i%4)/2*(-5)+(i/4)*(!((i%4)/2))+7)*8-0.5,200,0,4200);
      sstr_name.str("");
      sstr_title.str("");

      sstr_name<<"ADC vs strip number for gap "<<(!(i%2))*6+1+j<<" of C"<<i/4+2<<(char)((i%4)/2+88)<<(char)((i%4%2)*6+76)<<" after gain calibration";
      sstr_title<<"Gap "<<(!(i%2))*6+1+j<<";Strip Number;ADC";
      histADCVsStAfter[i][j]=new TH2D(sstr_name.str().c_str(),sstr_title.str().c_str(),((i%4)/2*(-5)+(i/4)*(!((i%4)/2))+7)*8,-0.5,((i%4)/2*(-5)+(i/4)*(!((i%4)/2))+7)*8-0.5,200,0,4200);
      sstr_name.str("");
      sstr_title.str("");

      sstr_name<<"Position for gap "<<(!(i%2))*6+1+j<<" of C"<<i/4+2<<(char)((i%4)/2+88)<<(char)((i%4%2)*6+76)<<" before gain calibration";
      sstr_title<<"Gap "<<(!(i%2))*6+1+j<<";Position;Counts";
      histPosBefore[i][j]=new TH1D(sstr_name.str().c_str(),sstr_title.str().c_str(),((i%4)/2*(-5)+(i/4)*(!((i%4)/2))+7)*8,0,((i%4)/2*(-5)+(i/4)*(!((i%4)/2))+7)*8);
      sstr_name.str("");
      sstr_title.str("");

      sstr_name<<"Position for gap "<<(!(i%2))*6+1+j<<" of C"<<i/4+2<<(char)((i%4)/2+88)<<(char)((i%4%2)*6+76)<<" after gain calibration";
      sstr_title<<"Gap "<<(!(i%2))*6+1+j<<";Position;Counts";
      histPosAfter[i][j]=new TH1D(sstr_name.str().c_str(),sstr_title.str().c_str(),((i%4)/2*(-5)+(i/4)*(!((i%4)/2))+7)*8,0,((i%4)/2*(-5)+(i/4)*(!((i%4)/2))+7)*8);
      sstr_name.str("");
      sstr_title.str("");
    }
  }

  for(int i=0;i<6;i++){
    for(int j=0;j<6;j++){
      sstr_name<<"Difference vs sum of XY total-charge for gap "<<(!(i%2))*6+1+j<<" of C"<<i/2+2<<(char)((i%2)*6+76)<<" before gain calibration";
      sstr_title<<"Gap "<<(!(i%2))*6+1+j<<";TC_{X}+TC_{Y};TC_{X}-TC_{Y}";
      histDiffVsSumBefore[i][j]=new TH2D(sstr_name.str().c_str(),sstr_title.str().c_str(),200,0,6000,200,-3000,3000);
      sstr_name.str("");
      sstr_title.str("");

      sstr_name<<"Difference vs sum of XY total-charge for gap "<<(!(i%2))*6+1+j<<" of C"<<i/2+2<<(char)((i%2)*6+76)<<" after gain calibration";
      sstr_title<<"Gap "<<(!(i%2))*6+1+j<<";TC_{X}+TC_{Y};TC_{X}-TC_{Y}";
      histDiffVsSumAfter[i][j]=new TH2D(sstr_name.str().c_str(),sstr_title.str().c_str(),200,0,6000,200,-3000,3000);
      sstr_name.str("");
      sstr_title.str("");
    }
  }

  return 0;
}

static vector <double> gain[12][6];
static double norm[6][6];
Long_t Det_MWPC::startup_testGain(){
  getBranchObject("CaliMwpcInfo",(TObject **) &treeCali);
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

  for(int i=0;i<6;i++){
    norm[0][i]=mwpcXYNormFactorC2L[i];
    norm[1][i]=mwpcXYNormFactorC2R[i];
    norm[2][i]=mwpcXYNormFactorC3L[i];
    norm[3][i]=mwpcXYNormFactorC3R[i];
    norm[4][i]=mwpcXYNormFactorC4L[i];
    norm[5][i]=mwpcXYNormFactorC4R[i];
  }

  return 0;
};

static vector<double> stCharge;
static double totalChargeX,totalChargeY;
Long_t Det_MWPC::process_testGain(){
  input>>runNum>>eventNum>>gapNumTOF1>>gapNumTOF2;
  if((UInt_t)runNum!=treeCali->run || (UInt_t)eventNum!=treeCali->event){
    cout<<"Error: Run number or event number are not consistent between two files!"<<endl;
    cout<<"For MWPC file: "<<"Run number - "<<treeCali->run<<"and  Event number - "<<treeCali->event<<endl;
    cout<<"For TRIUMF file: "<<"Run number - "<<runNum<<"and  Event number - "<<eventNum<<endl;
  }
 
  MwpcEvent *event[12];
  for(int i=0;i<12;i++){
    event[i]=new MwpcEvent(((i%4)/2*(-5)+(i/4)*(!((i%4)/2))+7)*8,!((i/4==0)&&((i%4/2))==0));
    if(gapNumTOF2!=-1){
      if((i==4 || i==6) && gapNumTOF2==7){
        if(treeCali->thTCCutNumClusterC3XLGap7Swap==1 && treeCali->thTCCutNumCluster[6]==1 && i==4){
          stCharge.clear();
          for(UInt_t j=0;j<treeCali->thTCCutClusterStNumC3XLGap7Swap[0].size();j++){
            histADCVsStBefore[4][0]->Fill(treeCali->thTCCutClusterStNumC3XLGap7Swap[0][j],treeCali->thTCCutClusterStChargeC3XLGap7Swap[0][j]);
            histADCVsStAfter[4][0]->Fill(treeCali->thTCCutClusterStNumC3XLGap7Swap[0][j],treeCali->thTCCutClusterStChargeC3XLGap7Swap[0][j]*gain[4][0][treeCali->thTCCutClusterStNumC3XLGap7Swap[0][j]]);
            stCharge.push_back(treeCali->thTCCutClusterStChargeC3XLGap7Swap[0][j]*gain[4][0][treeCali->thTCCutClusterStNumC3XLGap7Swap[0][j]]);
          }
          histPosBefore[4][0]->Fill(treeCali->thTCCutClusterPosCentroidC3XLGap7Swap[0]);
          histPosAfter[4][0]->Fill(event[4]->getEachClusterPositionCentroid(treeCali->thTCCutClusterStNumC3XLGap7Swap[0],stCharge));
          stCharge.clear();
          for(UInt_t j=0;j<treeCali->thTCCutClusterStNum[6][0].size();j++){
            histADCVsStBefore[6][0]->Fill(treeCali->thTCCutClusterStNum[6][0][j],treeCali->thTCCutClusterStCharge[6][0][j]);
            histADCVsStAfter[6][0]->Fill(treeCali->thTCCutClusterStNum[6][0][j],treeCali->thTCCutClusterStCharge[6][0][j]*gain[6][0][treeCali->thTCCutClusterStNum[6][0][j]]);
            stCharge.push_back(treeCali->thTCCutClusterStCharge[6][0][j]*gain[6][0][treeCali->thTCCutClusterStNum[6][0][j]]);
          }
          histPosBefore[6][0]->Fill(treeCali->thTCCutClusterPosCentroid[6][0]);
          histPosAfter[6][0]->Fill(event[i]->getEachClusterPositionCentroid(treeCali->thTCCutClusterStNum[6][0],stCharge));
        }
      }
      else if((i==8 || i==10) && gapNumTOF2==7){
        if(treeCali->thTCCutNumClusterC4XLGap7Swap==1 && treeCali->thTCCutNumCluster[10]==1 && i==8){
          stCharge.clear();
          for(UInt_t j=0;j<treeCali->thTCCutClusterStNumC4XLGap7Swap[0].size();j++){
            histADCVsStBefore[8][0]->Fill(treeCali->thTCCutClusterStNumC4XLGap7Swap[0][j],treeCali->thTCCutClusterStChargeC4XLGap7Swap[0][j]);
            histADCVsStAfter[8][0]->Fill(treeCali->thTCCutClusterStNumC4XLGap7Swap[0][j],treeCali->thTCCutClusterStChargeC4XLGap7Swap[0][j]*gain[8][0][treeCali->thTCCutClusterStNumC4XLGap7Swap[0][j]]);
            stCharge.push_back(treeCali->thTCCutClusterStChargeC4XLGap7Swap[0][j]*gain[8][0][treeCali->thTCCutClusterStNumC4XLGap7Swap[0][j]]);
          }
          histPosBefore[8][0]->Fill(treeCali->thTCCutClusterPosCentroidC4XLGap7Swap[0]);
          histPosAfter[8][0]->Fill(event[8]->getEachClusterPositionCentroid(treeCali->thTCCutClusterStNumC4XLGap7Swap[0],stCharge));
          stCharge.clear();
          for(UInt_t j=0;j<treeCali->thTCCutClusterStNum[10][0].size();j++){
            histADCVsStBefore[10][0]->Fill(treeCali->thTCCutClusterStNum[10][0][j],treeCali->thTCCutClusterStCharge[10][0][j]);
            histADCVsStAfter[10][0]->Fill(treeCali->thTCCutClusterStNum[10][0][j],treeCali->thTCCutClusterStCharge[10][0][j]*gain[10][0][treeCali->thTCCutClusterStNum[10][0][j]]);
            stCharge.push_back(treeCali->thTCCutClusterStCharge[10][0][j]*gain[10][0][treeCali->thTCCutClusterStNum[10][0][j]]);
          }
          histPosBefore[10][0]->Fill(treeCali->thTCCutClusterPosCentroid[10][0]);
          histPosAfter[10][0]->Fill(event[i]->getEachClusterPositionCentroid(treeCali->thTCCutClusterStNum[10][0],stCharge));
        }
      }
      else if((i==8 || i==10) && gapNumTOF2==11){
        if(treeCali->thTCCutNumClusterC4XLGap11Swap==1 && treeCali->thTCCutNumCluster[10]==1 && i==8){
          stCharge.clear();
          for(UInt_t j=0;j<treeCali->thTCCutClusterStNumC4XLGap11Swap[0].size();j++){
            histADCVsStBefore[8][4]->Fill(treeCali->thTCCutClusterStNumC4XLGap11Swap[0][j],treeCali->thTCCutClusterStChargeC4XLGap11Swap[0][j]);
            histADCVsStAfter[8][4]->Fill(treeCali->thTCCutClusterStNumC4XLGap11Swap[0][j],treeCali->thTCCutClusterStChargeC4XLGap11Swap[0][j]*gain[8][4][treeCali->thTCCutClusterStNumC4XLGap11Swap[0][j]]);
            stCharge.push_back(treeCali->thTCCutClusterStChargeC4XLGap11Swap[0][j]*gain[8][4][treeCali->thTCCutClusterStNumC4XLGap11Swap[0][j]]);
          } 
          histPosBefore[8][4]->Fill(treeCali->thTCCutClusterPosCentroidC4XLGap11Swap[0]);
          histPosAfter[8][4]->Fill(event[8]->getEachClusterPositionCentroid(treeCali->thTCCutClusterStNumC4XLGap11Swap[0],stCharge));
          stCharge.clear();
          for(UInt_t j=0;j<treeCali->thTCCutClusterStNum[10][0].size();j++){
            histADCVsStBefore[10][4]->Fill(treeCali->thTCCutClusterStNum[10][0][j],treeCali->thTCCutClusterStCharge[10][0][j]);
            histADCVsStAfter[10][4]->Fill(treeCali->thTCCutClusterStNum[10][0][j],treeCali->thTCCutClusterStCharge[10][0][j]*gain[10][4][treeCali->thTCCutClusterStNum[10][0][j]]);
            stCharge.push_back(treeCali->thTCCutClusterStCharge[10][0][j]*gain[10][4][treeCali->thTCCutClusterStNum[10][0][j]]);
          }
          histPosBefore[10][4]->Fill(treeCali->thTCCutClusterPosCentroid[10][0]);
          histPosAfter[10][4]->Fill(event[i]->getEachClusterPositionCentroid(treeCali->thTCCutClusterStNum[10][0],stCharge));
        } 
      }
      else{
        if(treeCali->thTCCutNumCluster[i]==1 && treeCali->thTCCutNumCluster[i+2-i%4/2*4]==1){ 
          if(gapNumTOF2>=(!(i%2))*6+1 && gapNumTOF2<=(!(i%2))*6+6){
            stCharge.clear();
            for(UInt_t j=0;j<treeCali->thTCCutClusterStNum[i][0].size();j++){
              histADCVsStBefore[i][gapNumTOF2-(!(i%2))*6-1]->Fill(treeCali->thTCCutClusterStNum[i][0][j],treeCali->thTCCutClusterStCharge[i][0][j]);
              histADCVsStAfter[i][gapNumTOF2-(!(i%2))*6-1]->Fill(treeCali->thTCCutClusterStNum[i][0][j],treeCali->thTCCutClusterStCharge[i][0][j]*gain[i][gapNumTOF2-(!(i%2))*6-1][treeCali->thTCCutClusterStNum[i][0][j]]);
              stCharge.push_back(treeCali->thTCCutClusterStCharge[i][0][j]*gain[i][gapNumTOF2-(!(i%2))*6-1][treeCali->thTCCutClusterStNum[i][0][j]]);
            }
            histPosBefore[i][gapNumTOF2-(!(i%2))*6-1]->Fill(treeCali->thTCCutClusterPosCentroid[i][0]);
            histPosAfter[i][gapNumTOF2-(!(i%2))*6-1]->Fill(event[i]->getEachClusterPositionCentroid(treeCali->thTCCutClusterStNum[i][0],stCharge));
          }
        }
      }
    }
  }
      



  for(int i=0;i<6;i++){
    if(gapNumTOF2!=-1){
      if(i==2 && gapNumTOF2==7){
        if(treeCali->thTCCutNumClusterC3XLGap7Swap==1 && treeCali->thTCCutNumCluster[i*2-i%2+2]==1){
          totalChargeX=0;
          for(UInt_t j=0;j<treeCali->thTCCutClusterStNumC3XLGap7Swap[0].size();j++){
            totalChargeX+=treeCali->thTCCutClusterStChargeC3XLGap7Swap[0][j]*gain[i*2-i%2][gapNumTOF2-(!(i%2))*6-1][treeCali->thTCCutClusterStNumC3XLGap7Swap[0][j]]/norm[i][gapNumTOF2-(!(i%2))*6-1];
          }
          totalChargeY=0;
          for(UInt_t j=0;j<treeCali->thTCCutClusterStNum[i*2-i%2+2][0].size();j++){
            totalChargeY+=treeCali->thTCCutClusterStCharge[i*2-i%2+2][0][j]*gain[i*2-i%2+2][gapNumTOF2-(!(i%2))*6-1][treeCali->thTCCutClusterStNum[i*2-i%2+2][0][j]];
          }
          histDiffVsSumBefore[i][gapNumTOF2-(!(i%2))*6-1]->Fill(treeCali->thTCCutClusterTCC3XLGap7Swap[0]+treeCali->thTCCutClusterTC[i*2-i%2+2][0],treeCali->thTCCutClusterTCC3XLGap7Swap[0]-treeCali->thTCCutClusterTC[i*2-i%2+2][0]);
          histDiffVsSumAfter[i][gapNumTOF2-(!(i%2))*6-1]->Fill(totalChargeX+totalChargeY,totalChargeX-totalChargeY);
        }
      }
      if(i==4 && gapNumTOF2==7){
        if(treeCali->thTCCutNumClusterC4XLGap7Swap==1 && treeCali->thTCCutNumCluster[i*2-i%2+2]==1){
          totalChargeX=0;
          for(UInt_t j=0;j<treeCali->thTCCutClusterStNumC4XLGap7Swap[0].size();j++){
            totalChargeX+=treeCali->thTCCutClusterStChargeC4XLGap7Swap[0][j]*gain[i*2-i%2][gapNumTOF2-(!(i%2))*6-1][treeCali->thTCCutClusterStNumC4XLGap7Swap[0][j]]/norm[i][gapNumTOF2-(!(i%2))*6-1];
          }
          totalChargeY=0;
          for(UInt_t j=0;j<treeCali->thTCCutClusterStNum[i*2-i%2+2][0].size();j++){
            totalChargeY+=treeCali->thTCCutClusterStCharge[i*2-i%2+2][0][j]*gain[i*2-i%2+2][gapNumTOF2-(!(i%2))*6-1][treeCali->thTCCutClusterStNum[i*2-i%2+2][0][j]];
          }
          histDiffVsSumBefore[i][gapNumTOF2-(!(i%2))*6-1]->Fill(treeCali->thTCCutClusterTCC4XLGap7Swap[0]+treeCali->thTCCutClusterTC[i*2-i%2+2][0],treeCali->thTCCutClusterTCC4XLGap7Swap[0]-treeCali->thTCCutClusterTC[i*2-i%2+2][0]);
          histDiffVsSumAfter[i][gapNumTOF2-(!(i%2))*6-1]->Fill(totalChargeX+totalChargeY,totalChargeX-totalChargeY);
        }
      }
      if(i==4 && gapNumTOF2==11){
        if(treeCali->thTCCutNumClusterC4XLGap11Swap==1 && treeCali->thTCCutNumCluster[i*2-i%2+2]==1){
          totalChargeX=0;
          for(UInt_t j=0;j<treeCali->thTCCutClusterStNumC4XLGap11Swap[0].size();j++){
            totalChargeX+=treeCali->thTCCutClusterStChargeC4XLGap11Swap[0][j]*gain[i*2-i%2][gapNumTOF2-(!(i%2))*6-1][treeCali->thTCCutClusterStNumC4XLGap11Swap[0][j]]/norm[i][gapNumTOF2-(!(i%2))*6-1];
          }
          totalChargeY=0;
          for(UInt_t j=0;j<treeCali->thTCCutClusterStNum[i*2-i%2+2][0].size();j++){
            totalChargeY+=treeCali->thTCCutClusterStCharge[i*2-i%2+2][0][j]*gain[i*2-i%2+2][gapNumTOF2-(!(i%2))*6-1][treeCali->thTCCutClusterStNum[i*2-i%2+2][0][j]];
          }
          histDiffVsSumBefore[i][gapNumTOF2-(!(i%2))*6-1]->Fill(treeCali->thTCCutClusterTCC4XLGap11Swap[0]+treeCali->thTCCutClusterTC[i*2-i%2+2][0],treeCali->thTCCutClusterTCC4XLGap11Swap[0]-treeCali->thTCCutClusterTC[i*2-i%2+2][0]);
          histDiffVsSumAfter[i][gapNumTOF2-(!(i%2))*6-1]->Fill(totalChargeX+totalChargeY,totalChargeX-totalChargeY);
        }
      }
      else{
        if(treeCali->thTCCutNumCluster[i*2-i%2]==1 && treeCali->thTCCutNumCluster[i*2-i%2+2]==1 && gapNumTOF2!=-1){
          if(gapNumTOF2>=(!(i%2))*6+1 && gapNumTOF2<=(!(i%2))*6+6){
            totalChargeX=0;
            for(UInt_t j=0;j<treeCali->thTCCutClusterStNum[i*2-i%2][0].size();j++){
              totalChargeX+=treeCali->thTCCutClusterStCharge[i*2-i%2][0][j]*gain[i*2-i%2][gapNumTOF2-(!(i%2))*6-1][treeCali->thTCCutClusterStNum[i*2-i%2][0][j]]/norm[i][gapNumTOF2-(!(i%2))*6-1];
            }
            totalChargeY=0;
            for(UInt_t j=0;j<treeCali->thTCCutClusterStNum[i*2-i%2+2][0].size();j++){
              totalChargeY+=treeCali->thTCCutClusterStCharge[i*2-i%2+2][0][j]*gain[i*2-i%2+2][gapNumTOF2-(!(i%2))*6-1][treeCali->thTCCutClusterStNum[i*2-i%2+2][0][j]];
            }
            histDiffVsSumBefore[i][gapNumTOF2-(!(i%2))*6-1]->Fill(treeCali->thTCCutClusterTC[i*2-i%2][0]+treeCali->thTCCutClusterTC[i*2-i%2+2][0],treeCali->thTCCutClusterTC[i*2-i%2][0]-treeCali->thTCCutClusterTC[i*2-i%2+2][0]);
            histDiffVsSumAfter[i][gapNumTOF2-(!(i%2))*6-1]->Fill(totalChargeX+totalChargeY,totalChargeX-totalChargeY);
          }
        }
      }
    }
  }
 

  return 0;
}

Long_t Det_MWPC::done_testGain(){
  ostringstream sstr;
  TCanvas *canvADCVsSt[12];
  TCanvas *canvPos[12];
  TLegend *legend=new TLegend(0.7,0.7,0.9,0.9);
  legend->AddEntry(histPosBefore[0][0],"Before","l");
  legend->AddEntry(histPosAfter[0][0],"After","l");
  for(int i=0;i<12;i++){
    sstr<<"ADC vs strip number for C"<<i/4+2<<(char)((i%4)/2+88)<<(char)((i%4%2)*6+76);
    canvADCVsSt[i]=new TCanvas(sstr.str().c_str(),sstr.str().c_str(),2000,900);
    sstr.str("");
    canvADCVsSt[i]->Divide(4,3);
    for(int j=0;j<6;j++){
      canvADCVsSt[i]->cd(j*2-j%2+1);
      gPad->SetLogz();
      histADCVsStBefore[i][j]->Draw("colz");
      canvADCVsSt[i]->cd(j*2-j%2+3);
      gPad->SetLogz();
      histADCVsStAfter[i][j]->Draw("colz");
    }
    sstr<<"ADC vs strip number for C"<<i/4+2<<(char)((i%4)/2+88)<<(char)((i%4%2)*6+76)<<".pdf";
    canvADCVsSt[i]->SaveAs(sstr.str().c_str());
    sstr.str("");

    sstr<<"Position for C"<<i/4+2<<(char)((i%4)/2+88)<<(char)((i%4%2)*6+76);
    canvPos[i]=new TCanvas(sstr.str().c_str(),sstr.str().c_str(),1000,900);
    sstr.str("");
    canvPos[i]->Divide(2,3);
    for(int j=0;j<6;j++){
      canvPos[i]->cd(j+1);
      histPosBefore[i][j]->Draw();
      histPosAfter[i][j]->SetLineColor(kRed);
      histPosAfter[i][j]->Draw("same");
      legend->Draw();
    }
    sstr<<"Position for C"<<i/4+2<<(char)((i%4)/2+88)<<(char)((i%4%2)*6+76)<<".pdf";
    canvPos[i]->SaveAs(sstr.str().c_str());
    sstr.str("");
  }

  TCanvas *canvDiffVsSum[6];
  for(int i=0;i<6;i++){
    sstr<<"Diff vs sum of XY total-charge for C"<<i/2+2<<(char)((i%2)*6+76);
    canvDiffVsSum[i]=new TCanvas(sstr.str().c_str(),sstr.str().c_str(),2000,900);
    sstr.str("");
    canvDiffVsSum[i]->Divide(4,3);
    for(int j=0;j<6;j++){
      canvDiffVsSum[i]->cd(j*2-j%2+1);
      gPad->SetLogz();
      histDiffVsSumBefore[i][j]->Draw("colz");
      canvDiffVsSum[i]->cd(j*2-j%2+3);
      gPad->SetLogz();
      histDiffVsSumAfter[i][j]->Draw("colz");
    }
    sstr<<"Diff vs sum of XY total-charge for C"<<i/2+2<<(char)((i%2)*6+76)<<".pdf";
    canvDiffVsSum[i]->SaveAs(sstr.str().c_str());
    sstr.str("");
  }


  return 0;
};

