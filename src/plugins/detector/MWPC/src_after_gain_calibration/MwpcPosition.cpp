#include "MwpcPosition.h"
#include "TF1.h"
#include <vector>
#include "TCanvas.h"
#include "TLine.h"
#include "TSpectrum.h"
#include "TLatex.h"
#include <algorithm>
#include <iostream>
#include <sstream>
using namespace std;

void MwpcPosition::defineHist(){
  ostringstream name,title;
  name<<detectorName<<"_Event"<<eventNumber;
  title<<detectorName<<"_Event"<<eventNumber<<";X;Y";
  positionHist=new TH2D(name.str().c_str(),title.str().c_str(),numBinX,0,numBinX,16,0,16);
  name.str("");
  title.str("");
  for(int i=0;i<numBinX;i++){
    for(int j=0;j<16;j++){
      positionHist->SetBinContent(i+1,j+1,0);
    }
  }
}

void MwpcPosition::setBinContent(vector<vector<int> > stNumX, vector<vector<double> > stChargeX, vector<double> pairTotalChargeX,vector<vector<int> > stNumY, vector<vector<double> > stChargeY, vector<double> pairTotalChargeY){
  for(UInt_t i=0;i<stNumX.size();i++){
    for(UInt_t m=0;m<stNumX[i].size();m++){
      for(UInt_t n=0;n<stNumY[i].size();n++){
        positionHist->SetBinContent(stNumX[i][m]+1,stNumY[i][n]+1,(pairTotalChargeX[i]+pairTotalChargeY[i])*stChargeX[i][m]/pairTotalChargeX[i]*stChargeY[i][n]/pairTotalChargeY[i]);
      }
     }
   }
}

vector<TMarker*> MwpcPosition::markPosition(vector<double> posX,vector<double> posErrX,vector<double> posY,vector<double> posErrY){
  TMarker *marker;
  vector<TMarker*> vect_marker;
  vect_marker.clear();
  for(UInt_t i=0;i<posX.size();i++){
    marker=new TMarker(posX[i],posY[i],20);
    vect_marker.push_back(marker);
  }
  return vect_marker;
}

vector<TLine*> MwpcPosition::markPositionErrorLeft(vector<double> posX,vector<double> posErrX,vector<double> posY,vector<double> posErrY){
  TLine *line;
  vector<TLine*> vect_line;
  vect_line.clear();
  for(UInt_t i=0;i<posX.size();i++){
    line=new TLine(posX[i]-posErrX[i],posY[i]-posErrY[i],posX[i]-posErrX[i],posY[i]+posErrY[i]);
    vect_line.push_back(line);
  }
  return vect_line;
}

vector<TLine*> MwpcPosition::markPositionErrorRight(vector<double> posX,vector<double> posErrX,vector<double> posY,vector<double> posErrY){
  TLine *line;
  vector<TLine*> vect_line;
  vect_line.clear();
  for(UInt_t i=0;i<posX.size();i++){
    line=new TLine(posX[i]+posErrX[i],posY[i]-posErrY[i],posX[i]+posErrX[i],posY[i]+posErrY[i]);
    vect_line.push_back(line);
  }
  return vect_line;
}

vector<TLine*> MwpcPosition::markPositionErrorTop(vector<double> posX,vector<double> posErrX,vector<double> posY,vector<double> posErrY){
  TLine *line;
  vector<TLine*> vect_line;
  vect_line.clear();
  for(UInt_t i=0;i<posX.size();i++){
    line=new TLine(posX[i]-posErrX[i],posY[i]+posErrY[i],posX[i]+posErrX[i],posY[i]+posErrY[i]);
    vect_line.push_back(line);
  }
  return vect_line;
}

vector<TLine*> MwpcPosition::markPositionErrorBottom(vector<double> posX,vector<double> posErrX,vector<double> posY,vector<double> posErrY){
  TLine *line;
  vector<TLine*> vect_line;
  vect_line.clear();
  for(UInt_t i=0;i<posX.size();i++){
    line=new TLine(posX[i]-posErrX[i],posY[i]-posErrY[i],posX[i]+posErrX[i],posY[i]-posErrY[i]);
    vect_line.push_back(line);
  }
  return vect_line;
}
