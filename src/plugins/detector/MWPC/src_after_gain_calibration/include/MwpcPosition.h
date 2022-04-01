#ifndef __MWPC_POSITION__
#define __MWPC_POSITION__
#include <string>
#include <TH2D.h>
#include <TDirectory.h>
#include <TList.h>
#include <vector>
#include "TMarker.h"
#include "TLine.h"
using namespace std;

class MwpcPosition{
 private:
  int eventNumber;
  std::string detectorName;
  TH2D* positionHist;
  int numBinX;
 public:
 MwpcPosition(int eventNum, std::string name,int nBinX):eventNumber(eventNum),detectorName(name),numBinX(nBinX),positionHist(0){} 
 void defineHist();
 void setBinContent(vector<vector<int> > stNumX, vector<vector<double> > stChargeX, vector<double> pairTotalChargeX, vector<vector<int> > stNumY, vector<vector<double> > stChargeY, vector<double> pairTotalChargeY);
 vector<TMarker*> markPosition(vector<double> posX,vector<double> posErrX,vector<double> posY,vector<double> posErrY);
 vector<TLine*> markPositionErrorLeft(vector<double> posX,vector<double> posErrX,vector<double> posY,vector<double> posErrY);
 vector<TLine*> markPositionErrorRight(vector<double> posX,vector<double> posErrX,vector<double> posY,vector<double> posErrY);
 vector<TLine*> markPositionErrorTop(vector<double> posX,vector<double> posErrX,vector<double> posY,vector<double> posErrY);
 vector<TLine*> markPositionErrorBottom(vector<double> posX,vector<double> posErrX,vector<double> posY,vector<double> posErrY);
 UInt_t getEventNumber(){return eventNumber;}  // Return event number
 int getNBinX(){return numBinX;}  //return number of bins of histogram in X-axis == number of strips in X direction
 TH2D* getHist(){return positionHist;} // Return histogram 


};
#endif
