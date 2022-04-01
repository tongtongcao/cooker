#ifndef __CHECK_TOF__
#define __CHECK_TOF__
#include "TObject.h"
#include "Plugin.h"
#include "TTree.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TF1.h"
#include "tof1tree.h"
#include <iostream>
#include <iomanip>
#include <vector>
#include <cmath>
#include <utility>
#include <map>
#include <fstream>
class infoBeamKeito{
 public:
  UInt_t indexEvent;//start from 0
  Double_t momentum;
  Double_t theta;
  UInt_t indexGap;//1 - 6
};
class CheckToF:public Plugin
{
 private:

  Tof1Info *treeTof1;			/// Input tree with TOF1 raw data
  Tof2Info *treeTof2;
  BeamInfo* treeBeam;
  //
  // Detector parameters set in init file
  //

  // constants 

 public:
  CheckToF(TTree *in_,TTree *out_,TFile *inf_, TFile * outf_,TObject *p_);
  virtual ~CheckToF();

  // Main analysis methods
  double tofDiff[24];//12*x for tof2 A and B
  Long_t set_tofDiff(int,double);
  Long_t histos();
  Long_t startup();
  Long_t process();
  Long_t done();
  std::vector<infoBeamKeito> listEvent;  
  Long_t set_track_file();
  TH2D* h2TimeVSAngle;
  TH1D* h1TimeDiff;
  TH2D* h2TimeDiffVSGap;
  TH2D* h2TOF1VSTOF2[2][12];
  TH2D* h2MomentumVSCz;
  TH2D* h2DeltaTVSCz[2][12];
  TH2D* h2DeltaTVSCzMu;
  TH2D* h2DeltaTVSCzPi;
  TH2D* h2DeltaTVSCzPositron;
  TString nameTrack;
  std::ifstream fileTrack;
  std::map<UInt_t,Double_t> mapMomentum;
  std::map<UInt_t,UInt_t> mapGap;
  std::map<UInt_t,Double_t> mapCz;








  virtual Long_t cmdline(char * cmd);
  
  
  ClassDef(CheckToF,1);
};

#endif
