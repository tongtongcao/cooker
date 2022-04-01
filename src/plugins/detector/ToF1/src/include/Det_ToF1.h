#ifndef __DET_TOF1__
#define __DET_TOF1__
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
#include "TrekAdcPedestal.h"
#include "TrekTimewalk.h"
class Det_ToF1:public Plugin
{
 private:

  CRTCaliTof1 *treeCali;//output tree
  Tof1Info *treeRaw;			/// Input tree with TOF1 raw data
  BeamInfo* treeBeam;
  //
  // Detector parameters set in init file
  //

  // constants 

 public:
  Det_ToF1(TTree *in_,TTree *out_,TFile *inf_, TFile * outf_,TObject *p_);
  virtual ~Det_ToF1();

  // Main ToF1 analysis methods
  Long_t histos();
  Long_t startup();
  Long_t process();
  Long_t done();

  TH2D* h2AdcVSTdc_u[12];
  TH2D* h2AdcVSTdc_d[12];

  //Fit Pedestal
  Long_t histos_ped();
  Long_t startup_ped();
  Long_t process_ped();
  Long_t done_ped();

  TrekAdcPedestal* adcTof1U[12];
  TrekAdcPedestal* adcTof1D[12];

  TH1D* h1Pedestal;

  //Check Pedestal
  Long_t histos_checkPed();
  Long_t startup_checkPed();
  Long_t process_checkPed();
  Long_t done_checkPed();
  TH1D* h1Adc[24];

  //Conversion
  Long_t histos_conversion();
  Long_t startup_conversion();
  Long_t process_conversion();
  Long_t done_conversion();
  TH1D* h1Multi;
  TH2D* h2MultiADCVSTDC;
  TH1D* h1TDC;
  TH1D* h1ADC;
  virtual Long_t cmdline(char * cmd);

  //Left and Right
  Long_t histos_lar();
  Long_t startup_lar();
  Long_t process_lar();
  Long_t done_lar();

  TH1D* h1Diff[12];
  TH2D* h2GapVSTDiff;
  TH2D* h2GapVSTDiffCorr;
  Double_t shiftLeftRight[12];

  //veff
  Long_t histos_veff();
  Long_t startup_veff();
  Long_t process_veff();
  Long_t done_veff();

  TH2D* h2HitVSDiff[12];
  Double_t zGap[500000];
  int indexGap[500000];
  Double_t veff[12];
  Double_t offset[12];
  unsigned countGap[12];  

  //time walk
  Long_t histos_timewalk();
  Long_t startup_timewalk();
  Long_t process_timewalk();
  Long_t done_timewalk();

  TH2D* h2DiffVSAdc[12][2];
  Double_t Timewalk[12][2];
  Double_t listPedestal[24];
  TrekTimewalk* tof1Timewalk[12][2];
  ClassDef(Det_ToF1,1);
};

#endif
