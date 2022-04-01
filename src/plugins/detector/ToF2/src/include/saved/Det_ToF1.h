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
class Det_ToF1:public Plugin
{
 private:

  CRTCaliTof1 *treeCali;
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


  virtual Long_t cmdline(char * cmd);

  
  
  ClassDef(Det_ToF1,1);
};

#endif
