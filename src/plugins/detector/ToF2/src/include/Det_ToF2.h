#ifndef __DET_TOF2__
#define __DET_TOF2__
#include "TObject.h"
#include "Plugin.h"
#include "TTree.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TF1.h"
#include "tof2tree.h"
#include <iostream>
#include <iomanip>
#include <vector>
#include <cmath>
#include <utility>
#include <map>
#include "TrekAdcPedestal.h"
class Det_ToF2:public Plugin
{
 private:

  CRTCaliTof2 *treeCali;
  Tof2Info *treeRaw;			/// Input tree with TOF1 raw data
  BeamInfo* treeBeam;
  //
  // Detector parameters set in init file
  //

  // constants 

 public:
  Det_ToF2(TTree *in_,TTree *out_,TFile *inf_, TFile * outf_,TObject *p_);
  virtual ~Det_ToF2();

  // Main ToF2 analysis methods
  Long_t histos();
  Long_t startup();
  Long_t process();
  Long_t done();

  //Fit Pedestal
  Long_t histos_ped();
  Long_t startup_ped();
  Long_t process_ped();
  Long_t done_ped();

  TrekAdcPedestal* adcTof2AO[12];
  TrekAdcPedestal* adcTof2AI[12];
  TrekAdcPedestal* adcTof2BO[12];
  TrekAdcPedestal* adcTof2BI[12];

  TH1D* h1Pedestal;

  //Check Pedestal
  Long_t histos_checkPed();
  Long_t startup_checkPed();
  Long_t process_checkPed();
  Long_t done_checkPed();
  TH1D* h1Adc[48];

  //Conversion
  Long_t histos_conversion();
  Long_t startup_conversion();
  Long_t process_conversion();
  Long_t done_conversion();
  TH1D* h1Multi;
  TH1D* h1TDC;
  TH1D* h1ADC;



  virtual Long_t cmdline(char * cmd);

  
  
  ClassDef(Det_ToF2,1);
};

#endif
