#ifndef __DET_AC__
#define __DET_AC__

//#include <bits/stdc++.h> //conflicts with cint/boost incompatibility

#include "TObject.h"
#include "Plugin.h"
#include "TTree.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TF1.h"
#include "TCanvas.h"
#include "actree.h"
#include <iostream>
#include <vector>
#include <cmath>
#include <iomanip>
#include <utility>
#include <map>
#include "TrekAcAdcPedestal.h"
#include "genericPedMon.h"
#include "peakSelect.h"
//#include "trackingE36.h"
using namespace std;

class Det_Ac:public Plugin{
private:
  CRTCaliAc *treeCali;
  AcInfo *treeRaw;			/// Input tree with AC raw data
  BeamInfo* treeBeam;
  // second rootfile:
//  trackingE36 *ptr;
  
//trackingE36 *ptr;
  //Double_t pVertPosition[50];

  //
  // Detector parameters set in init file
  //

  // constants 

  double mPedestal;
  double mError;
  double mSigma; 
  double mSigmaError;

public:
  Det_Ac(TTree *in_,TTree *out_,TFile *inf_, TFile * outf_,TObject *p_);
  virtual ~Det_Ac();

  // Main Ac analysis methods
  Long_t histos();
  Long_t startup();
  Long_t process();
  Long_t hist_calib();  // will use method to calibrate ac ADC
  Long_t done();

  TH1D* h1adc_ACu[12];
  TH1D* h1adc_ACd[12];

  // gain matching histos
  TH1D* gainM_ACu[12];
  TH1D* gainM_ACd[12];

  // calib matching histos
  TH1D* calib_ACu[12];
  TH1D* calib_ACd[12];

  // calib and gain-matching variables
  TTree* tree;
  Double_t fpedAdcU[12]; //tree variables
  Double_t fpedAdcD[12];
   TH1D *trackMom;

  /*double (&fpedAdcU_ref)[12] = fpedAdcU;
  double (&fpedAdcD_ref)[12] = fpedAdcD;
*/
  double fgainU[12];
  double fgainD[12];
  
  TH1D* h1tdc_ACu[12];
  TH1D* h1tdc_ACd[12];
  TH1D* h_sum;
  //Fit Pedestal
  Long_t histos_ped();
  Long_t startup_ped();
  Long_t start_calib();
  Long_t calib_peak();
  Long_t calib();
  Long_t process_ped();
  Long_t done_ped();

//conversionAc
  Long_t conversionAc_begin();
  Long_t conversionAc_histos();
  Long_t conversionAc_calib(); 
  Long_t conversionAc_done(); 

//Int_t Run, D1, D2, D3, D4, D5, D6, D7, D8, D9, D10, D11, D12;












  TrekAcAdcPedestal* adcAcU[12];
  TrekAcAdcPedestal* adcAcD[12];
  genericPedMon* adcAcU1[12];
  genericPedMon* adcAcD1[12];
  genericPedMon* adcAcU2[12];
  genericPedMon* adcAcD2[12];
  genericPedMon* adcAcU3[12];
  genericPedMon* adcAcD3[12];
  genericPedMon* adcAcU4[12];
  genericPedMon* adcAcD4[12];
  genericPedMon* adcAcU5[12];
  genericPedMon* adcAcD5[12];
  genericPedMon* adcAcU6[12];
  genericPedMon* adcAcD6[12];
  genericPedMon* adcAcU7[12];
  genericPedMon* adcAcD7[12];

  peakSelect* pedADC_AcU[12];
  peakSelect* pedADC_AcD[12];

  TH1D* h1Pedestal; 
  TH1D* pedSub;
  //monitor ped
  TH1D* mpedU[30][12];
  TH1D* mpedD[30][12];
  TH1F* mH[30][12];

  //Check Pedestal
  Long_t histos_checkPed();
  Long_t monitor_ped();
  Long_t monitor_begin();
  Long_t startup_checkPed();
  Long_t process_checkPed();
  Long_t process_monitor();
  Long_t done_checkPed();
  Long_t done_calib();
  Long_t monitor_done();
  TH1D* h1Adc[24];
  Float_t    myAdcAcU[12]; 
  Float_t    myAdcAcD[12]; 
 virtual Long_t cmdline(char * cmd);
  
  ClassDef(Det_Ac,1);
};

#endif
