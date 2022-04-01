#ifndef __DET_PGC__
#define __DET_PGC__

//#include <bits/stdc++.h> //conflicts with cint/boost incompatibility
#include <math.h>
#include "TObject.h"
#include "Plugin.h"
#include "TTree.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TF1.h"
#include "pgctree.h"
#include <iostream>
#include <vector>
#include <cmath>
#include <iomanip>
#include <utility>
#include <map>
#include "TCanvas.h"
#include "peakSelect.h"
#include "track.h"
#include "TProfile.h"
//#include "TClonesArray.h"
#include <stdio.h>
#include "trackingE36.h"
#include "pgctree.h"
using namespace std;


class Det_PgC:public Plugin{
 private:
  CRTCaliPgC *treeCali;
  PgcInfo *treeRaw;			/// Input tree with PGC raw data
  BeamInfo* treeBeam;
   trackingE36 *ptr;
public:
  Det_PgC(TTree *in_,TTree *out_,TFile *inf_, TFile * outf_,TObject *p_);
  virtual ~Det_PgC();

  // Main PgC analysis methods
  Long_t histos();
  Long_t startup();
  Long_t process();
  Long_t done();
  TH1D* h1Adc_PGC[12][7];

 TCanvas *c1[12];
TH1D *h1Mom_PgC[12];
TH1D *h1E_PgC[12];
TH1D *h1adc_PgC[12];
TH1D *ADC_PgC[12][7];
TH2D *h1Mom_ADC[12];
TProfile *h1Mom_ADC_TP[12];
//  void Add(const TH1D*);
  TH1D* h_sum;
Long_t start_calibS();
  Long_t calib_peakS();
  Long_t hist_calibS();
 Long_t done_calibS(); 
  //Fit Pedestal
  Long_t histos_ped();
  Long_t startup_ped();
  Long_t process_ped();
  Long_t done_ped();

Double_t pathlength(Double_t z0, Double_t z, Double_t n_z);
Double_t x_p(Double_t x0, Double_t n_x, Double_t t_p);
Double_t y_p(Double_t y0, Double_t n_y, Double_t t_p);
Double_t dp(Double_t p0, Double_t t_p, Double_t Z_A, Double_t I, Double_t rho);
Double_t E_loss(Double_t p0,Double_t p_loss);


  peakSelect* pedADC_PgC[12][7];
  //For ToF2

  Long_t conversionPgC_begin();
  Long_t conversionPgC_histos();
  Long_t conversionPgC_calib();
  Long_t conversionPgC_done();
double pgc_ped[12][7];
 Float_t  myAdcPgC[12][7];
 //double adcPgC[12][7];
  virtual Long_t cmdline(char * cmd);

  
  
  ClassDef(Det_PgC,1);
};

#endif
