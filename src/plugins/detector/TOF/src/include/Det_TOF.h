#ifndef __DET_TOF__
#define __DET_TOF__
#include "TObject.h"
#include "Plugin.h"
#include "TTree.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TF1.h"
#include "tof1tree.h"
#include "tof2tree.h"
#include <iostream>
#include <iomanip>
#include <vector>
#include <cmath>
#include <utility>
#include <map>
#include "trackingE36.h"
#include "TCanvas.h"

class Det_TOF:public Plugin
{
 private:

  CRTCaliTof1 *treeTof1;//input tree tof1
  CRTCaliTof2 *treeTof2;//input tree tof2
  BeamInfo* treeBeam;//input tree beam
  trackingE36 *ptr;

  //
  // Detector parameters set in init file
  //

  // constants 

 public:
  Det_TOF(TTree *in_,TTree *out_,TFile *inf_, TFile * outf_,TObject *p_);
  virtual ~Det_TOF();

  // Main TOF analysis methods
  Long_t histos();
  Long_t startup();
  Long_t process();
  Long_t done();

  //Conversion
  Long_t histos_conversion();
  Long_t startup_conversion();
  Long_t process_conversion();
  Long_t done_conversion();
  TH1D* h1Multi;
  TH1D* h1MultiLow;
  TH1D* h1MultiUp;
  TH1D* h1MultiTotal;
   TH1D *trackM;
   TH1D *trackBeta;

/*Double_t path_length;
Double_t path_time;
Double_t Beta;
Double_t Mass;
*/
Double_t path_length[12];
Double_t path_time[12];
Double_t Beta[12];
Double_t Mass[12];
TH1D* h1E_TOF[12];
//Conversion
 Long_t histos_tofCao1();
 Long_t startup_tofCao1();
 Long_t process_tofCao1();
 Long_t done_tofCao1();


  virtual Long_t cmdline(char * cmd);

  
  
  ClassDef(Det_TOF,1);
};

#endif
