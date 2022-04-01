#ifndef __TARGET__
#define __TARGET__
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
//#include "TrekAdcPedestal.h"


class Target:public Plugin
{
   private:
  
	//CRTCaliTof1 *treeCali;
  	TargetInfo *treeRaw;			/// Input tree with TARGET raw data
 	BeamInfo* treeBeam;


   public:
      Target(TTree *in, TTree *out,TFile *inf_, TFile * outf_, TObject *p);
      virtual ~Target();
      // add funtions with return value Long_t here:

      // Main ToF1 analysis methods
      Long_t histos();
      Long_t startup();
      Long_t process();
      Long_t done();
      //Long_t finalize();
 
      TH1D* hADC_HG[256];
      TH1D* hADC_LG[256];
      TH1D* hTDC_LE[256];
      TH1D* hTDC_TE[256];
 
      virtual Long_t cmdline(char * cmd);

      ClassDef(Target,1);
};

#endif
