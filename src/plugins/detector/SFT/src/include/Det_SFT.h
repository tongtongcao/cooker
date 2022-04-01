#ifndef __DET_SFT__
#define __DET_SFT__
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


class Det_SFT:public Plugin
{
   private:
	
	//CRTCaliTof1 *treeCali;
  	SftInfo *treeRaw;			/// Input tree with SFT raw data
 	BeamInfo* treeBeam;


   public:
      Det_SFT(TTree *in, TTree *out,TFile *inf_, TFile * outf_, TObject *p);
      virtual ~Det_SFT();
      // add funtions with return value Long_t here:

      // Main TARGET analysis methods
      Long_t histos();
      Long_t startup();
      Long_t process();
      Long_t done();
      //Long_t finalize();

      TH1D* hADC_HG_SFT[128];
      TH1D* hADC_LG_SFT[128];
      TH1D* hTDC_LE_SFT[128];
      TH1D* hTDC_TE_SFT[128];


      virtual Long_t cmdline(char * cmd);

      ClassDef(Det_SFT,1);
};

#endif
