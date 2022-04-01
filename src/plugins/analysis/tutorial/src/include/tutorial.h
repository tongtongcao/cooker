#ifndef __TUTORIAL__
#define __TUTORIAL__

#include "TObject.h"
#include "Plugin.h"
#include "TTree.h"
#include "actree.h"


#include "TH1D.h"
#include "TH2D.h"
#include "TF1.h"
#include "TCanvas.h"
#include <iostream>
#include <vector>
#include <cmath>
#include <iomanip>
#include <utility>
#include <map>
#include "peakSelect.h"

class CRTCaliAc;

class tutorial:public Plugin
{
	private:
		CRTCaliAc *treeAc;
		AcInfo *treeRaw;
		BeamInfo* treeBeam;


		double mPedestal;
		double mError;
		double mSigma;
		double mSigmaError;
	public:
		tutorial(TTree *in, TTree *out,TFile *inf_, TFile * outf_, TObject *p);
		virtual ~tutorial();
		// add funtions with return value Long_t here:


		TTree* tree;
		Long_t histos();
		Long_t startup();
		Long_t process();
		Long_t done();

		//Fit Pedestal
		Long_t start_calibS();
		Long_t hist_calibS(); 
		Long_t calib_peakS();
		Long_t done_calibS();	
		Long_t calib();


		peakSelect* pedADC_AcU[12];
		peakSelect* pedADC_AcD[12];
		TH1D* h1Pedestal;


		TH1D* h1adc_ACu[12];
		TH1D* h1adc_ACd[12];
                 double getFWHM();
		virtual Long_t cmdline(char * cmd);

		ClassDef(tutorial,1);
};

#endif
