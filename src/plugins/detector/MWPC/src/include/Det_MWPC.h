#ifndef __DET_MWPC__
#define __DET_MWPC__
#include "TObject.h"
#include "Plugin.h"
#include "TTree.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TF1.h"
#include "mwpctree.h"
#include <iostream>
#include <iomanip>
#include <vector>
#include <cmath>
#include <utility>
#include <map>
#include "MwpcAdcSpectra.h"
#include "MwpcPedestal.h"
#include "MwpcEvent.h"
#include "MwpcPosition.h"

class Det_MWPC:public Plugin
{
 private:

  CRTCaliMwpc *treeCali;                // Tree for calibration analysis
  CRTCookMwpcV1 *treeCookV1;            // Tree to save cooked data of version 1
  CRTCookMwpcV2 *treeCookV2;            // Tree to save cooked data of version 2
  MwpcInfo *treeRaw;			// Input tree with MWPC raw data
  BeamInfo* treeBeam;

  // Detector parameters set in init file

  // constants 

 public:
  Det_MWPC(TTree *in_,TTree *out_,TFile *inf_, TFile * outf_,TObject *p_);
  virtual ~Det_MWPC();

  //Preliminary analysis for ADC readout: Plot ADC vs. Strip, ADC Spectra, and Spectra of two dead strip
  Long_t histos();
  Long_t startup();
  Long_t process();
  Long_t done();

  //Pedestal Determination: Plot ADC spectra of pedestal runs, and then get pedestal and its sigma by guassian fit 
  Long_t histos_ped();
  Long_t startup_ped();
  Long_t process_ped();
  Long_t done_ped();


  //Event Analysis: Plot hit distribution, multiplet distribution; Divide all events into 8 types, and plot their occupation distribution
  Long_t histos_event();
  Long_t startup_event();
  Long_t process_event();
  Long_t done_event();

  //Cluster Analysis: Plot cluster distribution
  Long_t histos_cluster();
  Long_t startup_cluster();
  Long_t process_cluster();
  Long_t done_cluster();

  //Threshold analysis
  Long_t histos_eventthresholdFirstApproach();
  Long_t startup_eventthresholdFirstApproach();
  Long_t process_eventthresholdFirstApproach();
  Long_t done_eventthresholdFirstApproach();

  Long_t histos_eventthresholdSecondApproach();
  Long_t startup_eventthresholdSecondApproach();
  Long_t process_eventthresholdSecondApproach();
  Long_t done_eventthresholdSecondApproach();

  Long_t histos_clusterthresholdFirstApproach();
  Long_t startup_clusterthresholdFirstApproach();
  Long_t process_clusterthresholdFirstApproach();
  Long_t done_clusterthresholdFirstApproach();

  Long_t histos_clusterthresholdSecondApproach();
  Long_t startup_clusterthresholdSecondApproach();
  Long_t process_clusterthresholdSecondApproach();
  Long_t done_clusterthresholdSecondApproach();

  Long_t histos_xclusterdoub();
  Long_t startup_xclusterdoub();
  Long_t process_xclusterdoub();
  Long_t done_xclusterdoub();

  //Determine total-charege cut value 
  Long_t histos_clusterthresholdfromcharge();
  Long_t startup_clusterthresholdfromcharge();
  Long_t process_clusterthresholdfromcharge();
  Long_t done_clusterthresholdfromcharge();

  //Compare XY 
  Long_t histos_comparexy();
  Long_t startup_comparexy();
  Long_t process_comparexy();
  Long_t done_comparexy();

  //Extract clusters
  Long_t histos_extractCluster();
  Long_t startup_extractCluster();
  Long_t process_extractCluster();
  Long_t done_extractCluster();
      
  Long_t histos_analysisExtractedCluster();
  Long_t startup_analysisExtractedCluster();
  Long_t process_analysisExtractedCluster();
  Long_t done_analysisExtractedCluster();

  Long_t histos_analysisExtractedCluster2();
  Long_t startup_analysisExtractedCluster2();
  Long_t process_analysisExtractedCluster2();
  Long_t done_analysisExtractedCluster2();

  Long_t histos_analysisExtractedCluster3();
  Long_t startup_analysisExtractedCluster3();
  Long_t process_analysisExtractedCluster3();
  Long_t done_analysisExtractedCluster3();

  // Analysis noise
  Long_t histos_analysisNoise();
  Long_t startup_analysisNoise();
  Long_t process_analysisNoise();
  Long_t done_analysisNoise();

  //Check cable swapping
  Long_t histos_checkCableSwapping();
  Long_t startup_checkCableSwapping();
  Long_t process_checkCableSwapping();
  Long_t done_checkCableSwapping();

  //Gain Calibration
  Long_t histos_amplifierGainCalibration();
  Long_t startup_amplifierGainCalibration();
  Long_t process_amplifierGainCalibration();
  Long_t done_amplifierGainCalibration();

  Long_t histos_testGain();
  Long_t startup_testGain();
  Long_t process_testGain();
  Long_t done_testGain();

  //Handle dead strip
  Long_t histos_handleDeadStrip();
  Long_t startup_handleDeadStrip();
  Long_t process_handleDeadStrip();
  Long_t done_handleDeadStrip();
 
  //Get Position
  Long_t histos_demonstratePosition();
  Long_t startup_demonstratePosition();
  Long_t process_demonstratePosition();
  Long_t done_demonstratePosition();

  Long_t histos_posPuzzle();
  Long_t startup_posPuzzle();
  Long_t process_posPuzzle();
  Long_t done_posPuzzle();

  Long_t histos_test();
  Long_t startup_test();
  Long_t process_test();
  Long_t done_test();

  Long_t histos_compGlobalPosition();
  Long_t startup_compGlobalPosition();
  Long_t process_compGlobalPosition();
  Long_t done_compGlobalPosition();

  // Cooked Data
  Long_t histos_mwpcCookedDataV1();
  Long_t startup_mwpcCookedDataV1();
  Long_t process_mwpcCookedDataV1();
  Long_t done_mwpcCookedDataV1();

  Long_t histos_mwpcCookedDataV11();
  Long_t startup_mwpcCookedDataV11();
  Long_t process_mwpcCookedDataV11();
  Long_t done_mwpcCookedDataV11();

  Long_t histos_mwpcCookedDataV2();
  Long_t startup_mwpcCookedDataV2();
  Long_t process_mwpcCookedDataV2();
  Long_t done_mwpcCookedDataV2();

  Long_t histos_mwpcCookedDataV21();
  Long_t startup_mwpcCookedDataV21();
  Long_t process_mwpcCookedDataV21();
  Long_t done_mwpcCookedDataV21();

  // Position Resolution
  Long_t histos_posResolution();
  Long_t startup_posResolution();
  Long_t process_posResolution();
  Long_t done_posResolution();

  // Track
  Long_t histos_trackV1();
  Long_t startup_trackV1();
  Long_t process_trackV1();
  Long_t done_trackV1();

  Long_t histos_trackV2();
  Long_t startup_trackV2();
  Long_t process_trackV2();
  Long_t done_trackV2();

  virtual Long_t cmdline(char * cmd);

  ClassDef(Det_MWPC,1);
};

#endif
