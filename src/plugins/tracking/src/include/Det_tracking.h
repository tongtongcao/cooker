#ifndef __DET_TRACKING__
#define __DET_TRACKING__

#include "TObject.h"
#include "TTree.h"
#include "TClonesArray.h"

#include "Plugin.h"
#include "track.h"
#include "mergedCaliTargetTofMwpc.h"


#include "utility.h"
#include "matrix.h"
#include "track.h"
#include "fieldMap.h"
#include "fieldStepper.h"
#include "trackState.h"
#include "trackSite.h"
#include "trackSystem.h"
#include "trackFinder.h"
#include "particle.h"
#include "planeHit.h"
#include "plane.h"
#include "target.h"
#include "sft.h"
#include "ac.h"
#include "tof1.h"
#include "mwpc.h"
#include "tof2.h"
#include "actree.h"

#include <iostream>
using namespace std;
class Det_tracking:public Plugin{
  private:
    CRTCaliTargetTofMwpc *caliTargetTofMwpc;
    trackArray *trackArr;
    AcInfo *treeRaw;			/// Input tree with AC raw data

    sft *uSft;
    vector<double> sftX, sftY, sftZ;
    vector<double> sftXErr, sftYErr, sftZErr;
    ac *uAc;

    particleId partId;

  public:
    Det_tracking(TTree *in, TTree *out,TFile *inf_, TFile * outf_, TObject *p);
    virtual ~Det_tracking();
      // add funtions with return value Long_t here:
    Long_t histos();
    Long_t startup();
    Long_t process();
    Long_t done();

    Long_t histos_basicDataType();
    Long_t startup_basicDataType();
    Long_t process_basicDataType();
    Long_t done_basicDataType();


    virtual Long_t cmdline(char * cmd);

    ClassDef(Det_tracking,1);
};

#endif
