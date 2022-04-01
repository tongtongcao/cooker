#ifndef __MERGEDCALITARGETTOFMWPC_H_
#define __MERGEDCALITARGETTOFMWPC_H_
#include "cookerrawtree.h" // for CRTBase

class CRTCaliTargetTofMwpc:public CRTBase{
  public:
    UInt_t run,event;
    UInt_t nHitsC2X[12],nHitsC2Y[12],nHitsC2Z[12],nHitsC3X[12],nHitsC3Y[12],nHitsC3Z[12],nHitsC4X[12],nHitsC4Y[12],nHitsC4Z[12];
    Double_t posC2X[12][28], posC2Y[12][8], posC2Z[12][28], posC3X[12][32], posC3Y[12][8], posC3Z[12][1], posC4X[12][36], posC4Y[12][8], posC4Z[12][1];

    Int_t gapNumTof1;
    Int_t gapNumTof2;
    Double_t timeUTof1, timeDTof1;
    Double_t timeAITof2, timeAOTof2, timeBITof2, timeBOTof2;

    Int_t tag;
    Double_t vertX,vertY,vertPhi;
    Double_t vertXErr,vertYErr,vertPhiErr;

    CRTCaliTargetTofMwpc();
    virtual ~CRTCaliTargetTofMwpc();
    ClassDef(CRTCaliTargetTofMwpc,1);
};

#endif
