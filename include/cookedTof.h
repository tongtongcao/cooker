#ifndef cookedTof_h
#define cookedTof_h 1

#include <TObject.h>
#include <vector>
#include "cookerrawtree.h" // for CRTBase

class cookedTof : public CRTBase {
public:
  Int_t evtNum;
  Int_t gapNumToF1;
    Int_t gapNumToF2;
    Double_t timeUToF1, timeDToF1;
    Double_t timeAIToF2, timeAOToF2, timeBIToF2, timeBOToF2;
  double adcUToF1, adcDToF1;
  double adcAIToF2, adcAOToF2;
  double adcBIToF2, adcBOToF2;    

  cookedTof();
  ~cookedTof();
  ClassDef(cookedTof, 2);
};
#endif
