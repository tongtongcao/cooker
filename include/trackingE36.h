#ifndef trackingE36_h
#define trackingE36_h 1
#include "cookerrawtree.h" // for CRTBase

#include <TObject.h>
#include <vector>
//const int maxEvt=151212;
/*class CRTBase: public TObject
{
 public:
  CRTBase();
  virtual ~CRTBase();
  ClassDef(CRTBase,1);
};
*/
class trackingE36 : public CRTBase {
public:
  Int_t evtNum, fgapNumTof2, nTracks;
  Double_t fTof2SP,fChi2PerNDF;
  Double_t fTof2SX,fTof2SY,fTof2SZ;
  Double_t fTof2SNx,fTof2SNy,fTof2SNz;
Double_t fTof1PathLocation;
Double_t fTof2PathLocation;
Double_t fTof1PathTime;
Double_t fTof2PathTime;
Double_t fTof1SP;
Int_t fFlagTargetToAC;
Bool_t fFlagToFMatch;


  trackingE36();
  ~trackingE36();
  ClassDef(trackingE36, 2);
};
#endif
