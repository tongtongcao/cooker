#ifndef anaE36_h
#define anaE36_h 1
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
class anaE36 : public CRTBase {
public:
  Int_t evtNum, fgapNumTof2,nTracks;
  Double_t fTof2SP;
  Double_t fTof2SX,fTof2SY,fTof2SZ;
  Double_t fTof2SNx,fTof2SNy,fTof2SNz;
  Double_t fTof2PathLocation;
  Double_t fTof2PathTime;
  
  Int_t fgapNumTof1;
  Double_t fTof1SP,fVertSP;
  Double_t fTof1SX,fTof1SY,fTof1SZ;
  Double_t fTof1SNx,fTof1SNy,fTof1SNz;
  Double_t fTof1PathLocation;
  Double_t fTof1PathTime;
  Double_t theta, phi;
  Double_t p_zx;
  Double_t fMwpcTotalEloss; 
  Double_t fMwpc4SP; 
  Double_t fChi2PerNDF;
  Int_t fFlagTargetToAC;
  Bool_t fFlagToFMatch;

  anaE36();
  ~anaE36();
  ClassDef(anaE36, 4);
};
#endif
