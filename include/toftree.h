/** 
 * This file has the definition of the ToF tree branch objects 
 */

#ifndef __TOFTREE_H_
#define __TOFTREE_H_
#include "cookerrawtree.h" // for CRTBase
#include <vector>
class CRTCaliTof:public CRTBase
{
 public:
  Int_t run;
  Int_t event;
  std::vector<Int_t> gapToF1;//from 1 to 12
  std::vector<Int_t> gapToF2;//from 1 to 12
  std::vector<Float_t> adcUToF1;//ped-subtracted
  std::vector<Float_t>adcDToF1;
  std::vector<Float_t>timeUToF1;
  std::vector<Float_t>timeDToF1;//converted from tdc channel to ns
  std::vector<Float_t> adcUToF2;//ped-subtracted
  std::vector<Float_t>adcDToF2;
  std::vector<Float_t>timeUToF2;
  std::vector<Float_t>timeDToF2;//converted from tdc channel to ns
  CRTCaliTof();
  virtual ~CRTCaliTof();
  ClassDef(CRTCaliTof,2);
};

#endif
