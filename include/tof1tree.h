/** 
 * This file has the definition of the ToF tree branch objects 
 */

#ifndef __TOF1TREE_H_
#define __TOF1TREE_H_
#include "cookerrawtree.h" // for CRTBase
#include <vector>
class CRTCaliTof1:public CRTBase
{
 public:
  Int_t run;
  Int_t event;
  std::vector<Int_t> gap;//from 1 to 12
  std::vector<Float_t> adcU;//ped-subtracted
  std::vector<Float_t>adcD;
  std::vector<Float_t>timeU;
  std::vector<Float_t>timeD;//converted from tdc channel to ns
  CRTCaliTof1();
  virtual ~CRTCaliTof1();
  ClassDef(CRTCaliTof1,2);
};

#endif
