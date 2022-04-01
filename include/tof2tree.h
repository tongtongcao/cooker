/** 
 * This file has the definition of the ToF tree branch objects 
 */

#ifndef __TOF2TREE_H_
#define __TOF2TREE_H_
#include "cookerrawtree.h" // for CRTBase
#include <vector>
class CRTCaliTof2:public CRTBase
{
 public:
  Int_t run;
  Int_t event;
  std::vector<Int_t> gap;//from 1 to 12
  std::vector<Float_t> adcAI;//ped-subtracted
  std::vector<Float_t> adcAO;
  std::vector<Float_t> adcBI;//ped-subtracted
  std::vector<Float_t> adcBO;

  std::vector<Float_t> timeAI;
  std::vector<Float_t> timeAO;
  std::vector<Float_t> timeBI;
  std::vector<Float_t> timeBO;

  CRTCaliTof2();
  virtual ~CRTCaliTof2();
  ClassDef(CRTCaliTof2,1);
};

#endif
