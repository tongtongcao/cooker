/*
 *
 * This file contains the definition of AC tree branch obj 
 * Author: Bishoy Dongwi
 *
 */

#ifndef __PGCTREE__H
#define __PGCTREE__H

#include "cookerrawtree.h" // for CRTBase
#include <vector>
class CRTCaliPgC : public CRTBase {
public:
  Int_t run;
  Int_t event;
 // Double_t Ped_PGC;
//  std::vector<Float_t>ADC_PGC_Gap;
  Double_t TO_PGC[12];
  Double_t adcPgC[12][7];
//std::vector<Float_t>adcPgC;
//
  CRTCaliPgC();
  virtual ~CRTCaliPgC();
  ClassDef(CRTCaliPgC, 1);
};
#endif

