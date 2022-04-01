/*
 *
 * This file contains the definition of AC tree branch obj 
 * Author: Bishoy Dongwi
 *
 */

#ifndef __ACTREE__H
#define __ACTREE__H
#include "cookerrawtree.h" // for CRTBase
#include <vector>
class CRTCaliAc : public CRTBase {
public:
  Int_t run;
  Int_t event;
//Double_t pVertPosition[50];  
//Double_t Ped_ACU;
  //Double_t Ped_ACD;
  std::vector<Float_t>adcACU;
  std::vector<Float_t>adcACD;

  
//  std::vector<Float_t>cer_adc_sum;
  Double_t TO_ACU[12];
  Double_t TO_ACD[12];

  CRTCaliAc();
  virtual ~CRTCaliAc();
  ClassDef(CRTCaliAc, 1);
};



#endif
