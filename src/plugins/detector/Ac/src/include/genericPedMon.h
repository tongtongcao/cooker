#ifndef __GENERICPEDMON__
#define __GENERICPEDMON__
#include <string>
#include <TH1D.h>
#include "TF1.h"
#include <TStyle.h>
#include <TSpectrum.h>
#include <TDirectory.h>
#include <TList.h>
class genericPedMon{
 private:
  unsigned int mNBin;
  double mLowBound;
  double mUpBound;
  std::string mNameAdc;
  TH1D* mH1;
  double mPedestal;
  double xpos;
  double mError;
  double mSigma;
  double mSigmaError;
  double fwhm;
 public:
 genericPedMon():mNBin(100),mLowBound(0),mUpBound(1000),mNameAdc("unNamed"),mH1(0),mPedestal(0),mError(0){}
 genericPedMon(unsigned int nBin,double low,double up,std::string name):mNBin(nBin),mLowBound(low),mUpBound(up),mNameAdc(name),mH1(0),mPedestal(0),mError(0){}
  ~genericPedMon();
  void addData(double value);
  void setName(std::string name);
  double getPedestal(TSpectrum* s, bool writeImg=false);
  double getPosX();
  double getError();
  double getSigma();
  double getFWHM();
  double getSigmaError();
 TH1D* h1adc_ACu[12];
  TH1D* h1adc_ACd[12];

};

#endif
