#ifndef PEAKSELECT_H
#define PEAKSELECT_H 1

#include <string>
#include <TH1D.h>
#include <TDirectory.h>
#include <TSpectrum.h>
#include <TList.h>

class peakSelect{
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
  peakSelect():mNBin(100),mLowBound(0),mUpBound(1000),mNameAdc("unNamed"),mH1(0),mPedestal(0),mError(0){}
 peakSelect(unsigned int nBin,double low,double up,std::string name):mNBin(nBin),mLowBound(low),mUpBound(up),mNameAdc(name),mH1(0),mPedestal(0),mError(0){}
  void addData(double value);
  void setName(std::string name);
  double getPedestal(bool writeImg=false);
  double getPosX();
  double getError();
  double getSigma();
  double getSigmaError();
  double getFWHM();
  TH1D* h1Adc_PgC[12][7];



};

#endif
