#ifndef __MWPC_ADC_SPECTRA__
#define __MWPC_ADC_SPECTRA__
#include <string>
#include <TH1D.h>
#include <TDirectory.h>
#include <TList.h>
class MwpcAdcSpectra{
 private:
  unsigned int mNBin;
  double mLowBound;
  double mUpBound;
  std::string mNameAdc;
  TH1D* mH1;
 public:
 MwpcAdcSpectra():mNBin(100),mLowBound(0),mUpBound(1000),mNameAdc("unNamed"),mH1(0){}
 MwpcAdcSpectra(unsigned int nBin,double low,double up,std::string name):mNBin(nBin),mLowBound(low),mUpBound(up),mNameAdc(name),mH1(0){}
  void addData(double value);
  void setName(std::string name);
  void savePlot();

  virtual ~MwpcAdcSpectra();
  ClassDef(MwpcAdcSpectra,1);
};

#endif
