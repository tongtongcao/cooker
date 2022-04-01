#include "genericPedMon.h"
#include "TCanvas.h"
#include <algorithm>
#include <cmath>

using namespace std;

void genericPedMon::addData(double value){
  if(mH1==0){
    std::string name(mNameAdc);
    name+="Adc";
    while(gDirectory->GetListOfKeys()->Contains(name.c_str())){
      name+="_";
    }
    mH1=new TH1D(name.c_str(),"ADC spectrum;ADC [channel]",mNBin,mLowBound,mUpBound);
  }
  mH1->Fill(value);
}
void genericPedMon::setName(std::string name){
  mNameAdc=name;
}
genericPedMon::~genericPedMon(){
  delete mH1;
}
double genericPedMon::getPedestal(TSpectrum* s, bool doDraw){
  TF1* f1=mH1->GetFunction("gaus");
  if(f1==0){
    Int_t nfound = s->Search(mH1,1,"",0.10);
    //std::cout << " --- Check number of peaks found: " << nfound << std::endl;
    double lowRange, upRange;
    if(nfound == 1){
      Float_t *xpeaks = s->GetPositionX();
      std::sort(xpeaks,xpeaks+nfound);
      lowRange=mH1->GetBinLowEdge(mH1->GetMaximumBin()-37);  //xpeaks[0]/2.0;
      upRange=mH1->GetBinLowEdge(mH1->GetMaximumBin()+10);   //xpeaks[0]*3.0/2.0;
      mH1->Fit("gaus","Q","",lowRange,upRange);
 
      mPedestal=mH1->GetFunction("gaus")->GetParameter(1);
      mError=mH1->GetFunction("gaus")->GetParError(1);
      mSigma=mH1->GetFunction("gaus")->GetParameter(2);
      mSigmaError=mH1->GetFunction("gaus")->GetParError(2);
      fwhm=2*sqrt(2*log(2))*mSigma;
    }else if (nfound == 2){
      Float_t *xpeaks = s->GetPositionX();
      std::sort(xpeaks,xpeaks+nfound);
      lowRange=mH1->GetBinLowEdge(mH1->GetMaximumBin()-xpeaks[0]/2.0);
      upRange=mH1->GetBinLowEdge(mH1->GetMaximumBin()+11);   //xpeaks[0]*3.0/2.0;
      mH1->Fit("gaus","Q","",lowRange,upRange);
 
      mPedestal=mH1->GetFunction("gaus")->GetParameter(1);
      mError=mH1->GetFunction("gaus")->GetParError(1);
      mSigma=mH1->GetFunction("gaus")->GetParameter(2);
      mSigmaError=mH1->GetFunction("gaus")->GetParError(2);
      fwhm=2*sqrt(2*log(2))*mSigma;
    }else if (nfound > 2){
      Float_t *xpeaks = s->GetPositionX();
      std::sort(xpeaks,xpeaks+nfound);
      lowRange=mH1->GetBinLowEdge(mH1->GetMaximumBin()-21);  //xpeaks[0]/2.0;
      upRange=mH1->GetBinLowEdge(mH1->GetMaximumBin()+32);   //xpeaks[0]*3.0/2.0;
      mH1->Fit("gaus","Q","",lowRange,upRange);
 
      mPedestal=mH1->GetFunction("gaus")->GetParameter(1);
      mError=mH1->GetFunction("gaus")->GetParError(1);
      mSigma=mH1->GetFunction("gaus")->GetParameter(2);
      mSigmaError=mH1->GetFunction("gaus")->GetParError(2);
      fwhm=2*sqrt(2*log(2))*mSigma;
    }
  }
  if(doDraw){
    mH1->Draw();
  }
  return mPedestal;
}

double genericPedMon::getError(){
  return mError;
}

double genericPedMon::getSigma(){
  return mSigma;
}

double genericPedMon::getSigmaError(){
  return mSigmaError;
}

double genericPedMon::getFWHM(){
  return fwhm;
}

