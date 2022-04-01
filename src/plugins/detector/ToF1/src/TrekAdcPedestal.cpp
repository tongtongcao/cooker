#include "TrekAdcPedestal.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TSpectrum.h"
#include <algorithm>
void TrekAdcPedestal::addData(double value){
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
void TrekAdcPedestal::setName(std::string name){
  mNameAdc=name;
}
double TrekAdcPedestal::getPedestal(bool doDraw){
  TF1* f1=mH1->GetFunction("gaus");
  if(f1==0){
    TSpectrum *s = new TSpectrum(4);
    Int_t nfound = s->Search(mH1,2,"",0.10);
    double lowRange,upRange;
    Float_t *xpeaks = s->GetPositionX();
    std::sort(xpeaks,xpeaks+nfound);
    lowRange=xpeaks[0]/2.0;
    upRange=xpeaks[0]*3.0/2.0;
    mH1->Fit("gaus","Q","",lowRange,upRange);

    mPedestal=mH1->GetFunction("gaus")->GetParameter(1);
    mError=mH1->GetFunction("gaus")->GetParError(1);
  }
  if(doDraw){
    mH1->Draw();
  }
  return mPedestal;
}
double TrekAdcPedestal::getError(){
  return mError;

}
