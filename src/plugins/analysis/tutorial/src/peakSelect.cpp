#include "peakSelect.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TSpectrum.h"
#include <algorithm>

Int_t npeaks = 30;
void peakSelect::addData(double value){
  if(mH1==0){
    std::string name(mNameAdc);
    name+="Adc";
    while(gDirectory->GetListOfKeys()->Contains(name.c_str())){
      name+="_";
    }
    mH1=new TH1D(name.c_str(),"AC ADC Peak spectrum;ADC [channel]",mNBin,mLowBound,mUpBound);
  }
  mH1->Fill(value);
}
//histo name 
void peakSelect::setName(std::string name){
  mNameAdc=name;
}
//how to find peaks in histograms by generating a random number of gaussian peaks on top of a linear background (mH1-> orignianl histo)
//The background is computed and drawn on top of the original histogram.
double peakSelect::getPedestal(bool doDraw){
TF1* f1=mH1->GetFunction("gaus");
  if(f1==0){
    TSpectrum *s = new TSpectrum(4);
    Int_t nfound = s->Search(mH1,1,"",0.10);
    double lowRange,upRange;
    if(nfound == 1){
      Float_t *xpeaks = s->GetPositionX();
      std::sort(xpeaks,xpeaks+nfound);
      lowRange=mH1->GetBinLowEdge(mH1->GetMaximumBin()-10);  //xpeaks[0]/2.0;
      upRange=mH1->GetBinLowEdge(mH1->GetMaximumBin()+12);   //xpeaks[0]*3.0/2.0;
      mH1->Fit("gaus","Q","",lowRange,upRange);
 
      mPedestal=mH1->GetFunction("gaus")->GetParameter(1);
      mError=mH1->GetFunction("gaus")->GetParError(1);
      mSigma=mH1->GetFunction("gaus")->GetParameter(2);
      mSigmaError=mH1->GetFunction("gaus")->GetParError(2);
    }else if (nfound == 2){
      Float_t *xpeaks = s->GetPositionX();
      std::sort(xpeaks,xpeaks+nfound);
      lowRange=mH1->GetBinLowEdge(mH1->GetMaximumBin()-11);
      upRange=mH1->GetBinLowEdge(mH1->GetMaximumBin()+11);   //xpeaks[0]*3.0/2.0;
      mH1->Fit("gaus","Q","",lowRange,upRange);
 
      mPedestal=mH1->GetFunction("gaus")->GetParameter(1);
      mError=mH1->GetFunction("gaus")->GetParError(1);
      mSigma=mH1->GetFunction("gaus")->GetParameter(2);
      mSigmaError=mH1->GetFunction("gaus")->GetParError(2);
    }else if (nfound > 2){
      Float_t *xpeaks = s->GetPositionX();
      std::sort(xpeaks,xpeaks+nfound);
      lowRange=mH1->GetBinLowEdge(mH1->GetMaximumBin()-25);  //xpeaks[0]/2.0;
      upRange=mH1->GetBinLowEdge(mH1->GetMaximumBin()+23);   //xpeaks[0]*3.0/2.0;
      mH1->Fit("gaus","Q","",lowRange,upRange);
 
      mPedestal=mH1->GetFunction("gaus")->GetParameter(1);
      mError=mH1->GetFunction("gaus")->GetParError(1);
      mSigma=mH1->GetFunction("gaus")->GetParameter(2);
      mSigmaError=mH1->GetFunction("gaus")->GetParError(2);
    }
}
  if(doDraw){
    mH1->Draw();
  }
  return mPedestal;
}

double peakSelect::getPosX(){
  return xpos;
}

double peakSelect::getError(){
  return mError;
}

double peakSelect::getSigma(){
  return mSigma;
}

double peakSelect::getSigmaError(){
  return mSigmaError;
}
