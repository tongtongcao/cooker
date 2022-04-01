#include "TrekTimewalk.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TSpectrum.h"
#include <algorithm>
void TrekTimewalk::addData(const double& diff,const double& adc){
  if(mH2==0){
    std::string name(mName);
    //name+="Adc";
    while(gDirectory->GetListOfKeys()->Contains(name.c_str())){
      name+="_";
    }
    mH2=new TH2D(name.c_str(),"Diff VS Adc;Adc;Diff (ns)",mNBinsX,mLowBoundX,mUpBoundX,mNBinsY,mLowBoundY,mUpBoundY);
  }
  mH2->Fill(adc,diff);
}
void TrekTimewalk::setName(std::string name){
  mName=name;
}
bool TrekTimewalk::getValue(double& lambda,double& order,bool doDraw){
  if(mProf==0){
    mProf=mH2->ProfileX();
  }
  if(doDraw){
    mH2->Draw("colz");
    mProf->Draw("same");
  }
  return true;
}
void TrekTimewalk::getError(double& errLambda,double& errOrder){
  errLambda=mErrorConst;
  errOrder=mErrorOrder;

}
