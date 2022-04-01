#include "MwpcPedestal.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TLine.h"
#include "TSpectrum.h"
#include <algorithm>
#include <iostream>
#include <sstream>
using namespace std;

void MwpcPedestal::addData(double value){
  if(mH1==0){
    std::string name(mNameAdc);
    name+="";
    while(gDirectory->GetListOfKeys()->Contains(name.c_str())){
      name+="_";
    }
    ostringstream title;
    title<<name<<";ADC;Counts";
    mH1=new TH1D(name.c_str(),title.str().c_str(),mNBin,mLowBound,mUpBound);
    title.str("");
  }
  mH1->Fill(value);
}
void MwpcPedestal::setName(std::string name){
  mNameAdc=name;
}
double MwpcPedestal::getPedestal(bool doDraw){
  if(mH1==0){
    mPedestal=-100;
    mError=-100;
    return mPedestal;
  }
  else{
    TF1* f1=mH1->GetFunction("gaus");
    if(f1==0){
      double lowRange,upRange;
      lowRange=mH1->GetBinLowEdge(mH1->GetMaximumBin()-10);
      upRange=mH1->GetBinLowEdge(mH1->GetMaximumBin()+11);
      mH1->Fit("gaus","Q","",lowRange,upRange);

      mPedestal=mH1->GetFunction("gaus")->GetParameter(1);
      mError=mH1->GetFunction("gaus")->GetParError(1);
      mSigma=mH1->GetFunction("gaus")->GetParameter(2);
      mSigmaError=mH1->GetFunction("gaus")->GetParError(2);
    }
    if(doDraw){
      mH1->Draw();
    }
    return mPedestal;
  }
}
double MwpcPedestal::getError(){
  return mError;
}
double MwpcPedestal::getSigma(){
  return mSigma;
}
double MwpcPedestal::getSigmaError(){
  return mSigmaError;
} 
void MwpcPedestal::savePlot(){
  if(mH1==0){
    cout<<"Histogram is not ready."<<endl;
    exit;
  }
  else{
    TCanvas *canv=new TCanvas(mNameAdc.c_str(),mNameAdc.c_str(),500,300);
    //gPad->SetLogx();
    mH1->Draw();
    TLine *line=new TLine(mPedestal,0,mPedestal,mH1->GetMaximum());
    line->SetLineColor(kRed);
    //line->Draw();
    ostringstream filename;
    filename<<mNameAdc<<".pdf";
    //canv->SaveAs(filename.str().c_str());
    filename.str("");
    filename<<mNameAdc<<".png";
    canv->SaveAs(filename.str().c_str());
    filename.str("");
  }
}

MwpcPedestal::~MwpcPedestal(){}
ClassImp(MwpcPedestal);
