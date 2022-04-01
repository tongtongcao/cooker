#include "MwpcAdcSpectra.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TLine.h"
#include "TSpectrum.h"
#include <algorithm>
#include <iostream>
#include <sstream>
using namespace std;

void MwpcAdcSpectra::addData(double value){
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
void MwpcAdcSpectra::setName(std::string name){
  mNameAdc=name;
}
void MwpcAdcSpectra::savePlot(){
  ostringstream filename;
  if(mH1==0){
    cout<<"Histogram is not ready."<<endl;
    exit;
  }
  else{
    TCanvas *canv=new TCanvas(mNameAdc.c_str(),mNameAdc.c_str(),500,300);
    gPad->SetLogy();
    mH1->Draw();
    filename<<mNameAdc<<".png";
    canv->SaveAs(filename.str().c_str());
    filename.str("");
  }
}

MwpcAdcSpectra::~MwpcAdcSpectra(){}
ClassImp(MwpcAdcSpectra);
