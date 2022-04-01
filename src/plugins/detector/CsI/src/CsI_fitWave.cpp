#include <Det_CsI.h>
#include <singleCsI.h>
#include<cmath>
#include <TMath.h>
#include <fstream>
#include <iostream>
#include <sstream>
#include "TF1.h"
#include <cstdlib>
#include <TROOT.h>
#include <TCanvas.h>
#include <TLine.h>
#include <TStyle.h>
#include <time.h>
using namespace std;
Long_t Det_CsI::set_goodEvents(int run, int event)
{
  cout<<"set good events..."<<endl;
  return 0;
};

Long_t Det_CsI::histos_fit(){
  h2TdcVSCsI=dH2("TdcVSCsI","fitted Tdc;CsI index;tdc [tdc channel]",768,0.5,768.5,1001,-0.5,1000.5);
  h2AdcVSCsI=dH2("AdcVSCsI","Integrated Adc;CsI index;adc [adc channel]",768,0.5,768.5,200,0.0,2000.0);
  h1Skipped=dH1("SkippedCsI","Skipped CsI;CsI index",768,0.5,768.5);
  h2BaseVSCsI=dH2("BaseVSCsI","base line of CsI",768,0.5,768.5,100,0.0,1000.0);
  return 0;
}
Long_t Det_CsI::startup_fit(){
  getBranchObject("vf48",(TObject **) &treeRaw);
  treeCali=new CRTCaliCsI();
  // Make new branch on output tree                                                                                                                                                             
  makeBranch("csi_cali",(TObject **)&treeCali);                                                                                                                  gROOT->SetBatch(kTRUE);
  gStyle->SetOptStat(0);
  return 0;
}
Double_t firstDerive(Double_t *x,Double_t *par){


}
Double_t waveform(Double_t *x,Double_t *par){
  //  if(x[0]<1 && x[0]>250) return 0;
  Double_t termFirst=par[0]/(1+TMath::Exp(-par[2]*(x[0]-par[7])));
  Double_t termSecond=(x[0]-par[1])/(par[3]*par[3]);
  Double_t termThird=TMath::Exp(-(x[0]-par[1])/par[3])+par[4]*TMath::Exp(-(x[0]-par[1])/par[5]);
  Double_t value=termFirst*termSecond*termThird+par[6];
  if(value>1023.0) value=1023.0;
  return value;

}

Long_t Det_CsI::process_fit(){
  static unsigned int nEvent=0;
  if(nEvent%100==0){
    time_t t=time(0);
    struct tm* now=localtime(&t);
    std::cout<<nEvent<<" events processed at "<<asctime(now)<<endl;
  }
  
  treeCali->runNo=treeRaw->runNo;
  treeCali->eventNo=treeRaw->eventNo;
  treeCali->isBad=treeRaw->isBad;
  treeCali->indexCsI.clear();
  treeCali->tdc.clear();
  treeCali->adc.clear();
  treeCali->peak.clear();
  treeCali->baseline.clear();
  //process all events even slipped events
  /*
    if(treeRaw->isBad<0){
    cout<<"slipped event"<<endl;
    return 0;
    }
  */
  treeCali->nChannel=0;//counter
  for(UInt_t iCh=0;iCh<treeRaw->nChannel;iCh++){
    //    cout<<"processing "<<iCh+1<<"/"<<treeRaw->nChannel<<endl;
    UInt_t myEvent=treeRaw->eventNo;
    std::stringstream ss;
    IdCsI myIndex(treeRaw->nameCsI[iCh],treeRaw->indexCsI[iCh]);
    UInt_t iCsI=mapCsI[myIndex];
    if(iCsI==16||iCsI==272||iCsI==528||iCsI==144 || iCsI==400 || iCsI==656){
      h1Skipped->Fill(iCsI);
      continue;
    }
    SingleCsI myCsI(treeRaw->runNo,myEvent,iCsI);
    myCsI.setData(treeRaw->data[iCh]);
    if(nEvent<10)
      myCsI.quick_fit(true);
    else
      myCsI.quick_fit();
    for(UInt_t iWave=0,nWave=myCsI.numberWave();iWave<nWave;iWave++){
      h2TdcVSCsI->Fill(iCsI,myCsI.tdc(iWave));
      h2AdcVSCsI->Fill(iCsI,myCsI.adc(iWave));
      h2BaseVSCsI->Fill(iCsI,myCsI.baseline());
      treeCali->indexCsI.push_back(iCsI);
      treeCali->tdc.push_back(myCsI.tdc(iWave));
      treeCali->adc.push_back(myCsI.adc(iWave));
      treeCali->peak.push_back(myCsI.peak(iWave));
      treeCali->baseline.push_back(myCsI.baseline());
      treeCali->nChannel++;
    }      
  }
  nEvent++;
  return 0;
}
Long_t Det_CsI::finalize_fit(){
  cout<<"finalize_fit"<<endl;
  return 0;
}
