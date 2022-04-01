#include <Det_TOF.h>

#include<cmath>

#include <fstream>
#include <iostream>
#include <TStyle.h>
using namespace std;
Long_t Det_TOF::histos_tofCao1(){
  return 0;

}

static Int_t run,event;
static vector<Int_t> gapNumToF1, gapNumToF2;
static vector<Double_t> adcUToF1, adcDToF1, timeUToF1, timeDToF1;
static vector<Double_t> adcAIToF2, adcAOToF2, adcBIToF2, adcBOToF2, timeAIToF2, timeAOToF2, timeBIToF2, timeBOToF2; 
Long_t Det_TOF::startup_tofCao1(){
  getBranchObject("hitTof1",(TObject **) &treeTof1);
  getBranchObject("hitTof2",(TObject **) &treeTof2);

  out->Branch("run",&run,"run/I");
  out->Branch("event",&event,"event/I");
  out->Branch("gapNumToF1",&gapNumToF1);
  out->Branch("adcUToF1",&adcUToF1);
  out->Branch("adcDToF1",&adcDToF1);
  out->Branch("timeUToF1",&timeUToF1);
  out->Branch("timeDToF1",&timeDToF1);

  out->Branch("gapNumToF2",&gapNumToF2);
  out->Branch("adcAIToF2",&adcAIToF2);
  out->Branch("adcAOToF2",&adcAOToF2);
  out->Branch("adcBIToF2",&adcBIToF2);
  out->Branch("adcBOToF2",&adcBOToF2);
  out->Branch("timeAIToF2",&timeAIToF2);
  out->Branch("timeAOToF2",&timeAOToF2);
  out->Branch("timeBIToF2",&timeBIToF2);
  out->Branch("timeBOToF2",&timeBOToF2);

  gStyle->SetOptStat(0);

  return 0;
};


Long_t Det_TOF::process_tofCao1(){
  run=treeTof1->run;
  event=treeTof1->event;

  gapNumToF1.clear();
  adcUToF1.clear();
  adcDToF1.clear();
  timeUToF1.clear();
  timeDToF1.clear();

  gapNumToF2.clear();
  adcAIToF2.clear();
  adcAOToF2.clear();
  adcBIToF2.clear();
  adcBOToF2.clear();
  timeAIToF2.clear();
  timeAOToF2.clear();
  timeBIToF2.clear();
  timeBOToF2.clear();

//cout<<"event: "<<event<<endl;
//cout<<"ToF1:"<<endl;
  for(UInt_t i=0;i<treeTof1->gap.size();i++){
    gapNumToF1.push_back(treeTof1->gap[i]);
    adcUToF1.push_back(treeTof1->adcU[i]);
    adcDToF1.push_back(treeTof1->adcD[i]);
    timeUToF1.push_back(treeTof1->timeU[i]);
    timeDToF1.push_back(treeTof1->timeD[i]);
//cout<<treeTof1->gap[i]<<"  ";
  }
//cout<<endl;

//cout<<"ToF2:"<<endl;
  for(UInt_t i=0;i<treeTof2->gap.size();i++){
    gapNumToF2.push_back(treeTof2->gap[i]);
    adcAIToF2.push_back(treeTof2->adcAI[i]);
    adcAOToF2.push_back(treeTof2->adcAO[i]);
    adcBIToF2.push_back(treeTof2->adcBI[i]);
    adcBOToF2.push_back(treeTof2->adcBO[i]);

    timeAIToF2.push_back(treeTof2->timeAI[i]);
    timeAOToF2.push_back(treeTof2->timeAO[i]);
    timeBIToF2.push_back(treeTof2->timeBI[i]);
    timeBOToF2.push_back(treeTof2->timeBO[i]);
//cout<<treeTof2->gap[i]<<"  ";
  }
//cout<<endl;
 return 0;
};


Long_t Det_TOF::done_tofCao1(){

  return 0;
};



